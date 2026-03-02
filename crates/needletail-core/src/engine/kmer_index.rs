//! On-disk k-mer seed table: precomputed BWT [l, r] intervals for k-mers.
//!
//! Two storage strategies based on K:
//!
//! **Dense (K ≤ 12):** Flat array indexed by 2-bit rank. 4^K × 8 bytes.
//! O(1) lookup. Absent k-mers stored as `(u32::MAX, 0)`.
//!
//! **Sparse (K ≥ 13):** Only k-mers present in the genome are stored.
//! Sorted array of `(rank: u32, l: u32, r: u32)` triples, 12 bytes each.
//! O(log n) lookup via binary search. For K=14 on a 2 MB genome this is
//! ~24 MB instead of 2 GB — a ~85× reduction.
//!
//! At build time for sparse tables, we scan the genome text to collect
//! unique k-mer ranks, then perform backward search only for those k-mers.
//! This reduces build time from 268M searches to ~2M for a typical
//! microbial genome.
//!
//! File layouts:
//! ```text
//! Dense (K ≤ 12):
//!   offset 0:  [l_0: u32 LE] [r_0: u32 LE]   ← k-mer rank 0
//!   offset 8:  [l_1: u32 LE] [r_1: u32 LE]   ← k-mer rank 1
//!   ...
//!   Total: 4^K × 8 bytes
//!
//! Sparse (K ≥ 13):
//!   offset 0:  magic "NTKS" (4 bytes)
//!   offset 4:  K (u32 LE)
//!   offset 8:  n_entries (u32 LE)
//!   offset 12: padding (4 bytes, zero)
//!   offset 16: [rank_0: u32 LE] [l_0: u32 LE] [r_0: u32 LE]  (sorted by rank)
//!   offset 28: [rank_1: u32 LE] [l_1: u32 LE] [r_1: u32 LE]
//!   ...
//!   Total: 16 + n_entries × 12 bytes
//! ```

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;

use super::simd_search::FmOcc;

/// Small seed length. 4^10 = 1,048,576 entries × 8 bytes = 8 MiB.
/// Used for mm ≤ 2: non-overlapping segments, seed_mm=1, 43 variants max.
pub const SEED_K_SMALL: usize = 10;

/// Large seed length. Sparse format: only genome-present k-mers stored.
/// Used for mm = 3: deeper seed shortens BWT to 6 steps, seed_mm=2.
pub const SEED_K_LARGE: usize = 14;

/// Threshold: K ≤ this uses dense format, K > this uses sparse.
const SPARSE_THRESHOLD: usize = 12;

/// Number of entries in a dense table: 4^K.
const fn table_size(k: usize) -> usize {
    1usize << (2 * k) // 4^k
}

/// Sentinel value for absent k-mers (l > r → empty interval).
const ABSENT_L: u32 = u32::MAX;
const ABSENT_R: u32 = 0;

/// Magic bytes for sparse format header.
const SPARSE_MAGIC: [u8; 4] = *b"NTKS";

/// Sparse entry size: rank(4) + l(4) + r(4) = 12 bytes.
const SPARSE_ENTRY_SIZE: usize = 12;

/// Sparse header size: magic(4) + k(4) + n_entries(4) + padding(4) = 16 bytes.
const SPARSE_HEADER_SIZE: usize = 16;

/// Encode a DNA base to 2-bit rank: A=0, C=1, G=2, T=3.
/// Returns `None` for non-ACGT (N, $, etc.).
#[inline]
fn base_to_rank(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Decode 2-bit rank back to ASCII base.
#[inline]
pub fn rank_to_base(r: u8) -> u8 {
    [b'A', b'C', b'G', b'T'][r as usize]
}

/// Compute the rank (0-based table offset) of a k-mer sequence.
/// Returns `None` if any base is not in {A, C, G, T}.
pub fn kmer_to_rank(kmer: &[u8]) -> Option<usize> {
    let mut rank = 0usize;
    for &b in kmer {
        rank = (rank << 2) | base_to_rank(b)? as usize;
    }
    Some(rank)
}

/// Decode a rank back to a k-mer byte vector (MSB-first).
fn rank_to_kmer(mut rank: usize, k: usize) -> Vec<u8> {
    let mut kmer = vec![0u8; k];
    for i in (0..k).rev() {
        kmer[i] = rank_to_base((rank & 3) as u8);
        rank >>= 2;
    }
    kmer
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Build: dense (all 4^K) or sparse (genome-present only)
// ═══════════════════════════════════════════════════════════════════════════════

/// Build the k-mer seed table and write it to `out_path`.
/// Automatically selects dense or sparse format based on K.
pub fn build_kmer_index<I: FmOcc>(index: &I, k: usize, out_path: &Path) -> std::io::Result<()> {
    if k <= SPARSE_THRESHOLD {
        build_dense(index, k, out_path)
    } else {
        build_sparse(index, k, out_path)
    }
}

/// Build a k-mer seed table from genome text + FM-index.
/// Scans the text for present k-mers, then searches only those.
pub fn build_kmer_index_from_text<I: FmOcc>(
    index: &I,
    text: &[u8],
    k: usize,
    out_path: &Path,
) -> std::io::Result<()> {
    if k <= SPARSE_THRESHOLD {
        build_dense(index, k, out_path)
    } else {
        build_sparse_from_text(index, text, k, out_path)
    }
}

/// Dense build: enumerate all 4^K k-mers, exact backward search.
fn build_dense<I: FmOcc>(index: &I, k: usize, out_path: &Path) -> std::io::Result<()> {
    let n_entries = table_size(k);
    let file = File::create(out_path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    let sa_len = index.sa_len();

    for rank in 0..n_entries {
        let kmer = rank_to_kmer(rank, k);
        let (l, r) = exact_backward_search(index, &kmer, sa_len);
        w.write_all(&(l as u32).to_le_bytes())?;
        w.write_all(&(r as u32).to_le_bytes())?;
    }

    w.flush()?;
    Ok(())
}

/// Sparse build: enumerate all 4^K k-mers (no text available).
/// Falls back to full enumeration but stores only present entries.
fn build_sparse<I: FmOcc>(index: &I, k: usize, out_path: &Path) -> std::io::Result<()> {
    let n_entries = table_size(k);
    let sa_len = index.sa_len();

    // Collect present k-mers
    let mut entries: Vec<(u32, u32, u32)> = Vec::new();
    for rank in 0..n_entries {
        let kmer = rank_to_kmer(rank, k);
        let (l, r) = exact_backward_search(index, &kmer, sa_len);
        if l as u32 != ABSENT_L || r as u32 != ABSENT_R {
            entries.push((rank as u32, l as u32, r as u32));
        }
    }

    write_sparse_file(&entries, k, out_path)
}

/// Sparse build from genome text: scan for present k-mers first,
/// then backward-search only those. Much faster for large K.
fn build_sparse_from_text<I: FmOcc>(
    index: &I,
    text: &[u8],
    k: usize,
    out_path: &Path,
) -> std::io::Result<()> {
    // Scan text for unique k-mer ranks
    let mut seen = vec![false; table_size(k)];
    let mut unique_ranks: Vec<usize> = Vec::new();

    if text.len() >= k {
        for i in 0..=text.len() - k {
            if let Some(rank) = kmer_to_rank(&text[i..i + k]) {
                if !seen[rank] {
                    seen[rank] = true;
                    unique_ranks.push(rank);
                }
            }
        }
    }
    drop(seen);
    unique_ranks.sort_unstable();

    // Backward-search only present k-mers
    let sa_len = index.sa_len();
    let mut entries: Vec<(u32, u32, u32)> = Vec::with_capacity(unique_ranks.len());

    for &rank in &unique_ranks {
        let kmer = rank_to_kmer(rank, k);
        let (l, r) = exact_backward_search(index, &kmer, sa_len);
        if l as u32 != ABSENT_L || r as u32 != ABSENT_R {
            entries.push((rank as u32, l as u32, r as u32));
        }
    }

    write_sparse_file(&entries, k, out_path)
}

/// Write sparse entries to disk with header.
fn write_sparse_file(
    entries: &[(u32, u32, u32)],
    k: usize,
    out_path: &Path,
) -> std::io::Result<()> {
    let file = File::create(out_path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    // Header
    w.write_all(&SPARSE_MAGIC)?;
    w.write_all(&(k as u32).to_le_bytes())?;
    w.write_all(&(entries.len() as u32).to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?; // padding

    // Sorted entries
    for &(rank, l, r) in entries {
        w.write_all(&rank.to_le_bytes())?;
        w.write_all(&l.to_le_bytes())?;
        w.write_all(&r.to_le_bytes())?;
    }

    w.flush()?;
    Ok(())
}

/// Exact backward search for a k-mer. Returns inclusive BWT interval (l, r).
/// If absent, returns (usize::MAX, 0) so l > r signals empty.
fn exact_backward_search<I: FmOcc>(index: &I, kmer: &[u8], sa_len: usize) -> (usize, usize) {
    let mut l: usize = 0;
    let mut r: usize = sa_len - 1;

    // Consume right-to-left (same order as the width-first engine).
    for &base in kmer.iter().rev() {
        let less_b = index.less(base);
        let occ_l = if l > 0 { index.occ(l - 1, base) } else { 0 };
        let occ_r = index.occ(r, base);

        let nl = less_b + occ_l;
        let nr_exc = less_b + occ_r;

        if nl >= nr_exc {
            return (ABSENT_L as usize, ABSENT_R as usize);
        }
        l = nl;
        r = nr_exc - 1; // back to inclusive
    }

    (l, r)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Load: memory-mapped read-only table (auto-detects format)
// ═══════════════════════════════════════════════════════════════════════════════

/// Storage variant for the loaded table.
enum TableFormat {
    /// Flat array: offset = rank × 8.
    Dense,
    /// Sorted (rank, l, r) triples with binary search.
    Sparse { n_entries: usize },
}

/// Backing storage: either mmap (file-loaded) or in-memory (SA sweep).
enum TableData {
    Mmap(Mmap),
    Memory(Vec<u8>),
}

impl TableData {
    #[inline]
    fn as_bytes(&self) -> &[u8] {
        match self {
            TableData::Mmap(m) => m,
            TableData::Memory(v) => v,
        }
    }
}

/// K-mer seed table. Supports both dense and sparse formats,
/// backed by either mmap (file-loaded) or in-memory data (SA sweep).
pub struct KmerSeedTable {
    data: TableData,
    k: usize,
    format: TableFormat,
}

// SAFETY: Mmap is Send+Sync (read-only view of a file). Vec<u8> is trivially Send+Sync.
unsafe impl Send for KmerSeedTable {}
unsafe impl Sync for KmerSeedTable {}

impl KmerSeedTable {
    /// Open an existing k-mer index file via `mmap`. Auto-detects format.
    pub fn open(path: &Path, k: usize) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        let format = Self::detect_format(mmap.as_ref(), k)?;
        Ok(KmerSeedTable { data: TableData::Mmap(mmap), k, format })
    }

    /// Detect whether byte data is dense or sparse format.
    fn detect_format(bytes: &[u8], k: usize) -> std::io::Result<TableFormat> {
        // Check for sparse magic header
        if bytes.len() >= SPARSE_HEADER_SIZE && bytes[0..4] == SPARSE_MAGIC {
            let stored_k = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize;
            let n_entries =
                u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;

            if stored_k != k {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("sparse index k mismatch: file has k={}, expected k={}", stored_k, k),
                ));
            }

            let expected = SPARSE_HEADER_SIZE + n_entries * SPARSE_ENTRY_SIZE;
            if bytes.len() != expected {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "sparse index size mismatch: expected {} bytes for {} entries, got {}",
                        expected, n_entries, bytes.len()
                    ),
                ));
            }

            Ok(TableFormat::Sparse { n_entries })
        } else {
            // Dense format: validate size
            let expected = table_size(k) * 8;
            if bytes.len() != expected {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "dense index size mismatch: expected {} bytes for k={}, got {}",
                        expected, k, bytes.len()
                    ),
                ));
            }
            Ok(TableFormat::Dense)
        }
    }

    /// Build (if missing) then open.
    pub fn open_or_build<I: FmOcc>(
        index: &I,
        fasta_path: &str,
        k: usize,
    ) -> std::io::Result<Self> {
        let idx_path = kmer_index_path(fasta_path, k);
        if !idx_path.exists() {
            build_kmer_index(index, k, &idx_path)?;
        }
        Self::open(&idx_path, k)
    }

    /// Build from genome text (if missing) then open. Uses sparse scan for large K.
    pub fn open_or_build_from_text<I: FmOcc>(
        index: &I,
        text: &[u8],
        fasta_path: &str,
        k: usize,
    ) -> std::io::Result<Self> {
        let idx_path = kmer_index_path(fasta_path, k);
        if !idx_path.exists() {
            build_kmer_index_from_text(index, text, k, &idx_path)?;
        }
        Self::open(&idx_path, k)
    }

    /// Seed length K.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Look up the BWT interval for a k-mer by its 2-bit rank.
    /// Returns `None` if the k-mer is absent in the genome.
    #[inline]
    pub fn lookup_rank(&self, rank: usize) -> Option<(u32, u32)> {
        match self.format {
            TableFormat::Dense => self.lookup_rank_dense(rank),
            TableFormat::Sparse { n_entries } => self.lookup_rank_sparse(rank, n_entries),
        }
    }

    #[inline]
    fn lookup_rank_dense(&self, rank: usize) -> Option<(u32, u32)> {
        let offset = rank * 8;
        let bytes = &self.data.as_bytes()[offset..offset + 8];

        let l = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        let r = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);

        if l == ABSENT_L && r == ABSENT_R {
            None
        } else {
            Some((l, r))
        }
    }

    #[inline]
    fn lookup_rank_sparse(&self, rank: usize, n_entries: usize) -> Option<(u32, u32)> {
        let rank = rank as u32;
        let data = &self.data.as_bytes()[SPARSE_HEADER_SIZE..];

        // Binary search over sorted (rank, l, r) triples
        let mut lo = 0usize;
        let mut hi = n_entries;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let off = mid * SPARSE_ENTRY_SIZE;
            let entry_rank =
                u32::from_le_bytes([data[off], data[off + 1], data[off + 2], data[off + 3]]);

            if entry_rank < rank {
                lo = mid + 1;
            } else if entry_rank > rank {
                hi = mid;
            } else {
                let l = u32::from_le_bytes(
                    [data[off + 4], data[off + 5], data[off + 6], data[off + 7]],
                );
                let r = u32::from_le_bytes(
                    [data[off + 8], data[off + 9], data[off + 10], data[off + 11]],
                );
                return Some((l, r));
            }
        }
        None
    }

    /// Look up a k-mer byte slice. Returns `None` if absent or contains non-ACGT.
    #[inline]
    pub fn lookup(&self, kmer: &[u8]) -> Option<(u32, u32)> {
        let rank = kmer_to_rank(kmer)?;
        self.lookup_rank(rank)
    }

    /// Number of present k-mers in the table.
    pub fn n_entries(&self) -> usize {
        match self.format {
            TableFormat::Dense => {
                // Count non-absent entries (expensive, mainly for diagnostics)
                let n = table_size(self.k);
                let mut count = 0;
                for rank in 0..n {
                    if self.lookup_rank_dense(rank).is_some() {
                        count += 1;
                    }
                }
                count
            }
            TableFormat::Sparse { n_entries } => n_entries,
        }
    }

    /// Whether this table uses the sparse format.
    pub fn is_sparse(&self) -> bool {
        matches!(self.format, TableFormat::Sparse { .. })
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Position Table — maps k-mer rank → sorted genome positions (no SA needed)
// ═══════════════════════════════════════════════════════════════════════════════

/// In-memory table mapping each k-mer rank to a sorted list of genome positions
/// where that k-mer occurs. Used by the verify phase to avoid random SA lookups.
///
/// For K ≤ 12: flat offset array (4^K × 4 bytes) + position data.
/// For K ≥ 13: hash map from rank to position slice boundaries.
pub struct PosTable {
    inner: PosTableInner,
    k: usize,
}

enum PosTableInner {
    /// Dense: offset table (4^K + 1 entries) + position data.
    Dense {
        offsets: Vec<u32>,
        positions: Vec<u32>,
    },
    /// Sparse: sorted (rank, offset_start, offset_end) + position data.
    Sparse {
        /// Sorted by rank for binary search.
        index: Vec<(u32, u32, u32)>,
        positions: Vec<u32>,
    },
}

impl PosTable {
    /// Build from the concatenated genome text. Two-pass: count then fill.
    pub fn build(text: &[u8], k: usize) -> Self {
        if k <= SPARSE_THRESHOLD {
            Self::build_dense(text, k)
        } else {
            Self::build_sparse(text, k)
        }
    }

    fn build_dense(text: &[u8], k: usize) -> Self {
        let n_entries = table_size(k);

        // First pass: count occurrences per k-mer rank.
        let mut counts = vec![0u32; n_entries];
        if text.len() >= k {
            for i in 0..=text.len() - k {
                if let Some(rank) = kmer_to_rank(&text[i..i + k]) {
                    counts[rank] += 1;
                }
            }
        }

        // Build offset table (exclusive prefix sums).
        let mut offsets = vec![0u32; n_entries + 1];
        for i in 0..n_entries {
            offsets[i + 1] = offsets[i] + counts[i];
        }
        let total = offsets[n_entries] as usize;

        // Second pass: fill positions.
        let mut positions = vec![0u32; total];
        let mut cursors: Vec<u32> = offsets[..n_entries].to_vec();
        if text.len() >= k {
            for i in 0..=text.len() - k {
                if let Some(rank) = kmer_to_rank(&text[i..i + k]) {
                    let c = cursors[rank] as usize;
                    positions[c] = i as u32;
                    cursors[rank] += 1;
                }
            }
        }

        PosTable {
            inner: PosTableInner::Dense { offsets, positions },
            k,
        }
    }

    fn build_sparse(text: &[u8], k: usize) -> Self {
        use std::collections::HashMap;

        // Single pass: collect positions per rank using a hash map.
        let mut map: HashMap<u32, Vec<u32>> = HashMap::new();
        if text.len() >= k {
            for i in 0..=text.len() - k {
                if let Some(rank) = kmer_to_rank(&text[i..i + k]) {
                    map.entry(rank as u32).or_default().push(i as u32);
                }
            }
        }

        // Sort ranks and flatten into contiguous position array.
        let mut sorted_ranks: Vec<u32> = map.keys().copied().collect();
        sorted_ranks.sort_unstable();

        let total: usize = map.values().map(|v| v.len()).sum();
        let mut positions = Vec::with_capacity(total);
        let mut index = Vec::with_capacity(sorted_ranks.len());

        for &rank in &sorted_ranks {
            let start = positions.len() as u32;
            let poses = map.get(&rank).unwrap();
            positions.extend_from_slice(poses);
            let end = positions.len() as u32;
            index.push((rank, start, end));
        }

        PosTable {
            inner: PosTableInner::Sparse { index, positions },
            k,
        }
    }

    /// Genome positions where the k-mer with the given rank occurs.
    #[inline]
    pub fn positions_for_rank(&self, rank: usize) -> &[u32] {
        match &self.inner {
            PosTableInner::Dense { offsets, positions } => {
                let start = offsets[rank] as usize;
                let end = offsets[rank + 1] as usize;
                &positions[start..end]
            }
            PosTableInner::Sparse { index, positions } => {
                let rank = rank as u32;
                match index.binary_search_by_key(&rank, |&(r, _, _)| r) {
                    Ok(i) => {
                        let (_, start, end) = index[i];
                        &positions[start as usize..end as usize]
                    }
                    Err(_) => &[],
                }
            }
        }
    }

    /// Seed length K.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  SA Sweep: O(N) construction of both KmerSeedTable and PosTable
// ═══════════════════════════════════════════════════════════════════════════════

/// Build both KmerSeedTable and PosTable from a single O(N) sweep over the
/// suffix array. Eliminates all backward searches and HashMap allocations.
///
/// Since the SA is lexicographically sorted, all occurrences of a k-mer form
/// a contiguous run. One linear scan identifies every run boundary, yielding
/// the exact BWT interval [l, r] (inclusive) for each k-mer — the same values
/// that backward search would produce.
///
/// Genome positions are extracted directly from the SA slice within each run,
/// sorted per k-mer to match the PosTable contract.
pub fn sa_sweep_build(sa: &[usize], text: &[u8], k: usize) -> (KmerSeedTable, PosTable) {
    let n = sa.len();

    // Phase 1: Identify contiguous runs of identical k-mers in the sorted SA.
    // Each run → (kmer_rank, sa_start_inclusive, sa_end_exclusive).
    let mut runs: Vec<(usize, usize, usize)> = Vec::new();
    let mut cur_rank: Option<usize> = None;
    let mut run_start: usize = 0;

    for i in 0..n {
        let pos = sa[i];
        let this_rank = if pos + k <= text.len() {
            kmer_to_rank(&text[pos..pos + k])
        } else {
            None
        };

        match (cur_rank, this_rank) {
            (Some(cr), Some(tr)) if cr == tr => {} // same k-mer, extend run
            (Some(cr), _) => {
                // run ended
                runs.push((cr, run_start, i));
                if let Some(tr) = this_rank {
                    cur_rank = Some(tr);
                    run_start = i;
                } else {
                    cur_rank = None;
                }
            }
            (None, Some(tr)) => {
                cur_rank = Some(tr);
                run_start = i;
            }
            (None, None) => {}
        }
    }
    if let Some(cr) = cur_rank {
        runs.push((cr, run_start, n));
    }

    // Phase 2: Build both tables from runs.
    let is_dense = k <= SPARSE_THRESHOLD;

    // ── KmerSeedTable ──────────────────────────────────────────────────────
    let kmer_table = if is_dense {
        let n_entries = table_size(k);
        let mut bytes = vec![0u8; n_entries * 8];
        // Initialize all to ABSENT
        for rank_idx in 0..n_entries {
            let off = rank_idx * 8;
            bytes[off..off + 4].copy_from_slice(&ABSENT_L.to_le_bytes());
            bytes[off + 4..off + 8].copy_from_slice(&ABSENT_R.to_le_bytes());
        }
        // Fill from runs: (l, r) are SA indices, inclusive
        for &(rank, sa_start, sa_end) in &runs {
            let off = rank * 8;
            bytes[off..off + 4].copy_from_slice(&(sa_start as u32).to_le_bytes());
            bytes[off + 4..off + 8].copy_from_slice(&((sa_end - 1) as u32).to_le_bytes());
        }
        KmerSeedTable {
            data: TableData::Memory(bytes),
            k,
            format: TableFormat::Dense,
        }
    } else {
        // Sparse: NTKS header + sorted (rank, l, r) triples
        let mut bytes = Vec::with_capacity(SPARSE_HEADER_SIZE + runs.len() * SPARSE_ENTRY_SIZE);
        bytes.extend_from_slice(&SPARSE_MAGIC);
        bytes.extend_from_slice(&(k as u32).to_le_bytes());
        bytes.extend_from_slice(&(runs.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes()); // padding
        // Runs are already sorted by rank (SA is sorted → k-mers appear in lex order)
        for &(rank, sa_start, sa_end) in &runs {
            bytes.extend_from_slice(&(rank as u32).to_le_bytes());
            bytes.extend_from_slice(&(sa_start as u32).to_le_bytes());
            bytes.extend_from_slice(&((sa_end - 1) as u32).to_le_bytes());
        }
        KmerSeedTable {
            data: TableData::Memory(bytes),
            k,
            format: TableFormat::Sparse { n_entries: runs.len() },
        }
    };

    // ── PosTable ───────────────────────────────────────────────────────────
    let total_positions: usize = runs.iter().map(|&(_, s, e)| e - s).sum();

    let pos_table = if is_dense {
        let n_entries = table_size(k);
        // Build offsets: prefix sum of run lengths, slotted by rank
        let mut offsets = vec![0u32; n_entries + 1];
        for &(rank, sa_start, sa_end) in &runs {
            offsets[rank] = (sa_end - sa_start) as u32;
        }
        // Convert counts to prefix sums
        let mut acc = 0u32;
        for i in 0..n_entries {
            let count = offsets[i];
            offsets[i] = acc;
            acc += count;
        }
        offsets[n_entries] = acc;

        // Fill positions from SA, sorted per k-mer
        let mut positions = vec![0u32; total_positions];
        for &(rank, sa_start, sa_end) in &runs {
            let dst_start = offsets[rank] as usize;
            let run_len = sa_end - sa_start;
            for j in 0..run_len {
                positions[dst_start + j] = sa[sa_start + j] as u32;
            }
            // Sort this k-mer's positions (genome order)
            positions[dst_start..dst_start + run_len].sort_unstable();
        }

        PosTable {
            inner: PosTableInner::Dense { offsets, positions },
            k,
        }
    } else {
        // Sparse: (rank, offset_start, offset_end) index + flat positions
        let mut index = Vec::with_capacity(runs.len());
        let mut positions = Vec::with_capacity(total_positions);

        for &(rank, sa_start, sa_end) in &runs {
            let dst_start = positions.len() as u32;
            let run_len = sa_end - sa_start;
            // Copy genome positions from SA
            let pos_start = positions.len();
            for j in 0..run_len {
                positions.push(sa[sa_start + j] as u32);
            }
            // Sort genome positions for this k-mer
            positions[pos_start..].sort_unstable();
            let dst_end = positions.len() as u32;
            index.push((rank as u32, dst_start, dst_end));
        }

        PosTable {
            inner: PosTableInner::Sparse { index, positions },
            k,
        }
    };

    (kmer_table, pos_table)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Utility
// ═══════════════════════════════════════════════════════════════════════════════

/// Derive the k-mer index path from a FASTA path: `genome.fa` → `genome.10mer.idx`.
pub fn kmer_index_path(fasta_path: &str, k: usize) -> PathBuf {
    let p = Path::new(fasta_path);
    let stem = p.file_stem().unwrap_or_default().to_string_lossy();
    let parent = p.parent().unwrap_or_else(|| Path::new("."));
    parent.join(format!("{}.{}mer.idx", stem, k))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_rank_roundtrip() {
        for k in [4, 10, 14] {
            for rank in [0, 1, 42, table_size(k) - 1] {
                let kmer = rank_to_kmer(rank, k);
                assert_eq!(kmer.len(), k);
                assert_eq!(kmer_to_rank(&kmer), Some(rank));
            }
        }
    }

    #[test]
    fn test_sparse_write_read() {
        let dir = std::env::temp_dir().join("needletail_test_sparse");
        let _ = fs::create_dir_all(&dir);
        let path = dir.join("test.14mer.idx");

        let entries = vec![
            (100u32, 10u32, 20u32),
            (200, 30, 40),
            (500, 50, 60),
        ];
        write_sparse_file(&entries, 14, &path).unwrap();

        let table = KmerSeedTable::open(&path, 14).unwrap();
        assert!(table.is_sparse());
        assert_eq!(table.n_entries(), 3);

        assert_eq!(table.lookup_rank(100), Some((10, 20)));
        assert_eq!(table.lookup_rank(200), Some((30, 40)));
        assert_eq!(table.lookup_rank(500), Some((50, 60)));
        assert_eq!(table.lookup_rank(0), None);
        assert_eq!(table.lookup_rank(150), None);
        assert_eq!(table.lookup_rank(999), None);

        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_pos_table_sparse() {
        // Simulate a small genome: ACGTACGT (8 bp)
        let text = b"ACGTACGT";
        let k = 4;

        let dense = PosTable::build_dense(text, k);
        let sparse = PosTable::build_sparse(text, k);

        // Both should give same results for all k-mers
        let n_entries = table_size(k);
        for rank in 0..n_entries {
            let d = dense.positions_for_rank(rank);
            let s = sparse.positions_for_rank(rank);
            assert_eq!(d, s, "mismatch at rank {}", rank);
        }
    }

    #[test]
    fn test_sa_sweep_vs_old_dense() {
        // Build a small genome with sentinel
        let text = b"ACGTACGTCCGGTTAA$";
        let k = 4;

        // Build suffix array
        let sa = bio::data_structures::suffix_array::suffix_array(text);

        // SA sweep
        let (sweep_kmer, sweep_pos) = sa_sweep_build(&sa, text, k);

        // Old builder: PosTable
        let old_pos = PosTable::build_dense(text, k);

        // Compare PosTable: all ranks should match
        let n_entries = table_size(k);
        for rank in 0..n_entries {
            let sweep_p = sweep_pos.positions_for_rank(rank);
            let old_p = old_pos.positions_for_rank(rank);
            assert_eq!(
                sweep_p, old_p,
                "PosTable mismatch at rank {} (k-mer {:?})",
                rank,
                rank_to_kmer(rank, k)
            );
        }

        // Verify KmerSeedTable: SA intervals should be consistent with positions
        for rank in 0..n_entries {
            let positions = sweep_pos.positions_for_rank(rank);
            match sweep_kmer.lookup_rank(rank) {
                Some((l, r)) => {
                    let interval_size = (r - l + 1) as usize;
                    assert_eq!(
                        interval_size,
                        positions.len(),
                        "KmerSeedTable interval size mismatch at rank {}",
                        rank
                    );
                }
                None => {
                    assert!(
                        positions.is_empty(),
                        "KmerSeedTable says absent but PosTable has {} positions at rank {}",
                        positions.len(),
                        rank
                    );
                }
            }
        }
    }

    #[test]
    fn test_sa_sweep_vs_old_sparse() {
        // Build a genome large enough for sparse (K=14 uses sparse format)
        // Use a repeated pattern to get enough length
        let mut text = Vec::new();
        for _ in 0..100 {
            text.extend_from_slice(b"ACGTACGTCCGGTTAANNACGTGCTAGCTAGC");
        }
        text.push(b'$');
        let k = 14;

        let sa = bio::data_structures::suffix_array::suffix_array(&text);

        let (sweep_kmer, sweep_pos) = sa_sweep_build(&sa, &text, k);
        let old_pos = PosTable::build_sparse(&text, k);

        // Collect all non-empty ranks from old builder
        for rank in 0..table_size(k).min(268_435_456) {
            let old_p = old_pos.positions_for_rank(rank);
            let sweep_p = sweep_pos.positions_for_rank(rank);

            if !old_p.is_empty() || !sweep_p.is_empty() {
                // Sort both for comparison (old sparse doesn't guarantee order)
                let mut old_sorted = old_p.to_vec();
                let mut sweep_sorted = sweep_p.to_vec();
                old_sorted.sort_unstable();
                sweep_sorted.sort_unstable();
                assert_eq!(
                    old_sorted, sweep_sorted,
                    "PosTable mismatch at rank {} (k={})",
                    rank, k
                );
            }
        }

        // Verify KmerSeedTable consistency
        assert!(sweep_kmer.is_sparse());
    }
}
