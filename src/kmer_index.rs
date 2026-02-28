//! On-disk k-mer seed table: precomputed BWT [l, r] intervals for all 4^K k-mers.
//!
//! At build time, we enumerate every possible k-mer of length K (default 10),
//! perform an exact backward search in the FM-Index, and write the resulting
//! BWT interval `(l: u32, r: u32)` to a flat binary file at a deterministic
//! offset: `kmer_rank * 8`. Invalid intervals (absent k-mers) are stored as
//! `(u32::MAX, 0)`.
//!
//! At query time, we `mmap` this file read-only. For each query, we extract
//! its trailing K-mer (the suffix consumed first by backward search), compute
//! mismatch neighborhood variations, look up each variation's BWT interval in
//! O(1), and seed the frontier directly at depth K — bypassing the exponential
//! branching of the first K BWT steps entirely.
//!
//! File layout:
//! ```text
//! offset 0:  [l_0: u32 LE] [r_0: u32 LE]   ← k-mer rank 0 (AAA...A)
//! offset 8:  [l_1: u32 LE] [r_1: u32 LE]   ← k-mer rank 1 (AAA...C)
//! ...
//! offset 8*(4^K - 1): [l_last] [r_last]    ← k-mer rank 4^K-1 (TTT...T)
//! ```
//!
//! Total size for K=10: 4^10 * 8 = 8,388,608 bytes (8 MiB).

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;

use crate::simd_search::FmOcc;

/// Small seed length. 4^10 = 1,048,576 entries × 8 bytes = 8 MiB.
/// Used for mm ≤ 2: non-overlapping segments, seed_mm=1, 43 variants max.
pub const SEED_K_SMALL: usize = 10;

/// Large seed length. 4^14 = 268,435,456 entries × 8 bytes = 2.1 GiB.
/// Used for mm = 3: deeper seed shortens BWT to 6 steps, seed_mm=2.
pub const SEED_K_LARGE: usize = 14;

/// Number of entries in the table: 4^K.
const fn table_size(k: usize) -> usize {
    1usize << (2 * k) // 4^k
}

/// Sentinel value for absent k-mers (l > r → empty interval).
const ABSENT_L: u32 = u32::MAX;
const ABSENT_R: u32 = 0;

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
fn kmer_to_rank(kmer: &[u8]) -> Option<usize> {
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
//  Build: enumerate all 4^K k-mers, exact backward search, write binary file
// ═══════════════════════════════════════════════════════════════════════════════

/// Build the k-mer seed table and write it to `out_path`.
///
/// Performs 4^K exact backward searches. For K=10 on SacCer3 this takes ~2s.
pub fn build_kmer_index<I: FmOcc>(index: &I, k: usize, out_path: &Path) -> std::io::Result<()> {
    let n_entries = table_size(k);
    let file = File::create(out_path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    let sa_len = index.sa_len();

    for rank in 0..n_entries {
        let kmer = rank_to_kmer(rank, k);

        // Exact backward search: consume k-mer right-to-left.
        let (l, r) = exact_backward_search(index, &kmer, sa_len);

        w.write_all(&(l as u32).to_le_bytes())?;
        w.write_all(&(r as u32).to_le_bytes())?;
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
//  Load: memory-mapped read-only table
// ═══════════════════════════════════════════════════════════════════════════════

/// Memory-mapped k-mer seed table. 8 bytes per entry, read-only.
pub struct KmerSeedTable {
    mmap: Mmap,
    k: usize,
}

// SAFETY: Mmap is Send+Sync (read-only view of a file). usize is trivially so.
unsafe impl Send for KmerSeedTable {}
unsafe impl Sync for KmerSeedTable {}

impl KmerSeedTable {
    /// Open an existing k-mer index file via `mmap`.
    pub fn open(path: &Path, k: usize) -> std::io::Result<Self> {
        let expected = table_size(k) * 8;
        let meta = fs::metadata(path)?;
        if meta.len() != expected as u64 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "kmer index size mismatch: expected {} bytes for k={}, got {}",
                    expected,
                    k,
                    meta.len()
                ),
            ));
        }
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        Ok(KmerSeedTable { mmap, k })
    }

    /// Build (if missing) then open. Returns the table and the path used.
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

    /// Seed length K.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Look up the BWT interval for a k-mer by its 2-bit rank.
    /// Returns `None` if the k-mer is absent in the genome.
    #[inline]
    pub fn lookup_rank(&self, rank: usize) -> Option<(u32, u32)> {
        let offset = rank * 8;
        let bytes = &self.mmap[offset..offset + 8];

        let l = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        let r = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]);

        if l == ABSENT_L && r == ABSENT_R {
            None
        } else {
            Some((l, r))
        }
    }

    /// Look up a k-mer byte slice. Returns `None` if absent or contains non-ACGT.
    #[inline]
    pub fn lookup(&self, kmer: &[u8]) -> Option<(u32, u32)> {
        let rank = kmer_to_rank(kmer)?;
        self.lookup_rank(rank)
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Mismatch neighborhood seeding
// ═══════════════════════════════════════════════════════════════════════════════

/// A seed: BWT interval at depth K with mismatch cost already incurred.
#[derive(Clone, Copy)]
pub struct KmerSeed {
    pub l: u32,
    pub r: u32,
    /// Mismatches consumed by the k-mer prefix variations.
    pub mm_used: u8,
}

/// Generate all valid seeds for a query's trailing K-mer, allowing up to
/// `max_mm` mismatches within the K-mer itself.
///
/// The trailing K-mer is `query[query_len - K ..]` (the suffix consumed first
/// by backward search). We enumerate all base substitution combinations up to
/// `max_mm` mismatches, look each up in the seed table, and collect non-absent
/// results.
///
/// For K=10, max_mm=3: at most C(10,0) + C(10,1)*3 + C(10,2)*9 + C(10,3)*27
/// = 1 + 30 + 270 + 3240 = 3541 lookups (each is a single array dereference).
pub fn seed_with_mismatches(
    table: &KmerSeedTable,
    kmer: &[u8],
    max_mm: u8,
    out: &mut Vec<KmerSeed>,
) {
    debug_assert_eq!(kmer.len(), table.k());
    let k = table.k();

    // Stack-based enumeration: (position, current_rank, mm_so_far).
    // We build the rank incrementally left-to-right (MSB first).
    struct Frame {
        pos: usize,
        rank: usize,
        mm: u8,
    }

    let mut stack = Vec::with_capacity(4096);
    stack.push(Frame { pos: 0, rank: 0, mm: 0 });

    while let Some(Frame { pos, rank, mm }) = stack.pop() {
        if pos == k {
            // Complete k-mer: look up in table.
            if let Some((l, r)) = table.lookup_rank(rank) {
                out.push(KmerSeed { l, r, mm_used: mm });
            }
            continue;
        }

        let orig_base = kmer[pos];
        for &sub_rank in &[0u8, 1, 2, 3] {
            let sub_base = rank_to_base(sub_rank);
            let cost = if sub_base == orig_base { 0u8 } else { 1 };
            let new_mm = mm + cost;

            // Prune: can't exceed budget even if all remaining positions match.
            if new_mm > max_mm {
                continue;
            }

            // Even stronger prune: if we've already used all mismatches,
            // only the exact base is viable for remaining positions.
            // (Handled naturally by the budget check above, but this comment
            // documents the pruning property.)

            stack.push(Frame {
                pos: pos + 1,
                rank: (rank << 2) | sub_rank as usize,
                mm: new_mm,
            });
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Position Table — maps k-mer rank → sorted genome positions (no SA needed)
// ═══════════════════════════════════════════════════════════════════════════════

/// In-memory table mapping each k-mer rank to a sorted list of genome positions
/// where that k-mer occurs. Used by the verify phase to avoid random SA lookups.
///
/// Memory: offset table (4^K × 4 bytes) + position data (~text_len × 4 bytes).
/// For SacCer3 with K=10: ~4MB + ~48MB = ~52MB.
pub struct PosTable {
    /// `offsets[rank]` = start index in `positions` for this k-mer rank.
    /// `offsets[4^K]` = total number of positions (sentinel).
    offsets: Vec<u32>,
    /// All genome positions, partitioned by k-mer rank.
    positions: Vec<u32>,
    k: usize,
}

impl PosTable {
    /// Build from the concatenated genome text. Two-pass: count then fill.
    pub fn build(text: &[u8], k: usize) -> Self {
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

        PosTable { offsets, positions, k }
    }

    /// Genome positions where the k-mer with the given rank occurs.
    #[inline]
    pub fn positions_for_rank(&self, rank: usize) -> &[u32] {
        let start = self.offsets[rank] as usize;
        let end = self.offsets[rank + 1] as usize;
        &self.positions[start..end]
    }

    /// Seed length K.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  K-mer mismatch variant enumeration (by rank, for PosTable)
// ═══════════════════════════════════════════════════════════════════════════════

/// Enumerate all k-mer variants with up to `max_mm` mismatches.
/// Outputs `(rank, mm_used)` pairs into `out`.
///
/// Same combinatorial enumeration as `seed_with_mismatches`, but produces
/// k-mer ranks (for PosTable lookup) instead of BWT intervals.
pub fn enumerate_kmer_variants(
    kmer: &[u8],
    max_mm: u8,
    out: &mut Vec<(usize, u8)>,
) {
    let k = kmer.len();

    struct Frame {
        pos: usize,
        rank: usize,
        mm: u8,
    }

    let mut stack = Vec::with_capacity(if max_mm == 0 { 16 } else { 4096 });
    stack.push(Frame { pos: 0, rank: 0, mm: 0 });

    while let Some(Frame { pos, rank, mm }) = stack.pop() {
        if pos == k {
            out.push((rank, mm));
            continue;
        }

        let orig_base = kmer[pos];
        for &sub_rank in &[0u8, 1, 2, 3] {
            let sub_base = rank_to_base(sub_rank);
            let cost = if sub_base == orig_base { 0u8 } else { 1 };
            let new_mm = mm + cost;
            if new_mm > max_mm {
                continue;
            }
            stack.push(Frame {
                pos: pos + 1,
                rank: (rank << 2) | sub_rank as usize,
                mm: new_mm,
            });
        }
    }
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
