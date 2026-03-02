//! FmIndexSearcher — Cache-Line Aligned Rank Blocks over a rust-bio BWT.
//!
//! Builds a BWT + suffix array from a FASTA file (via rust-bio), then
//! constructs **cache-line-aligned rank blocks** for the hot-path LF-mapping.
//!
//! ## Memory layout: BlockRank
//!
//! Each `RankBlock` is exactly 64 bytes (1 cache line), covering 64 BWT positions:
//! ```text
//! Offset  0: counts[4]  — absolute cumulative A,C,G,T before block start
//! Offset 16: bv_lo[4]   — per-base bitvectors for positions 0..31
//! Offset 32: bv_hi[4]   — per-base bitvectors for positions 32..63
//! Offset 48: mid[4]     — absolute cumulative A,C,G,T at block_start + 32
//! ```
//!
//! Memory: `ceil(bwt_len / 64) × 64 bytes`. For SacCer3 (12M bp): ~12 MB.
//! Fits entirely in L3 cache, reducing rank latency from ~100ns (DRAM) to ~10ns.
//!
//! Compare to the previous `InterleavedRank` (16 bytes/pos): 192 MB for SacCer3.
//! **16× compression**, same numeric results.

use bio::alphabets::Alphabet;
use bio::data_structures::bwt::{bwt, less};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::io::fasta;

use crate::error::SearchError;
use super::simd_search::{ChromGeometry, FmOcc, BASES};

// ═══════════════════════════════════════════════════════════════════════════════
//  Cache-Line Aligned Rank Blocks
// ═══════════════════════════════════════════════════════════════════════════════

/// Number of BWT positions per rank block.
pub const BLOCK_SIZE: usize = 64;

/// Cache-line-aligned rank block covering 64 consecutive BWT positions.
///
/// Every 16-byte field is a natural `__m128i` SIMD vector with one value per
/// DNA base (A=0, C=1, G=2, T=3). The struct is exactly 64 bytes, aligned to
/// a cache-line boundary, so a single memory fetch loads all data needed for
/// a rank query at any position within the block.
///
/// ```text
/// Field    Type       Bytes  Offset   Content
/// ─────    ─────────  ─────  ──────   ───────
/// counts   [u32; 4]   16     0        abs cumulative A,C,G,T before block start
/// bv_lo    [u32; 4]   16     16       bitvectors for positions 0..31
/// bv_hi    [u32; 4]   16     32       bitvectors for positions 32..63
/// mid      [u32; 4]   16     48       abs cumulative A,C,G,T at block_start + 32
/// ─────────────────────────────────
/// Total:               64 bytes = 1 cache line
/// ```
#[repr(C, align(64))]
#[derive(Clone, Copy)]
pub struct RankBlock {
    /// Absolute cumulative counts before this block:
    /// `counts[b] = count of base b in BWT[0..block_start]`.
    pub counts: [u32; 4],

    /// Per-base bitvectors for positions 0..31 within the block.
    /// `bv_lo[b]` bit `j` = 1 iff `BWT[block_start + j] == base b`.
    pub bv_lo: [u32; 4],

    /// Per-base bitvectors for positions 32..63 within the block.
    /// `bv_hi[b]` bit `j` = 1 iff `BWT[block_start + 32 + j] == base b`.
    pub bv_hi: [u32; 4],

    /// Absolute cumulative counts at block_start + 32:
    /// `mid[b] = counts[b] + popcount(bv_lo[b])`.
    pub mid: [u32; 4],
}

impl RankBlock {
    const fn zeroed() -> Self {
        RankBlock {
            counts: [0; 4],
            bv_lo: [0; 4],
            bv_hi: [0; 4],
            mid: [0; 4],
        }
    }
}

/// Cache-optimized rank structure using 64-byte aligned blocks.
///
/// Replaces the per-position `InterleavedRank` (16 bytes/pos = 192 MB for SacCer3)
/// with sampled blocks (1 byte/pos = 12 MB), fitting entirely in L3 cache.
pub struct BlockRank {
    blocks: Vec<RankBlock>,
    less: [usize; 256],
    bwt_len: usize,
}

impl BlockRank {
    /// Build from a BWT byte sequence and the `less` table produced by rust-bio.
    fn from_bwt_and_less(bwt_seq: &[u8], less_tbl: &[usize]) -> Self {
        let n = bwt_seq.len();
        let n_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
        let mut blocks = vec![RankBlock::zeroed(); n_blocks];

        let mut cumul = [0u32; 4];

        for b in 0..n_blocks {
            let block_start = b * BLOCK_SIZE;
            let block_end = (block_start + BLOCK_SIZE).min(n);

            // Store absolute cumulative counts before this block.
            blocks[b].counts = cumul;

            // Fill bv_lo (positions 0..31 within block) and accumulate.
            let lo_end = (block_start + 32).min(block_end);
            for pos in block_start..lo_end {
                let j = pos - block_start;
                let bi = base_idx(bwt_seq[pos]);
                if bi < 4 {
                    blocks[b].bv_lo[bi] |= 1u32 << j;
                    cumul[bi] += 1;
                }
            }

            // Store mid-block cumulative counts (at block_start + 32).
            blocks[b].mid = cumul;

            // Fill bv_hi (positions 32..63 within block) and accumulate.
            let hi_start = block_start + 32;
            if hi_start < block_end {
                for pos in hi_start..block_end {
                    let j = pos - hi_start;
                    let bi = base_idx(bwt_seq[pos]);
                    if bi < 4 {
                        blocks[b].bv_hi[bi] |= 1u32 << j;
                        cumul[bi] += 1;
                    }
                }
            }
        }

        let mut less_arr = [0usize; 256];
        for (i, &v) in less_tbl.iter().enumerate() {
            less_arr[i] = v;
        }

        BlockRank {
            blocks,
            less: less_arr,
            bwt_len: n,
        }
    }

    /// Compute rank for all 4 bases at a given position.
    /// Returns `[count_A, count_C, count_G, count_T]` in `BWT[0..=pos]`.
    #[inline]
    pub fn rank_all(&self, pos: usize) -> [u32; 4] {
        let block = &self.blocks[pos / BLOCK_SIZE];
        let offset = pos % BLOCK_SIZE;

        if offset < 32 {
            let mask = inclusive_mask(offset);
            [
                block.counts[0] + (block.bv_lo[0] & mask).count_ones(),
                block.counts[1] + (block.bv_lo[1] & mask).count_ones(),
                block.counts[2] + (block.bv_lo[2] & mask).count_ones(),
                block.counts[3] + (block.bv_lo[3] & mask).count_ones(),
            ]
        } else {
            let mask = inclusive_mask(offset - 32);
            [
                block.mid[0] + (block.bv_hi[0] & mask).count_ones(),
                block.mid[1] + (block.bv_hi[1] & mask).count_ones(),
                block.mid[2] + (block.bv_hi[2] & mask).count_ones(),
                block.mid[3] + (block.bv_hi[3] & mask).count_ones(),
            ]
        }
    }

    /// Single occ lookup: count of `base` in BWT[0..=pos].
    #[inline]
    pub fn occ(&self, pos: usize, base: u8) -> u32 {
        let bi = base_idx(base);
        if bi >= 4 {
            return 0;
        }
        let block = &self.blocks[pos / BLOCK_SIZE];
        let offset = pos % BLOCK_SIZE;

        if offset < 32 {
            let mask = inclusive_mask(offset);
            block.counts[bi] + (block.bv_lo[bi] & mask).count_ones()
        } else {
            let mask = inclusive_mask(offset - 32);
            block.mid[bi] + (block.bv_hi[bi] & mask).count_ones()
        }
    }

    /// Batch LF-mapping: compute all 4 `(nl, nr_exclusive)` for A, C, G, T.
    #[inline]
    pub fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        let occ_l = if l > 0 {
            self.rank_all(l - 1)
        } else {
            [0u32; 4]
        };
        let occ_r = self.rank_all(r);

        let mut result = [(0u32, 0u32); 4];
        for bi in 0..4 {
            let less_b = self.less[BASES[bi] as usize] as u32;
            result[bi] = (less_b + occ_l[bi], less_b + occ_r[bi]);
        }
        result
    }

    /// Prefetch the cache line(s) needed for a future `lf_map(l, r)` call.
    /// Each block is exactly 1 cache line, so one prefetch covers all 4 bases.
    #[inline]
    pub fn prefetch_lf(&self, l: usize, r: usize) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            if l > 0 {
                let ptr = &self.blocks[(l - 1) / BLOCK_SIZE] as *const RankBlock as *const i8;
                _mm_prefetch(ptr, _MM_HINT_T0);
            }
            let ptr = &self.blocks[r / BLOCK_SIZE] as *const RankBlock as *const i8;
            _mm_prefetch(ptr, _MM_HINT_T0);
        }
    }

    /// Blocks slice for SIMD access.
    #[inline]
    pub fn blocks(&self) -> &[RankBlock] {
        &self.blocks
    }

    /// Less table reference.
    #[inline]
    pub fn less_table(&self) -> &[usize; 256] {
        &self.less
    }

    /// BWT length.
    #[inline]
    pub fn bwt_len(&self) -> usize {
        self.bwt_len
    }

    /// Base pointer for AVX2 gather instructions.
    /// Each `RankBlock` is 16 × i32 = 64 bytes. Gather index for field `f`
    /// in block `b` is `b * 16 + f`.
    #[inline]
    pub fn blocks_i32_ptr(&self) -> *const i32 {
        self.blocks.as_ptr() as *const i32
    }

    /// Convert blocks back to flat interleaved format for serialization.
    /// Returns `data[pos * 4 + base_idx]` = count of base in `BWT[0..=pos]`.
    /// O(N) — only used during index persistence, never on the hot path.
    pub fn to_interleaved_data(&self) -> Vec<u32> {
        let n = self.bwt_len;
        let mut data = vec![0u32; n * 4];

        for b_idx in 0..self.blocks.len() {
            let block = &self.blocks[b_idx];
            let block_start = b_idx * BLOCK_SIZE;
            let block_end = (block_start + BLOCK_SIZE).min(n);

            // Low half: positions 0..31
            let lo_end = (block_start + 32).min(block_end);
            for pos in block_start..lo_end {
                let j = pos - block_start;
                let mask = inclusive_mask(j);
                let base = pos * 4;
                data[base] = block.counts[0] + (block.bv_lo[0] & mask).count_ones();
                data[base + 1] = block.counts[1] + (block.bv_lo[1] & mask).count_ones();
                data[base + 2] = block.counts[2] + (block.bv_lo[2] & mask).count_ones();
                data[base + 3] = block.counts[3] + (block.bv_lo[3] & mask).count_ones();
            }

            // High half: positions 32..63
            let hi_start = block_start + 32;
            if hi_start < block_end {
                for pos in hi_start..block_end {
                    let j = pos - hi_start;
                    let mask = inclusive_mask(j);
                    let base = pos * 4;
                    data[base] = block.mid[0] + (block.bv_hi[0] & mask).count_ones();
                    data[base + 1] = block.mid[1] + (block.bv_hi[1] & mask).count_ones();
                    data[base + 2] = block.mid[2] + (block.bv_hi[2] & mask).count_ones();
                    data[base + 3] = block.mid[3] + (block.bv_hi[3] & mask).count_ones();
                }
            }
        }

        data
    }
}

/// Inclusive bitmask: bits 0 through `offset` are set (inclusive).
/// `offset` must be in 0..32.
#[inline]
fn inclusive_mask(offset: usize) -> u32 {
    debug_assert!(offset < 32);
    // (2u64 << offset) - 1 handles offset=31 without overflow.
    ((2u64 << offset) - 1) as u32
}

/// Map base byte to interleaved index: A=0, C=1, G=2, T=3.
/// Returns 4 for non-ACGT (sentinels, N).
#[inline]
fn base_idx(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // $, N — not counted in rank
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Chromosome metadata
// ═══════════════════════════════════════════════════════════════════════════════

/// Metadata for a single chromosome / FASTA record.
pub struct ChromInfo {
    pub name: String,
    /// Start offset of this chromosome in the concatenated text.
    pub start: usize,
    /// Length of this chromosome's sequence (excluding the '$' sentinel).
    pub len: usize,
}

// ═══════════════════════════════════════════════════════════════════════════════
//  FmIndexSearcher
// ═══════════════════════════════════════════════════════════════════════════════

/// FM-Index over a concatenated multi-FASTA genome.
///
/// Uses `BlockRank` (64-byte cache-line-aligned blocks) for the hot-path
/// LF-mapping. 16× smaller than the previous per-position layout, fits in L3.
pub struct FmIndexSearcher {
    rank: BlockRank,
    sa: RawSuffixArray,
    chroms: Vec<ChromInfo>,
    /// Concatenated genome text (with '$' separators), kept for verify phase.
    text: Vec<u8>,
}

impl FmIndexSearcher {
    /// Build an FM-Index from a pre-built master buffer (zero-copy from Genome).
    ///
    /// Takes ownership of the text buffer. Builds SA, BWT, BlockRank directly.
    /// No file I/O, no uppercasing, no concatenation.
    pub fn from_text(text: Vec<u8>, chroms: Vec<ChromInfo>) -> Result<Self, SearchError> {
        use std::time::Instant;

        if chroms.is_empty() {
            return Err(SearchError::Other(anyhow::anyhow!(
                "No chromosomes provided"
            )));
        }

        let t0 = Instant::now();
        eprintln!("[fm_index] from_text: {} bp, {} chroms",
            text.len(), chroms.len());

        let t1 = Instant::now();
        let alphabet = Alphabet::new(b"$ACGTN");
        let sa = suffix_array(&text);
        eprintln!("[fm_index] suffix_array (SA-IS): {:.3}s", t1.elapsed().as_secs_f64());

        let t2 = Instant::now();
        let bwt_seq = bwt(&text, &sa);
        let less_tbl = less(&bwt_seq, &alphabet);
        eprintln!("[fm_index] BWT + less: {:.3}s", t2.elapsed().as_secs_f64());

        let t3 = Instant::now();
        let rank = BlockRank::from_bwt_and_less(&bwt_seq, &less_tbl);
        eprintln!("[fm_index] BlockRank: {:.3}s ({} blocks)",
            t3.elapsed().as_secs_f64(), rank.blocks.len());

        eprintln!("[fm_index] total: {:.3}s", t0.elapsed().as_secs_f64());

        Ok(FmIndexSearcher {
            rank,
            sa,
            chroms,
            text,
        })
    }

    /// Build an FM-Index from a FASTA file.
    pub fn from_fasta(path: &str) -> Result<Self, SearchError> {
        use std::time::Instant;

        let t0 = Instant::now();
        let reader =
            fasta::Reader::from_file(path).map_err(SearchError::Other)?;

        let mut text: Vec<u8> = Vec::new();
        let mut chroms: Vec<ChromInfo> = Vec::new();

        for result in reader.records() {
            let record = result.map_err(SearchError::Io)?;

            let start = text.len();
            let seq = record.seq();

            let seq_up: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();

            chroms.push(ChromInfo {
                name: record.id().to_string(),
                start,
                len: seq_up.len(),
            });

            text.extend_from_slice(&seq_up);
            text.push(b'$');
        }

        if chroms.is_empty() {
            return Err(SearchError::Other(anyhow::anyhow!(
                "FASTA file contains no records: {}",
                path
            )));
        }

        let t_fasta = t0.elapsed();
        eprintln!("[fm_index] FASTA parse: {:.3}s ({} bp, {} chroms)",
            t_fasta.as_secs_f64(), text.len(), chroms.len());

        Self::from_text(text, chroms)
    }

    /// Flat interleaved rank data for persistence (reconstructed from blocks).
    pub(crate) fn to_interleaved_rank_data(&self) -> Vec<u32> {
        self.rank.to_interleaved_data()
    }

    /// Raw less table for serialization.
    pub(crate) fn less_raw(&self) -> &[usize; 256] {
        self.rank.less_table()
    }

    /// Raw suffix array for serialization.
    pub(crate) fn sa_raw(&self) -> &[usize] {
        &self.sa
    }

    /// Chromosome info for serialization.
    pub(crate) fn chroms(&self) -> &[ChromInfo] {
        &self.chroms
    }

    /// Chromosome names in FASTA order (index = chrom_id).
    pub fn chrom_names(&self) -> Vec<&str> {
        self.chroms.iter().map(|c| c.name.as_str()).collect()
    }

    /// Concatenated genome text for the verify phase.
    pub fn text(&self) -> &[u8] {
        &self.text
    }

    /// Extract chromosome geometry for the width-first search engine.
    pub fn chrom_geometry(&self) -> ChromGeometry {
        ChromGeometry {
            ranges: self.chroms.iter().map(|c| (c.start, c.len)).collect(),
        }
    }
}

// ─── FmOcc bridge ────────────────────────────────────────────────────────────

impl FmOcc for FmIndexSearcher {
    #[inline]
    fn less(&self, c: u8) -> usize {
        self.rank.less[c as usize]
    }

    #[inline]
    fn occ(&self, pos: usize, c: u8) -> usize {
        self.rank.occ(pos, c) as usize
    }

    #[inline]
    fn sa(&self, idx: usize) -> usize {
        self.sa[idx]
    }

    #[inline]
    fn sa_len(&self) -> usize {
        self.sa.len()
    }

    #[inline]
    fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        self.rank.lf_map(l, r)
    }

    #[inline]
    fn prefetch_lf(&self, l: usize, r: usize) {
        self.rank.prefetch_lf(l, r);
    }

    #[inline]
    fn rank_data(&self) -> Option<(&[u32], &[usize; 256])> {
        None // No flat interleaved data — use rank_blocks() for SIMD
    }

    #[inline]
    fn rank_blocks(&self) -> Option<(&[RankBlock], &[usize; 256])> {
        Some((self.rank.blocks(), self.rank.less_table()))
    }

    #[inline]
    fn sa_slice(&self) -> Option<&[usize]> {
        Some(&self.sa)
    }
}
