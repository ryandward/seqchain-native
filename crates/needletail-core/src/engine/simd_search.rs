//! Width-first SIMD FM-Index search with mmap-backed out-of-core frontier.
//!
//! Three domains:
//!
//! 1. **MmapFrontier** — SoA columns backed by memory-mapped temp files.
//!    The OS page cache manages physical RAM; our footprint stays flat even
//!    if the frontier balloons to 80M+ rows.
//!
//! 2. **KmerSeedTable** — Precomputed 4^K BWT intervals on disk (see kmer_index.rs).
//!    Queries seed the frontier at depth K instead of depth 0, bypassing the
//!    exponential branching of the first K BWT steps.
//!
//! 3. **step_depth** — Branchless 4-way LF expansion streaming from one mmap
//!    frontier into the next.
//!
//! Invariant: all queries must be the same length.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Seek, SeekFrom, Write};
use std::slice;

use memmap2::Mmap;
use tempfile::tempfile;

use super::fm_index::{RankBlock, BLOCK_SIZE};
use super::kmer_index::{KmerSeedTable, PosTable};

// Precomputed score table: score[mm] = 1.0 / (1.0 + mm)
const SCORE_LUT: [f32; 4] = [1.0, 0.5, 1.0 / 3.0, 0.25];

/// Reusable workspace to eliminate per-BWT-step heap allocations.
/// Allocated once per search pipeline, cleared and reused at every depth step.
pub(crate) struct StepWorkspace {
    // Non-dedup path
    pub qchar: Vec<u8>,
    pub c_l: Vec<u32>,
    pub c_r: Vec<u32>,
    pub c_bud: Vec<u8>,
    pub c_qid: Vec<u32>,
    pub c_str: Vec<u8>,
    pub c_valid: Vec<u8>,
    // Dedup path
    pub order: Vec<u32>,
    pub group_starts: Vec<usize>,
}

impl StepWorkspace {
    pub fn new() -> Self {
        StepWorkspace {
            qchar: Vec::new(),
            c_l: Vec::new(),
            c_r: Vec::new(),
            c_bud: Vec::new(),
            c_qid: Vec::new(),
            c_str: Vec::new(),
            c_valid: Vec::new(),
            order: Vec::new(),
            group_starts: Vec::new(),
        }
    }
}

pub const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

// ═══════════════════════════════════════════════════════════════════════════════
//  Hit Accumulator
// ═══════════════════════════════════════════════════════════════════════════════

pub struct HitAccumulator {
    pub query_id: Vec<u32>,
    pub position: Vec<u32>,
    pub strand: Vec<bool>,
    pub score: Vec<f32>,
}

impl HitAccumulator {
    pub fn new() -> Self {
        Self {
            query_id: Vec::new(),
            position: Vec::new(),
            strand: Vec::new(),
            score: Vec::new(),
        }
    }

    #[inline]
    pub fn push(&mut self, query_id: u32, position: u32, strand: bool, score: f32) {
        self.query_id.push(query_id);
        self.position.push(position);
        self.strand.push(strand);
        self.score.push(score);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Index abstraction
// ═══════════════════════════════════════════════════════════════════════════════

pub trait FmOcc {
    fn less(&self, c: u8) -> usize;
    fn occ(&self, pos: usize, c: u8) -> usize;
    fn sa(&self, idx: usize) -> usize;
    fn sa_len(&self) -> usize;

    /// Batch LF-mapping: compute all 4 `(nl, nr_exclusive)` pairs for A, C, G, T.
    /// Default: 8 scalar occ() calls. Override with interleaved rank for 2 cache lines.
    #[inline]
    fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        let mut result = [(0u32, 0u32); 4];
        for (bi, &base) in BASES.iter().enumerate() {
            let less_b = self.less(base);
            let occ_l = if l > 0 { self.occ(l - 1, base) } else { 0 };
            let occ_r = self.occ(r, base);
            result[bi] = ((less_b + occ_l) as u32, (less_b + occ_r) as u32);
        }
        result
    }

    /// Prefetch cache lines for a future `lf_map(l, r)` call.
    /// No-op by default; interleaved rank implementation issues 2 prefetches.
    #[inline]
    fn prefetch_lf(&self, _l: usize, _r: usize) {}

    /// Raw interleaved rank data for SIMD access (legacy flat layout).
    /// Returns `(data_slice, less_table)` where `data_slice[pos*4 + base_idx]` = occ count.
    /// `None` = no flat layout available (use `rank_blocks()` or scalar `lf_map` instead).
    fn rank_data(&self) -> Option<(&[u32], &[usize; 256])> { None }

    /// Cache-line-aligned rank blocks for SIMD access.
    /// Returns `(blocks, less_table)` where each `RankBlock` is 64 bytes covering
    /// 64 BWT positions. Preferred over `rank_data()` — 16× smaller, fits in L3.
    fn rank_blocks(&self) -> Option<(&[RankBlock], &[usize; 256])> { None }

    /// Raw suffix array as a contiguous slice, if available.
    /// Enables O(N) SA sweep for seed table construction.
    /// Default: None (use `sa(idx)` for individual lookups).
    fn sa_slice(&self) -> Option<&[usize]> { None }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Chromosome geometry
// ═══════════════════════════════════════════════════════════════════════════════

pub struct ChromGeometry {
    pub ranges: Vec<(usize, usize)>,
}

impl ChromGeometry {
    #[inline]
    pub fn is_valid(&self, pos: usize, query_len: usize) -> bool {
        let mut lo = 0usize;
        let mut hi = self.ranges.len();
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let (start, len) = self.ranges[mid];
            if pos < start {
                hi = mid;
            } else if pos >= start + len {
                lo = mid + 1;
            } else {
                return pos + query_len <= start + len;
            }
        }
        false
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Lineage Table — maps lineage_id → set of virtual query IDs
// ═══════════════════════════════════════════════════════════════════════════════

/// Flat-buffer lineage table. Each lineage entry is a contiguous slice of vqids
/// in `data`, addressed by `(offset, len)` in `spans`.
///
/// Virtual query ID encoding: `vqid = query_id * 2 + strand_bit`
///   - strand_bit 0 = forward (reads queries_fwd)
///   - strand_bit 1 = reverse complement (reads queries_rc)
pub struct LineageTable {
    data: Vec<u32>,
    spans: Vec<(u32, u32)>, // (offset, len) per lineage_id
}

impl LineageTable {
    pub fn new() -> Self {
        Self { data: Vec::new(), spans: Vec::new() }
    }

    pub fn with_capacity(n_entries: usize, n_data: usize) -> Self {
        Self {
            data: Vec::with_capacity(n_data),
            spans: Vec::with_capacity(n_entries),
        }
    }

    /// Allocate a new lineage entry. Returns the lineage_id.
    #[inline]
    pub fn alloc(&mut self, vqids: &[u32]) -> u32 {
        let offset = self.data.len() as u32;
        let len = vqids.len() as u32;
        self.data.extend_from_slice(vqids);
        let id = self.spans.len() as u32;
        self.spans.push((offset, len));
        id
    }

    /// Allocate from an iterator (avoids intermediate Vec for merges).
    pub fn alloc_iter(&mut self, iter: impl Iterator<Item = u32>, hint: usize) -> u32 {
        let offset = self.data.len() as u32;
        self.data.reserve(hint);
        let mut len = 0u32;
        for v in iter {
            self.data.push(v);
            len += 1;
        }
        let id = self.spans.len() as u32;
        self.spans.push((offset, len));
        id
    }

    /// Get the vqid slice for a lineage.
    #[inline]
    pub fn get(&self, id: u32) -> &[u32] {
        let (off, len) = self.spans[id as usize];
        &self.data[off as usize..(off + len) as usize]
    }

    #[inline]
    pub fn entry_len(&self, id: u32) -> usize {
        self.spans[id as usize].1 as usize
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Mmap-Backed SoA Frontier
// ═══════════════════════════════════════════════════════════════════════════════

/// One SoA column backed by a temporary file + mmap.
///
/// Write path: append via `File::write_all` into a temp file.
/// Read path: `mmap` the file read-only, cast to `&[T]`.
struct MmapColumn {
    file: File,
    len: usize, // number of elements (not bytes)
}

impl MmapColumn {
    fn new() -> io::Result<Self> {
        Ok(MmapColumn {
            file: tempfile()?,
            len: 0,
        })
    }

    /// Reset for reuse as a write target. Truncates the file to zero.
    fn reset(&mut self) -> io::Result<()> {
        self.file.set_len(0)?;
        self.file.seek(SeekFrom::Start(0))?;
        self.len = 0;
        Ok(())
    }

    /// Append a single value (little-endian bytes).
    #[inline]
    fn push_u32(&mut self, v: u32) -> io::Result<()> {
        self.file.write_all(&v.to_le_bytes())?;
        self.len += 1;
        Ok(())
    }

    #[inline]
    fn push_u8(&mut self, v: u8) -> io::Result<()> {
        self.file.write_all(&[v])?;
        self.len += 1;
        Ok(())
    }

    /// Flush and memory-map the file for reading.
    fn freeze_u32(&mut self) -> io::Result<MmapSliceU32> {
        self.file.flush()?;
        if self.len == 0 {
            return Ok(MmapSliceU32 { _mmap: None, len: 0 });
        }
        self.file.seek(SeekFrom::Start(0))?;
        let mmap = unsafe { Mmap::map(&self.file)? };
        Ok(MmapSliceU32 { _mmap: Some(mmap), len: self.len })
    }

    fn freeze_u8(&mut self) -> io::Result<MmapSliceU8> {
        self.file.flush()?;
        if self.len == 0 {
            return Ok(MmapSliceU8 { _mmap: None, len: 0 });
        }
        self.file.seek(SeekFrom::Start(0))?;
        let mmap = unsafe { Mmap::map(&self.file)? };
        Ok(MmapSliceU8 { _mmap: Some(mmap), len: self.len })
    }
}

/// Read-only mmap'd slice of u32 values.
struct MmapSliceU32 {
    _mmap: Option<Mmap>,
    len: usize,
}

impl MmapSliceU32 {
    #[inline]
    fn as_slice(&self) -> &[u32] {
        match &self._mmap {
            None => &[],
            Some(m) => unsafe {
                slice::from_raw_parts(m.as_ptr() as *const u32, self.len)
            },
        }
    }
}

/// Read-only mmap'd slice of u8 values.
struct MmapSliceU8 {
    _mmap: Option<Mmap>,
    len: usize,
}

impl MmapSliceU8 {
    #[inline]
    fn as_slice(&self) -> &[u8] {
        match &self._mmap {
            None => &[],
            Some(m) => &m[..self.len],
        }
    }
}

/// Mmap-backed SoA frontier. Write columns append to temp files;
/// read columns are frozen mmaps. Ping/pong by swapping roles.
pub struct MmapFrontier {
    col_l: MmapColumn,
    col_r: MmapColumn,
    col_budget: MmapColumn,
    col_qid: MmapColumn,
    col_strand: MmapColumn,
}

impl MmapFrontier {
    pub fn new() -> io::Result<Self> {
        Ok(MmapFrontier {
            col_l: MmapColumn::new()?,
            col_r: MmapColumn::new()?,
            col_budget: MmapColumn::new()?,
            col_qid: MmapColumn::new()?,
            col_strand: MmapColumn::new()?,
        })
    }

    /// Number of rows written so far (based on the `l` column).
    #[inline]
    pub fn len(&self) -> usize {
        self.col_l.len
    }

    /// Reset all columns for writing.
    pub fn reset(&mut self) -> io::Result<()> {
        self.col_l.reset()?;
        self.col_r.reset()?;
        self.col_budget.reset()?;
        self.col_qid.reset()?;
        self.col_strand.reset()?;
        Ok(())
    }

    /// Append one row.
    #[inline]
    pub fn push(&mut self, l: u32, r: u32, budget: u8, qid: u32, strand: u8) -> io::Result<()> {
        self.col_l.push_u32(l)?;
        self.col_r.push_u32(r)?;
        self.col_budget.push_u8(budget)?;
        self.col_qid.push_u32(qid)?;
        self.col_strand.push_u8(strand)?;
        Ok(())
    }

    /// Freeze all columns into read-only mmap slices.
    pub fn freeze(&mut self) -> io::Result<FrozenFrontier> {
        Ok(FrozenFrontier {
            l: self.col_l.freeze_u32()?,
            r: self.col_r.freeze_u32()?,
            budget: self.col_budget.freeze_u8()?,
            qid: self.col_qid.freeze_u32()?,
            strand: self.col_strand.freeze_u8()?,
        })
    }
}

/// Immutable snapshot of a frontier — all columns are mmap'd read-only slices.
pub struct FrozenFrontier {
    l: MmapSliceU32,
    r: MmapSliceU32,
    budget: MmapSliceU8,
    qid: MmapSliceU32,
    strand: MmapSliceU8,
}

impl FrozenFrontier {
    #[inline]
    pub fn len(&self) -> usize {
        self.l.len
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  In-memory frontier (small batches bypass mmap overhead)
// ═══════════════════════════════════════════════════════════════════════════════

/// Vec-backed SoA frontier for small batches where mmap overhead isn't justified.
#[repr(align(32))]
pub struct SearchFrontier {
    pub l: Vec<u32>,
    pub r: Vec<u32>,
    pub budget: Vec<u8>,
    pub query_id: Vec<u32>,
    pub strand: Vec<u8>,
}

impl SearchFrontier {
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            l: Vec::with_capacity(cap),
            r: Vec::with_capacity(cap),
            budget: Vec::with_capacity(cap),
            query_id: Vec::with_capacity(cap),
            strand: Vec::with_capacity(cap),
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.l.len()
    }

    #[inline]
    pub fn clear(&mut self) {
        self.l.clear();
        self.r.clear();
        self.budget.clear();
        self.query_id.clear();
        self.strand.clear();
    }

    #[inline]
    pub fn push(&mut self, l: u32, r: u32, budget: u8, query_id: u32, strand: u8) {
        self.l.push(l);
        self.r.push(r);
        self.budget.push(budget);
        self.query_id.push(query_id);
        self.strand.push(strand);
    }

    pub fn seed(n_queries: usize, sa_len: u32, max_mm: u8) -> Self {
        let n = n_queries * 2;
        let mut f = Self::with_capacity(n);
        for qid in 0..n_queries as u32 {
            f.push(0, sa_len - 1, max_mm, qid, 1);
            f.push(0, sa_len - 1, max_mm, qid, 0);
        }
        f
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Branchless step_depth — streams from frozen mmap into MmapFrontier writer
// ═══════════════════════════════════════════════════════════════════════════════

/// Process one depth level in streaming chunks from a frozen (read) frontier
/// into an mmap writer (next) frontier. Terminal hits go to the accumulator.
///
/// Reads from `frozen` (mmap'd read-only), writes survivors to `next_writer`.
/// Chunk size bounds peak RAM for scratch buffers.
const CHUNK_SIZE: usize = 1 << 16; // 64K rows per chunk

pub fn step_depth_mmap<I: FmOcc>(
    index: &I,
    frozen: &FrozenFrontier,
    depth: usize,
    query_len: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
    next_writer: &mut MmapFrontier,
) -> io::Result<()> {
    let n = frozen.len();
    if n == 0 {
        return Ok(());
    }

    let fl = frozen.l.as_slice();
    let fr = frozen.r.as_slice();
    let fb = frozen.budget.as_slice();
    let fq = frozen.qid.as_slice();
    let fs = frozen.strand.as_slice();

    let col = query_len - 1 - depth;
    let is_terminal = depth + 1 == query_len;

    // Process in fixed-size chunks to bound scratch memory.
    let mut chunk_start = 0;
    while chunk_start < n {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(n);
        let chunk_n = chunk_end - chunk_start;

        // Precompute qchar for this chunk.
        let mut qchar = Vec::with_capacity(chunk_n);
        for i in chunk_start..chunk_end {
            let qid = fq[i] as usize;
            qchar.push(if fs[i] != 0 {
                queries_fwd[qid][col]
            } else {
                queries_rc[qid][col]
            });
        }

        // Scratch buffers for 4 × chunk_n candidates.
        let cap = 4 * chunk_n;
        let mut c_l = Vec::with_capacity(cap);
        let mut c_r = Vec::with_capacity(cap);
        let mut c_bud = Vec::with_capacity(cap);
        let mut c_qid = Vec::with_capacity(cap);
        let mut c_str = Vec::with_capacity(cap);
        let mut c_valid = Vec::with_capacity(cap);

        // ── Phase 1: EXPAND (branch-free) ───────────────────────────────
        for &base in &BASES {
            let less_b = index.less(base);

            for ci in 0..chunk_n {
                let gi = chunk_start + ci; // global index
                let l = fl[gi] as usize;
                let r = fr[gi] as usize;

                let occ_l = if l > 0 { index.occ(l - 1, base) } else { 0 };
                let occ_r = index.occ(r, base);

                let nl = (less_b + occ_l) as u32;
                let nr_exc = (less_b + occ_r) as u32;

                let mm_cost = (base != qchar[ci]) as u8;
                let new_bud = fb[gi].saturating_sub(mm_cost);
                let valid = (nl < nr_exc) & (fb[gi] >= mm_cost);

                c_l.push(nl);
                c_r.push(nr_exc.wrapping_sub(1));
                c_bud.push(new_bud);
                c_qid.push(fq[gi]);
                c_str.push(fs[gi]);
                c_valid.push(valid as u8);
            }
        }

        // ── Phase 2: COMPRESS + ACCUMULATE ──────────────────────────────
        for j in 0..cap {
            if c_valid[j] == 0 {
                continue;
            }

            if is_terminal {
                let mm_used = max_mm - c_bud[j];
                let score = SCORE_LUT[mm_used as usize];
                let lo = c_l[j] as usize;
                let hi = c_r[j] as usize + 1;
                let fwd = c_str[j] != 0;

                for sa_idx in lo..hi {
                    let sa_pos = index.sa(sa_idx);
                    if chroms.is_valid(sa_pos, query_len) {
                        hits.push(c_qid[j], sa_pos as u32, fwd, score);
                    }
                }
            } else {
                next_writer.push(c_l[j], c_r[j], c_bud[j], c_qid[j], c_str[j])?;
            }
        }

        chunk_start = chunk_end;
    }

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
//  In-memory step_depth (unchanged — used for small-batch fast path)
// ═══════════════════════════════════════════════════════════════════════════════

pub fn step_depth<I: FmOcc>(
    index: &I,
    frontier: &SearchFrontier,
    depth: usize,
    query_len: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    max_mm: u8,
    max_width: u32,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
    next: &mut SearchFrontier,
    ws: &mut StepWorkspace,
) {
    let n = frontier.len();
    if n == 0 {
        return;
    }
    next.clear();

    let col = query_len - 1 - depth;

    // Take workspace buffers (zero-alloc after first step — retains capacity).
    let mut qchar = std::mem::take(&mut ws.qchar);
    let mut c_l = std::mem::take(&mut ws.c_l);
    let mut c_r = std::mem::take(&mut ws.c_r);
    let mut c_bud = std::mem::take(&mut ws.c_bud);
    let mut c_qid = std::mem::take(&mut ws.c_qid);
    let mut c_str = std::mem::take(&mut ws.c_str);
    let mut c_valid = std::mem::take(&mut ws.c_valid);

    qchar.clear();
    for i in 0..n {
        let qid = frontier.query_id[i] as usize;
        qchar.push(if frontier.strand[i] != 0 {
            queries_fwd[qid][col]
        } else {
            queries_rc[qid][col]
        });
    }

    let cap = 4 * n;
    c_l.clear();
    c_r.clear();
    c_bud.clear();
    c_qid.clear();
    c_str.clear();
    c_valid.clear();

    for &base in &BASES {
        let less_b = index.less(base);
        for i in 0..n {
            let l = frontier.l[i] as usize;
            let r = frontier.r[i] as usize;
            let occ_l = if l > 0 { index.occ(l - 1, base) } else { 0 };
            let occ_r = index.occ(r, base);
            let nl = (less_b + occ_l) as u32;
            let nr_exc = (less_b + occ_r) as u32;
            let mm_cost = (base != qchar[i]) as u8;
            let new_bud = frontier.budget[i].saturating_sub(mm_cost);
            let width = nr_exc.wrapping_sub(nl);
            // Prune mismatch branches in repetitive regions (width > cap).
            // Exact-match branches (cost=0) always survive regardless of width.
            let valid = (nl < nr_exc)
                & (frontier.budget[i] >= mm_cost)
                & (mm_cost == 0 || width <= max_width);
            c_l.push(nl);
            c_r.push(nr_exc.wrapping_sub(1));
            c_bud.push(new_bud);
            c_qid.push(frontier.query_id[i]);
            c_str.push(frontier.strand[i]);
            c_valid.push(valid as u8);
        }
    }

    let is_terminal = depth + 1 == query_len;
    for j in 0..cap {
        if c_valid[j] == 0 {
            continue;
        }
        if is_terminal {
            let mm_used = max_mm - c_bud[j];
            let score = SCORE_LUT[mm_used as usize];
            let lo = c_l[j] as usize;
            let hi = c_r[j] as usize + 1;
            let fwd = c_str[j] != 0;
            for sa_idx in lo..hi {
                let sa_pos = index.sa(sa_idx);
                if chroms.is_valid(sa_pos, query_len) {
                    hits.push(c_qid[j], sa_pos as u32, fwd, score);
                }
            }
        } else {
            next.push(c_l[j], c_r[j], c_bud[j], c_qid[j], c_str[j]);
        }
    }

    // Return buffers to workspace (capacity preserved for next step).
    ws.qchar = qchar;
    ws.c_l = c_l;
    ws.c_r = c_r;
    ws.c_bud = c_bud;
    ws.c_qid = c_qid;
    ws.c_str = c_str;
    ws.c_valid = c_valid;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Sort-Merge Frontier Deduplication
// ═══════════════════════════════════════════════════════════════════════════════

/// Sort the frontier by `(l, r, budget)` and merge identical rows.
/// When merging, union their lineage vqid sets in the lineage table.
fn dedup_frontier(frontier: &mut SearchFrontier, lineage: &mut LineageTable) {
    let n = frontier.len();
    if n <= 1 {
        return;
    }

    // Sort an index array by (l, r, budget).
    let mut idx: Vec<u32> = (0..n as u32).collect();
    idx.sort_unstable_by(|&a, &b| {
        let a = a as usize;
        let b = b as usize;
        frontier.l[a].cmp(&frontier.l[b])
            .then(frontier.r[a].cmp(&frontier.r[b]))
            .then(frontier.budget[a].cmp(&frontier.budget[b]))
    });

    let mut new_l = Vec::with_capacity(n);
    let mut new_r = Vec::with_capacity(n);
    let mut new_bud = Vec::with_capacity(n);
    let mut new_lin = Vec::with_capacity(n);
    let mut new_str = Vec::with_capacity(n);

    let mut i = 0usize;
    while i < n {
        let ii = idx[i] as usize;
        let l = frontier.l[ii];
        let r = frontier.r[ii];
        let b = frontier.budget[ii];

        // Find run of identical (l, r, budget).
        let mut j = i + 1;
        while j < n {
            let jj = idx[j] as usize;
            if frontier.l[jj] != l || frontier.r[jj] != r || frontier.budget[jj] != b {
                break;
            }
            j += 1;
        }

        if j == i + 1 {
            // Single element — no merge needed.
            new_l.push(l);
            new_r.push(r);
            new_bud.push(b);
            new_lin.push(frontier.query_id[ii]); // lineage_id stored in query_id field
            new_str.push(0);
        } else {
            // Merge: union all lineages into a temp buffer, then alloc.
            let mut merged_vqids: Vec<u32> = Vec::new();
            for k in i..j {
                let lin = frontier.query_id[idx[k] as usize];
                merged_vqids.extend_from_slice(lineage.get(lin));
            }
            let merged_id = lineage.alloc(&merged_vqids);
            new_l.push(l);
            new_r.push(r);
            new_bud.push(b);
            new_lin.push(merged_id);
            new_str.push(0);
        }
        i = j;
    }

    frontier.l = new_l;
    frontier.r = new_r;
    frontier.budget = new_bud;
    frontier.query_id = new_lin;
    frontier.strand = new_str;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Block-Based SIMD Rank Helper
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute rank for all 4 bases at `pos` using SSE2 and cache-line-aligned blocks.
///
/// Returns `__m128i` with `[rank_A, rank_C, rank_G, rank_T]` as i32 lanes.
/// Each block is 64 bytes (1 cache line): a single memory fetch provides all data
/// needed for the rank query. The popcount is computed via 4 scalar POPCNT
/// instructions (~1 cycle each), negligible vs the ~90ns DRAM savings.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn rank_all_blocks_sse(blocks: &[RankBlock], pos: usize) -> std::arch::x86_64::__m128i {
    use std::arch::x86_64::*;

    let block = &blocks[pos / BLOCK_SIZE];
    let offset = pos % BLOCK_SIZE;

    // Select counts and bitvector half based on offset within block.
    let (base_counts, bv) = if offset < 32 {
        (
            _mm_load_si128(block.counts.as_ptr() as *const __m128i),
            _mm_load_si128(block.bv_lo.as_ptr() as *const __m128i),
        )
    } else {
        (
            _mm_load_si128(block.mid.as_ptr() as *const __m128i),
            _mm_load_si128(block.bv_hi.as_ptr() as *const __m128i),
        )
    };

    // Inclusive mask: bits 0..=(offset % 32) are set.
    let shift = (offset % 32) + 1;
    let mask_val = if shift >= 32 { u32::MAX } else { (1u32 << shift) - 1 };
    let mask = _mm_set1_epi32(mask_val as i32);
    let masked = _mm_and_si128(bv, mask);

    // Popcount per 32-bit lane via scalar POPCNT.
    let mut bv_arr = [0u32; 4];
    _mm_storeu_si128(bv_arr.as_mut_ptr() as *mut __m128i, masked);
    let pop = _mm_set_epi32(
        bv_arr[3].count_ones() as i32,
        bv_arr[2].count_ones() as i32,
        bv_arr[1].count_ones() as i32,
        bv_arr[0].count_ones() as i32,
    );

    _mm_add_epi32(base_counts, pop)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  AVX2 Wavefront Primitives — 8-wide vectorized rank via gather
// ═══════════════════════════════════════════════════════════════════════════════

/// Vectorized 32-bit lane popcount using nibble-lookup.
///
/// Each 32-bit lane's population count is computed entirely in SIMD:
///   1. Split each byte into low/high nibbles
///   2. Use `vpshufb` as a 4-bit → popcount LUT (0→0, 1→1, ..., F→4)
///   3. Sum byte popcounts within each 32-bit lane via `vpmaddubsw` + `vpmaddwd`
///
/// Zero branches, zero scalar extraction. ~3 cycles throughput.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn avx2_popcount_epi32(v: std::arch::x86_64::__m256i) -> std::arch::x86_64::__m256i {
    use std::arch::x86_64::*;

    let low_mask = _mm256_set1_epi8(0x0F);
    let lookup = _mm256_setr_epi8(
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, // lane 0
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, // lane 1
    );
    let lo = _mm256_shuffle_epi8(lookup, _mm256_and_si256(v, low_mask));
    let hi = _mm256_shuffle_epi8(
        lookup,
        _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask),
    );
    let byte_pop = _mm256_add_epi8(lo, hi);
    // Reduce bytes→32-bit: maddubs pairs u8→i16, madd pairs i16→i32
    let word_pop = _mm256_maddubs_epi16(byte_pop, _mm256_set1_epi8(1));
    _mm256_madd_epi16(word_pop, _mm256_set1_epi16(1))
}

/// 8-position vectorized rank computation using AVX2 gathers.
///
/// For each of 4 bases, computes `rank(pos)` = count of that base in `BWT[0..=pos]`
/// for all 8 positions simultaneously. Uses `_mm256_i32gather_epi32` to load
/// counts and bitvectors from the cache-line-aligned `RankBlock` array.
///
/// **Memory layout**: Each `RankBlock` is 16 × i32 (64 bytes):
/// ```text
/// [0..3]   counts[A,C,G,T]    — absolute cumulative before block start
/// [4..7]   bv_lo[A,C,G,T]     — bitvectors for positions 0..31
/// [8..11]  bv_hi[A,C,G,T]     — bitvectors for positions 32..63
/// [12..15] mid[A,C,G,T]       — absolute cumulative at block_start + 32
/// ```
///
/// **Sentinel handling**: positions equal to `u32::MAX` (the `l=0` case where
/// `rank(l-1)` should be zero) produce all-zero results.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn rank_8way_avx2(
    blocks_ptr: *const i32,
    positions: [u32; 8],
) -> [[u32; 4]; 8] {
    use std::arch::x86_64::*;

    let pos_v = _mm256_loadu_si256(positions.as_ptr() as *const __m256i);

    // Sentinel mask: lanes where pos == u32::MAX (l=0 case)
    let sentinel = _mm256_cmpeq_epi32(pos_v, _mm256_set1_epi32(-1i32));
    let valid_mask = _mm256_andnot_si256(sentinel, _mm256_set1_epi32(-1i32));

    // Clamp sentinels to 0 to avoid wild gather addresses
    let pos_safe = _mm256_andnot_si256(sentinel, pos_v);

    // block_idx = pos >> 6
    let block_idx = _mm256_srli_epi32(pos_safe, 6);
    // offset = pos & 63
    let offset = _mm256_and_si256(pos_safe, _mm256_set1_epi32(63));
    // is_hi: offset >= 32 → all-ones mask per lane
    let is_hi = _mm256_cmpgt_epi32(offset, _mm256_set1_epi32(31));

    // off_in_half = offset & 31 (0..31 within the lo or hi bitvector)
    let off_in_half = _mm256_and_si256(offset, _mm256_set1_epi32(31));
    // Inclusive mask: (1 << (off_in_half + 1)) - 1
    // When off_in_half = 31, shift = 32, sllv produces 0, then 0 - 1 = 0xFFFFFFFF. Correct.
    let one = _mm256_set1_epi32(1);
    let shift = _mm256_add_epi32(off_in_half, one);
    let inc_mask = _mm256_sub_epi32(_mm256_sllv_epi32(one, shift), one);

    // base_offset = block_idx * 16 (16 i32s per RankBlock)
    let base_off = _mm256_slli_epi32(block_idx, 4);

    let mut result = [[0u32; 4]; 8];

    for b in 0..4u32 {
        let b_i32 = b as i32;

        // lo path: count at base_off + b, bv at base_off + 4 + b
        let count_idx_lo = _mm256_add_epi32(base_off, _mm256_set1_epi32(b_i32));
        let bv_idx_lo = _mm256_add_epi32(base_off, _mm256_set1_epi32(4 + b_i32));

        // hi path: count at base_off + 12 + b (mid), bv at base_off + 8 + b (bv_hi)
        let count_idx_hi = _mm256_add_epi32(base_off, _mm256_set1_epi32(12 + b_i32));
        let bv_idx_hi = _mm256_add_epi32(base_off, _mm256_set1_epi32(8 + b_i32));

        // Blend: select hi path where offset >= 32
        let count_idx = _mm256_blendv_epi8(count_idx_lo, count_idx_hi, is_hi);
        let bv_idx = _mm256_blendv_epi8(bv_idx_lo, bv_idx_hi, is_hi);

        // Gather counts and bitvectors (2 gathers per base)
        let counts = _mm256_i32gather_epi32::<4>(blocks_ptr, count_idx);
        let bvs = _mm256_i32gather_epi32::<4>(blocks_ptr, bv_idx);

        // Mask bitvectors, popcount, accumulate
        let masked_bvs = _mm256_and_si256(bvs, inc_mask);
        let popcounts = avx2_popcount_epi32(masked_bvs);
        let rank = _mm256_add_epi32(counts, popcounts);

        // Zero out sentinel lanes
        let rank_masked = _mm256_and_si256(rank, valid_mask);

        // Extract to result array
        let mut rank_arr = [0u32; 8];
        _mm256_storeu_si256(rank_arr.as_mut_ptr() as *mut __m256i, rank_masked);
        for i in 0..8 {
            result[i][b as usize] = rank_arr[i];
        }
    }

    result
}

/// 16-position interleaved rank computation using AVX2 gathers.
///
/// Computes ranks for two independent batches of 8 positions simultaneously,
/// interleaving gather instructions across the l-batch and r-batch so that the
/// 11-12 cycle gather latency of one batch overlaps with the arithmetic of the
/// other. This saturates Port 5 (shuffle) and Port 0/1 (ALU) simultaneously,
/// improving throughput over two sequential `rank_8way_avx2` calls.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn rank_16way_interleaved(
    blocks_ptr: *const i32,
    positions_a: [u32; 8],
    positions_b: [u32; 8],
) -> ([[u32; 4]; 8], [[u32; 4]; 8]) {
    use std::arch::x86_64::*;

    let neg1 = _mm256_set1_epi32(-1i32);
    let one = _mm256_set1_epi32(1);
    let mask63 = _mm256_set1_epi32(63);
    let mask31 = _mm256_set1_epi32(31);
    let thresh31 = _mm256_set1_epi32(31);

    // ── Setup batch A (l positions) ──────────────────────────────────────
    let pos_a = _mm256_loadu_si256(positions_a.as_ptr() as *const __m256i);
    let sentinel_a = _mm256_cmpeq_epi32(pos_a, neg1);
    let valid_mask_a = _mm256_andnot_si256(sentinel_a, neg1);
    let pos_safe_a = _mm256_andnot_si256(sentinel_a, pos_a);
    let block_idx_a = _mm256_srli_epi32(pos_safe_a, 6);
    let offset_a = _mm256_and_si256(pos_safe_a, mask63);
    let is_hi_a = _mm256_cmpgt_epi32(offset_a, thresh31);
    let off_in_half_a = _mm256_and_si256(offset_a, mask31);
    let inc_mask_a = _mm256_sub_epi32(
        _mm256_sllv_epi32(one, _mm256_add_epi32(off_in_half_a, one)),
        one,
    );
    let base_off_a = _mm256_slli_epi32(block_idx_a, 4);

    // ── Setup batch B (r positions) ──────────────────────────────────────
    let pos_b = _mm256_loadu_si256(positions_b.as_ptr() as *const __m256i);
    let sentinel_b = _mm256_cmpeq_epi32(pos_b, neg1);
    let valid_mask_b = _mm256_andnot_si256(sentinel_b, neg1);
    let pos_safe_b = _mm256_andnot_si256(sentinel_b, pos_b);
    let block_idx_b = _mm256_srli_epi32(pos_safe_b, 6);
    let offset_b = _mm256_and_si256(pos_safe_b, mask63);
    let is_hi_b = _mm256_cmpgt_epi32(offset_b, thresh31);
    let off_in_half_b = _mm256_and_si256(offset_b, mask31);
    let inc_mask_b = _mm256_sub_epi32(
        _mm256_sllv_epi32(one, _mm256_add_epi32(off_in_half_b, one)),
        one,
    );
    let base_off_b = _mm256_slli_epi32(block_idx_b, 4);

    let mut result_a = [[0u32; 4]; 8];
    let mut result_b = [[0u32; 4]; 8];

    for b in 0..4u32 {
        let b_i32 = b as i32;

        // ── 1. ADDRESS COMPUTATION — both batches (independent) ──────────
        // Batch A: lo/hi index selection
        let count_idx_a = _mm256_blendv_epi8(
            _mm256_add_epi32(base_off_a, _mm256_set1_epi32(b_i32)),       // lo: base_off + b
            _mm256_add_epi32(base_off_a, _mm256_set1_epi32(12 + b_i32)), // hi: base_off + 12 + b
            is_hi_a,
        );
        let bv_idx_a = _mm256_blendv_epi8(
            _mm256_add_epi32(base_off_a, _mm256_set1_epi32(4 + b_i32)),  // lo: base_off + 4 + b
            _mm256_add_epi32(base_off_a, _mm256_set1_epi32(8 + b_i32)),  // hi: base_off + 8 + b
            is_hi_a,
        );

        // Batch B: lo/hi index selection
        let count_idx_b = _mm256_blendv_epi8(
            _mm256_add_epi32(base_off_b, _mm256_set1_epi32(b_i32)),
            _mm256_add_epi32(base_off_b, _mm256_set1_epi32(12 + b_i32)),
            is_hi_b,
        );
        let bv_idx_b = _mm256_blendv_epi8(
            _mm256_add_epi32(base_off_b, _mm256_set1_epi32(4 + b_i32)),
            _mm256_add_epi32(base_off_b, _mm256_set1_epi32(8 + b_i32)),
            is_hi_b,
        );

        // ── 2. INTERLEAVED GATHERS — overlap gather latency ─────────────
        let counts_a = _mm256_i32gather_epi32::<4>(blocks_ptr, count_idx_a); // A count gather fires
        let counts_b = _mm256_i32gather_epi32::<4>(blocks_ptr, count_idx_b); // B dispatches while A waits
        let bvs_a = _mm256_i32gather_epi32::<4>(blocks_ptr, bv_idx_a);      // A's counts likely ready
        let bvs_b = _mm256_i32gather_epi32::<4>(blocks_ptr, bv_idx_b);      // overlap continues

        // ── 3. INTERLEAVED MATH — arithmetic while late gathers land ────
        let masked_a = _mm256_and_si256(bvs_a, inc_mask_a);
        let pop_a = avx2_popcount_epi32(masked_a);       // Port 5 shuffle chain
        let masked_b = _mm256_and_si256(bvs_b, inc_mask_b);
        let pop_b = avx2_popcount_epi32(masked_b);       // Port 0/1 ALU overlaps

        // ── 4. ACCUMULATE + sentinel mask ────────────────────────────────
        let rank_a = _mm256_and_si256(_mm256_add_epi32(counts_a, pop_a), valid_mask_a);
        let rank_b = _mm256_and_si256(_mm256_add_epi32(counts_b, pop_b), valid_mask_b);

        // ── 5. EXTRACT — interleaved stores to overlap SLF stalls ───────
        let mut arr_a = [0u32; 8];
        let mut arr_b = [0u32; 8];
        _mm256_storeu_si256(arr_a.as_mut_ptr() as *mut __m256i, rank_a);
        _mm256_storeu_si256(arr_b.as_mut_ptr() as *mut __m256i, rank_b);
        for i in 0..8 {
            result_a[i][b as usize] = arr_a[i];
            result_b[i][b as usize] = arr_b[i];
        }
    }

    (result_a, result_b)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Vectorized Pigeonhole Filter + SIMD Hamming Verification
// ═══════════════════════════════════════════════════════════════════════════════

/// Exact-match pigeonhole filter: splits each query into `max_mm + 1` segments
/// and performs a 0-mismatch BWT backward search on each segment using AVX2.
///
/// By the pigeonhole principle, any alignment with ≤ max_mm mismatches must have
/// at least one segment that matches exactly. Seeds with non-empty BWT intervals
/// (survivors) are returned for SA resolution and Hamming verification.
///
/// Returns `(surv_l, surv_r, surv_qid, surv_seg, surv_strand)`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn exact_pigeonhole_filter(
    blocks_ptr: *const i32,
    less_arr: &[u32; 4],
    bwt_len: u32,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u8>) {
    use std::arch::x86_64::*;

    let n_queries = queries_fwd.len();
    let n_seg = max_mm as usize + 1;
    let n_qs = n_queries * 2; // number of query-strand pairs

    // Dynamic segment partition.
    let base_len = query_len / n_seg;
    let remainder = query_len % n_seg;
    let mut seg_bounds: Vec<(usize, usize)> = Vec::with_capacity(n_seg);
    {
        let mut pos = 0;
        for s in 0..n_seg {
            let slen = base_len + if s < remainder { 1 } else { 0 };
            seg_bounds.push((pos, pos + slen));
            pos += slen;
        }
    }

    // SoA flat arrays. Layout: seed_idx = seg_id * n_qs + (qid * 2 + strand)
    // Seeds for the same segment are contiguous → each segment's loop touches
    // only its own seeds.
    let n_seeds = n_seg * n_qs;
    let mut seed_l = vec![0u32; n_seeds];
    let mut seed_r = vec![bwt_len; n_seeds];

    // Process segment by segment — only touch this segment's contiguous block.
    for seg_id in 0..n_seg {
        let (seg_start, seg_end) = seg_bounds[seg_id];
        let seg_len = seg_end - seg_start;
        let seg_offset = seg_id * n_qs; // start of this segment's seeds

        for depth in 0..seg_len {
            let col = seg_end - 1 - depth; // BWT backward: right-to-left

            // Process this segment's seeds in batches of 8 for AVX2.
            let mut batch_start = 0;
            while batch_start < n_qs {
                let batch_end = (batch_start + 8).min(n_qs);
                let batch_len = batch_end - batch_start;

                // Collect (l-1, r) positions for rank computation.
                let mut pos_l = [u32::MAX; 8];
                let mut pos_r = [0u32; 8];
                for k in 0..batch_len {
                    let si = seg_offset + batch_start + k;
                    let l = seed_l[si];
                    let r = seed_r[si];
                    pos_l[k] = if l > 0 { l - 1 } else { u32::MAX };
                    pos_r[k] = if r > 0 { r - 1 } else { 0 }; // exclusive→inclusive for rank
                }

                // 16-wide interleaved rank.
                let (ranks_l, ranks_r) = rank_16way_interleaved(blocks_ptr, pos_l, pos_r);

                // Per-seed base selection and interval update.
                for k in 0..batch_len {
                    let si = seg_offset + batch_start + k;
                    let l = seed_l[si];
                    let r = seed_r[si];

                    // Dead seed check: l >= r means empty interval.
                    if l >= r {
                        continue;
                    }

                    // Decode seed index → qid, strand.
                    let qs_idx = batch_start + k; // = si - seg_offset
                    let qid = qs_idx / 2;
                    let strand = qs_idx % 2; // 0 = fwd, 1 = rc

                    let c = if strand == 0 {
                        queries_fwd[qid][col]
                    } else {
                        queries_rc[qid][col]
                    };

                    let bi = match c {
                        b'A' => 0usize,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => { seed_l[si] = bwt_len; seed_r[si] = 0; continue; }
                    };

                    let nl = less_arr[bi] + ranks_l[k][bi];
                    let nr = less_arr[bi] + ranks_r[k][bi];

                    if nl >= nr {
                        seed_l[si] = bwt_len;
                        seed_r[si] = 0;
                    } else {
                        seed_l[si] = nl;
                        seed_r[si] = nr; // exclusive
                    }
                }

                batch_start = batch_end;
            }
        }
    }

    // AVX2 existence check + sparse collection with width cap.
    // Seeds wider than MAX_SEED_WIDTH are in repetitive regions and discarded.
    const MAX_SEED_WIDTH: u32 = 128;

    let mut surv_l = Vec::new();
    let mut surv_r = Vec::new();
    let mut surv_qid = Vec::new();
    let mut surv_seg = Vec::new();
    let mut surv_strand = Vec::new();

    let mut batch_start = 0;
    while batch_start < n_seeds {
        let batch_end = (batch_start + 8).min(n_seeds);
        let batch_len = batch_end - batch_start;

        let mut l_buf = [0u32; 8];
        let mut r_buf = [0u32; 8];
        for k in 0..batch_len {
            l_buf[k] = seed_l[batch_start + k];
            r_buf[k] = seed_r[batch_start + k];
        }

        let l_v = _mm256_loadu_si256(l_buf.as_ptr() as *const __m256i);
        let r_v = _mm256_loadu_si256(r_buf.as_ptr() as *const __m256i);
        let alive_v = _mm256_cmpgt_epi32(
            _mm256_xor_si256(r_v, _mm256_set1_epi32(i32::MIN)),
            _mm256_xor_si256(l_v, _mm256_set1_epi32(i32::MIN)),
        );
        let mask = _mm256_movemask_ps(_mm256_castsi256_ps(alive_v)) as u32;

        let mut m = mask & ((1u32 << batch_len) - 1);
        while m != 0 {
            let bit = m.trailing_zeros() as usize;
            let si = batch_start + bit;
            let l = seed_l[si];
            let r = seed_r[si];
            let width = r - l;

            // Width cap: discard seeds in repetitive regions.
            if width <= MAX_SEED_WIDTH {
                let seg_id = si / n_qs;
                let qs_idx = si % n_qs;
                let qid = qs_idx / 2;
                let strand = (qs_idx % 2) as u8;

                surv_l.push(l);
                surv_r.push(r);
                surv_qid.push(qid as u32);
                surv_seg.push(seg_id as u32);
                surv_strand.push(strand);
            }

            m &= m - 1;
        }

        batch_start = batch_end;
    }

    (surv_l, surv_r, surv_qid, surv_seg, surv_strand)
}

/// SIMD Hamming distance for queries ≤ 32 bytes.
///
/// Compares `query[0..query_len]` against `genome[0..query_len]` using a single
/// AVX2 compare + popcount. Returns the number of mismatches.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn simd_hamming_le32(query: &[u8], genome: &[u8], query_len: usize) -> u32 {
    use std::arch::x86_64::*;

    let mut q_buf = [0u8; 32];
    let mut g_buf = [0u8; 32];
    q_buf[..query_len].copy_from_slice(&query[..query_len]);
    g_buf[..query_len].copy_from_slice(&genome[..query_len]);

    let q_v = _mm256_loadu_si256(q_buf.as_ptr() as *const __m256i);
    let g_v = _mm256_loadu_si256(g_buf.as_ptr() as *const __m256i);
    let eq_v = _mm256_cmpeq_epi8(q_v, g_v);
    let eq_bits = _mm256_movemask_epi8(eq_v) as u32;

    let valid = if query_len >= 32 { u32::MAX } else { (1u32 << query_len) - 1 };
    let matches = (eq_bits & valid).count_ones();
    query_len as u32 - matches
}

/// Hamming distance for queries of any length, with early exit.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn simd_hamming(query: &[u8], genome: &[u8], query_len: usize, max_mm: u8) -> u32 {
    if query_len <= 32 {
        return simd_hamming_le32(query, genome, query_len);
    }
    // Process in 32-byte chunks with early exit.
    let mut total_mm = 0u32;
    let mut offset = 0;
    while offset < query_len {
        let chunk = (query_len - offset).min(32);
        let mm = simd_hamming_le32(&query[offset..], &genome[offset..], chunk);
        total_mm += mm;
        if total_mm > max_mm as u32 {
            return total_mm;
        }
        offset += 32;
    }
    total_mm
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Unified Cost-Based Depth Step with Lineage Dedup
// ═══════════════════════════════════════════════════════════════════════════════

/// Map DNA base to BASES index: A=0, C=1, G=2, T=3.
#[inline]
fn base_char_to_idx(b: u8) -> usize {
    match b {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 0,
    }
}

/// Flattened work item for unified cost computation within an (l,r) group.
/// Pre-resolves the target character so the ACGT loop is pure integer math.
struct WorkItem {
    vqid: u32,
    budget: u8,
    target_char: u8,
}

/// Cardinality clog threshold: child intervals wider than this with remaining
/// mismatch budget are pruned. Prevents frontier explosion from repetitive
/// elements. Budget-0 branches (exact-match-only) are always kept.
///
const CLOG_THRESHOLD: u32 = 500;

/// Dispatch per-base survivors: resolve terminal SA positions or bucket into next frontier.
/// Applies cost-aware interval width pruning: only MISMATCH branches (cost=1) in
/// child intervals wider than `CLOG_THRESHOLD` are dropped. Exact-match branches
/// (cost=0) always pass through, even in wide intervals with remaining budget.
///
/// When a mismatch branch is clog-pruned, the affected vqids are recorded in
/// `clog_vqids` so a targeted verify pass can recover those hits.
#[inline]
fn dispatch_survivors<I: FmOcc>(
    base_survivors: &[Vec<(u32, u8, u8)>; 4],
    nl_arr: &[u32; 4],
    nr_arr: &[u32; 4],
    is_terminal: bool,
    query_len: usize,
    max_mm: u8,
    index: &I,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
    next: &mut SearchFrontier,
    lineage: &mut LineageTable,
    bkt: &mut [Vec<u32>; 4],
    clog_vqids: &mut Vec<u32>,
) {
    for bi in 0..4 {
        if base_survivors[bi].is_empty() { continue; }
        let nl = nl_arr[bi];
        let nr_exc = nr_arr[bi];
        if nl >= nr_exc { continue; }

        if is_terminal {
            for &(vqid, new_bud, _cost) in &base_survivors[bi] {
                let mm_used = max_mm - new_bud;
                let score = SCORE_LUT[mm_used as usize];
                let fwd = vqid & 1 == 0;
                for sa_idx in nl as usize..nr_exc as usize {
                    let sa_pos = index.sa(sa_idx);
                    if chroms.is_valid(sa_pos, query_len) {
                        hits.push(vqid >> 1, sa_pos as u32, fwd, score);
                    }
                }
            }
        } else {
            let child_width = nr_exc - nl;
            let clogged = child_width > CLOG_THRESHOLD;

            for b in bkt.iter_mut() { b.clear(); }
            for &(vqid, new_bud, cost) in &base_survivors[bi] {
                // Cost-aware clog: only prune MISMATCH branches (cost=1) in
                // wide intervals. Exact-match branches (cost=0) always survive.
                if clogged && cost > 0 {
                    clog_vqids.push(vqid);
                    continue;
                }
                bkt[new_bud as usize].push(vqid);
            }
            let nr = nr_exc - 1;
            for (bud_val, bucket) in bkt.iter().enumerate() {
                if bucket.is_empty() { continue; }
                let lin_id = lineage.alloc(bucket);
                next.push(nl, nr, bud_val as u8, lin_id, 0);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Vertical SIMD Automaton — AVX2 256-bit mismatch budget tracking
// ═══════════════════════════════════════════════════════════════════════════════

/// Vertical SIMD Automaton: processes up to 256 work items per cycle using AVX2.
/// Collapses the O(W) scalar per-item branch into O(1) SIMD transitions per
/// child base.
///
/// **Match Mask Precomputation**: Four 256-bit registers (M_A, M_C, M_G, M_T)
/// encode which query's character at the current depth is each nucleotide.
/// Bit i is set in M_x iff `work_items[i].target_char == x`.
///
/// **Budget State**: Four 256-bit registers (R_0..R_3) encode mismatch budget.
/// Bit i is set in R_j iff `work_items[i].budget >= j`.
///
/// **Recurrence** (per child base c, budget-remaining semantics):
///   R'_3 = R_3 & M_c                        (budget=3 survives only on match)
///   R'_2 = (R_2 & M_c) | (R_3 & ~M_c)      (match keeps 2, mismatch demotes 3→2)
///   R'_1 = (R_1 & M_c) | (R_2 & ~M_c)      (match keeps 1, mismatch demotes 2→1)
///   R'_0 = (R_0 & M_c) | (R_1 & ~M_c)      (match keeps 0, mismatch demotes 1→0)
///
/// **13-register map** (of 16 YMM):
///   ymm0-3:  M_A, M_C, M_G, M_T  (match masks, read-only after setup)
///   ymm4-7:  R_0, R_1, R_2, R_3  (budget state, read-only per group)
///   ymm8:    all-ones constant    (for NOT via XOR)
///   ymm9:    M_c                  (selected match mask for current base)
///   ymm10:   ~M_c                 (complement)
///   ymm11:   scratch              (recurrence temporaries)
///   ymm12:   scratch / any        (union for testz pruning)
///
/// **Branchless pruning**: `_mm256_testz_si256(any, any)` prunes dead children
/// in one instruction — no branch misprediction.
///
/// **Survivor extraction**: store 256-bit register to `[u64; 4]`, iterate only
/// set bits via `trailing_zeros()` + `word &= word - 1`. Loop count equals
/// survivor count, not 256.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn vertical_automaton_avx2(
    work_items: &[WorkItem],
    alive_mask: i32,
    base_survivors: &mut [Vec<(u32, u8, u8)>; 4],
) {
    use std::arch::x86_64::*;

    let w = work_items.len();
    if w == 0 { return; }

    let mut chunk_base = 0;
    while chunk_base < w {
        let chunk_end = (chunk_base + 256).min(w);
        let chunk_len = chunk_end - chunk_base;
        let chunk = &work_items[chunk_base..chunk_end];

        // ── Match Mask Precomputation ────────────────────────────────
        let mut ma = [0u64; 4];
        let mut mc = [0u64; 4];
        let mut mg = [0u64; 4];
        let mut mt = [0u64; 4];
        let mut b0 = [0u64; 4];
        let mut b1 = [0u64; 4];
        let mut b2 = [0u64; 4];
        let mut b3 = [0u64; 4];

        for (i, item) in chunk.iter().enumerate() {
            let word = i >> 6;
            let bit = 1u64 << (i & 63);
            match item.target_char {
                b'A' => ma[word] |= bit,
                b'C' => mc[word] |= bit,
                b'G' => mg[word] |= bit,
                _    => mt[word] |= bit,
            }
            // R_j: bit set iff budget >= j. max_mm ≤ 3, so 4 levels.
            // When max_mm == 0, only b0 is populated; b1..b3 stay zero.
            b0[word] |= bit;
            if item.budget >= 1 { b1[word] |= bit; }
            if item.budget >= 2 { b2[word] |= bit; }
            if item.budget >= 3 { b3[word] |= bit; }
        }

        // ymm0-3: match masks
        let m_masks: [__m256i; 4] = [
            _mm256_loadu_si256(ma.as_ptr() as *const __m256i),
            _mm256_loadu_si256(mc.as_ptr() as *const __m256i),
            _mm256_loadu_si256(mg.as_ptr() as *const __m256i),
            _mm256_loadu_si256(mt.as_ptr() as *const __m256i),
        ];
        // ymm4-7: budget state
        let r0 = _mm256_loadu_si256(b0.as_ptr() as *const __m256i);
        let r1 = _mm256_loadu_si256(b1.as_ptr() as *const __m256i);
        let r2 = _mm256_loadu_si256(b2.as_ptr() as *const __m256i);
        let r3 = _mm256_loadu_si256(b3.as_ptr() as *const __m256i);
        // ymm8: all-ones for NOT
        let all_ones = _mm256_set1_epi64x(-1i64);

        // ── Per-base recurrence + branchless pruning ─────────────────
        for bi in 0..4usize {
            if alive_mask & (1 << bi) == 0 { continue; }

            let m_c = m_masks[bi];                              // ymm9
            let not_mc = _mm256_xor_si256(m_c, all_ones);      // ymm10

            // Budget-remaining recurrence: mismatch demotes from j+1.
            // R'_3 = R_3 & M_c
            let r3_new = _mm256_and_si256(r3, m_c);
            // R'_2 = (R_2 & M_c) | (R_3 & ~M_c)
            let r2_new = _mm256_or_si256(
                _mm256_and_si256(r2, m_c),
                _mm256_and_si256(r3, not_mc),
            );
            // R'_1 = (R_1 & M_c) | (R_2 & ~M_c)
            let r1_new = _mm256_or_si256(
                _mm256_and_si256(r1, m_c),
                _mm256_and_si256(r2, not_mc),
            );
            // R'_0 = (R_0 & M_c) | (R_1 & ~M_c)
            let r0_new = _mm256_or_si256(
                _mm256_and_si256(r0, m_c),
                _mm256_and_si256(r1, not_mc),
            );

            // Branchless prune: any survivors?
            let any = _mm256_or_si256(
                _mm256_or_si256(r0_new, r1_new),
                _mm256_or_si256(r2_new, r3_new),
            );
            if _mm256_testz_si256(any, any) != 0 { continue; }

            // ── Survivor extraction via trailing_zeros bit scan ──
            // Extract qwords directly from YMM registers to avoid
            // Store-to-Load Forwarding stalls. The _mm256_storeu_si256
            // → u64 load pattern forces a pipeline drain (~10 cycles)
            // because the store buffer hasn't committed the 32-byte
            // YMM write before the 8-byte scalar read at the same
            // address. _mm256_extract_epi64 reads the register file
            // directly — zero memory round-trip.
            let surv = [
                _mm256_extract_epi64::<0>(any) as u64,
                _mm256_extract_epi64::<1>(any) as u64,
                _mm256_extract_epi64::<2>(any) as u64,
                _mm256_extract_epi64::<3>(any) as u64,
            ];
            let not_mc_bits = [
                _mm256_extract_epi64::<0>(not_mc) as u64,
                _mm256_extract_epi64::<1>(not_mc) as u64,
                _mm256_extract_epi64::<2>(not_mc) as u64,
                _mm256_extract_epi64::<3>(not_mc) as u64,
            ];

            for qword in 0..4usize {
                let mut bits = surv[qword];
                while bits != 0 {
                    let bit_pos = bits.trailing_zeros() as usize;
                    let idx = qword * 64 + bit_pos;
                    if idx >= chunk_len { break; }
                    let item = &chunk[idx];
                    let cost = ((not_mc_bits[qword] >> bit_pos) & 1) as u8;
                    base_survivors[bi].push((item.vqid, item.budget - cost, cost));
                    bits &= bits - 1;
                }
            }
        }

        chunk_base = chunk_end;
    }
}

/// One depth step with lineage-based deduplication — SIMD vectorized.
///
/// On x86_64 with AVX2: dispatches to `vertical_automaton_avx2` which processes
/// up to 256 work items per cycle using the vertical SIMD mismatch recurrence.
/// Falls back to per-item SSE2 survivor computation without AVX2.
///
/// LF mapping uses SSE2 (128-bit) or BlockRank for cache-line-aligned rank queries.
/// Double-lookahead prefetch hides DRAM latency behind compute.
///
/// Falls back to scalar path on non-x86_64 or if `rank_data()` returns `None`.
pub fn step_depth_dedup<I: FmOcc>(
    index: &I,
    frontier: &SearchFrontier,
    depth: usize,
    query_len: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
    next: &mut SearchFrontier,
    lineage: &mut LineageTable,
    clog_vqids: &mut Vec<u32>,
    ws: &mut StepWorkspace,
) {
    let n = frontier.len();
    if n == 0 {
        return;
    }
    next.clear();

    let col = query_len - 1 - depth;
    let is_terminal = depth + 1 == query_len;
    let seqs: [&[Vec<u8>]; 2] = [queries_fwd, queries_rc];

    // ── Sort frontier by (l, r) — zero-alloc via workspace reuse ────────
    let mut order = std::mem::take(&mut ws.order);
    order.clear();
    order.extend(0..n as u32);
    order.sort_unstable_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        frontier.l[a].cmp(&frontier.l[b])
            .then(frontier.r[a].cmp(&frontier.r[b]))
    });

    // ── Precompute group boundaries — zero-alloc via workspace reuse ────
    let mut group_starts = std::mem::take(&mut ws.group_starts);
    group_starts.clear();
    group_starts.push(0);
    for i in 1..n {
        let prev = order[i - 1] as usize;
        let curr = order[i] as usize;
        if frontier.l[curr] != frontier.l[prev] || frontier.r[curr] != frontier.r[prev] {
            group_starts.push(i);
        }
    }
    let n_groups = group_starts.len();

    // Reusable buffers — allocated once, cleared per group.
    let mut work_items: Vec<WorkItem> = Vec::new();
    let mut bkt: [Vec<u32>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    let mut base_survivors: [Vec<(u32, u8, u8)>; 4] =
        [Vec::new(), Vec::new(), Vec::new(), Vec::new()];

    // ══════════════════════════════════════════════════════════════════════
    //  SIMD fast path (x86_64 only)
    //  Supports two rank layouts via RankMode:
    //    - Blocks: 64-byte cache-line-aligned RankBlocks (12 MB, L3-resident)
    //    - Flat: legacy interleaved u32 array (192 MB, DRAM-bound)
    // ══════════════════════════════════════════════════════════════════════
    #[cfg(target_arch = "x86_64")]
    {
        enum RankMode<'a> {
            Blocks(&'a [RankBlock]),
            Flat(*const u32),
        }

        let mode_and_less: Option<(RankMode<'_>, &[usize; 256])> =
            if let Some((blocks, less_tbl)) = index.rank_blocks() {
                Some((RankMode::Blocks(blocks), less_tbl))
            } else if let Some((data, less_tbl)) = index.rank_data() {
                Some((RankMode::Flat(data.as_ptr()), less_tbl))
            } else {
                None
            };

        if let Some((ref rank_mode, less_tbl)) = mode_and_less {
            unsafe {
                use std::arch::x86_64::*;

                // SIMD constants — computed once, reused every group.
                let less_v = _mm_set_epi32(
                    less_tbl[b'T' as usize] as i32,
                    less_tbl[b'G' as usize] as i32,
                    less_tbl[b'C' as usize] as i32,
                    less_tbl[b'A' as usize] as i32,
                );
                let bases_idx_v = _mm_set_epi32(3, 2, 1, 0);
                let ones_v = _mm_set1_epi32(1);
                let has_avx2 = is_x86_feature_detected!("avx2");

                // ── Software-Pipelined Prefetch (depth = 2) ─────────
                //
                // Pipeline structure:
                //   Prologue: issue T0 prefetches for groups 0 and 1.
                //   Steady state: at group g, issue T0 for group g+2.
                //   By the time we compute group g, its RankBlocks were
                //   prefetched into L1 two iterations ago (or in the
                //   prologue). ~100–300 cycles of compute per group
                //   exceeds the ~40–60 cycle L3→L1 promotion latency.
                //
                // Address speculation (BlockRank layout):
                //   block_idx  = pos / BLOCK_SIZE = pos >> 6
                //   byte_addr  = block_idx × 64   = (pos >> 6) << 6
                //   Each RankBlock is repr(C, align(64)) — prefetch
                //   targets are always cache-line aligned (no splits).
                //
                // Why T0 everywhere (not T1 for g+2):
                //   With depth 2, the data has 2 full iterations to
                //   arrive. T0 places it directly in L1 (0-cycle access
                //   on hit) vs T1 which stops at L2 (~4-cycle penalty).
                //   L1 pressure is negligible: 2–4 cache lines (128–256
                //   bytes) vs 32 KB L1D capacity.

                // Prologue: prime the first 2 groups into L1.
                for p in 0..n_groups.min(2) {
                    let pgi = group_starts[p];
                    let pl = frontier.l[order[pgi] as usize] as usize;
                    let pr = frontier.r[order[pgi] as usize] as usize;
                    match rank_mode {
                        RankMode::Blocks(blocks) => {
                            if pl > 0 {
                                _mm_prefetch(
                                    &blocks[(pl - 1) / BLOCK_SIZE] as *const _ as *const i8,
                                    _MM_HINT_T0,
                                );
                            }
                            _mm_prefetch(
                                &blocks[pr / BLOCK_SIZE] as *const _ as *const i8,
                                _MM_HINT_T0,
                            );
                        }
                        RankMode::Flat(rank_ptr) => {
                            if pl > 0 {
                                _mm_prefetch(rank_ptr.add((pl - 1) * 4) as *const i8, _MM_HINT_T0);
                            }
                            _mm_prefetch(rank_ptr.add(pr * 4) as *const i8, _MM_HINT_T0);
                        }
                    }
                }

                for g in 0..n_groups {
                    let gi = group_starts[g];
                    let gj = if g + 1 < n_groups { group_starts[g + 1] } else { n };
                    let l = frontier.l[order[gi] as usize] as usize;
                    let r = frontier.r[order[gi] as usize] as usize;

                    // ── Steady state: prefetch g+2 into L1 ──────────
                    if g + 2 < n_groups {
                        let fgi = group_starts[g + 2];
                        let fl = frontier.l[order[fgi] as usize] as usize;
                        let fr = frontier.r[order[fgi] as usize] as usize;
                        match rank_mode {
                            RankMode::Blocks(blocks) => {
                                // block_addr = (pos >> 6) << 6, 64-byte aligned
                                if fl > 0 {
                                    _mm_prefetch(
                                        &blocks[(fl - 1) / BLOCK_SIZE] as *const _ as *const i8,
                                        _MM_HINT_T0,
                                    );
                                }
                                _mm_prefetch(
                                    &blocks[fr / BLOCK_SIZE] as *const _ as *const i8,
                                    _MM_HINT_T0,
                                );
                            }
                            RankMode::Flat(rank_ptr) => {
                                if fl > 0 {
                                    _mm_prefetch(rank_ptr.add((fl - 1) * 4) as *const i8, _MM_HINT_T0);
                                }
                                _mm_prefetch(rank_ptr.add(fr * 4) as *const i8, _MM_HINT_T0);
                            }
                        }
                    }

                    // ── Build work_items ─────────────────────────────────
                    work_items.clear();
                    for gk in gi..gj {
                        let k = order[gk] as usize;
                        let budget = frontier.budget[k];
                        let lin_id = frontier.query_id[k];
                        for &vqid in lineage.get(lin_id) {
                            let target_char =
                                seqs[(vqid & 1) as usize][(vqid >> 1) as usize][col];
                            work_items.push(WorkItem { vqid, budget, target_char });
                        }
                    }

                    // ── SIMD LF mapping ──────────────────────────────────
                    let (occ_l_v, occ_r_v) = match rank_mode {
                        RankMode::Blocks(blocks) => {
                            let occ_l = if l > 0 {
                                rank_all_blocks_sse(blocks, l - 1)
                            } else {
                                _mm_setzero_si128()
                            };
                            let occ_r = rank_all_blocks_sse(blocks, r);
                            (occ_l, occ_r)
                        }
                        RankMode::Flat(rank_ptr) => {
                            let occ_l = if l > 0 {
                                _mm_loadu_si128(rank_ptr.add((l - 1) * 4) as *const __m128i)
                            } else {
                                _mm_setzero_si128()
                            };
                            let occ_r = _mm_loadu_si128(rank_ptr.add(r * 4) as *const __m128i);
                            (occ_l, occ_r)
                        }
                    };

                    let nl_v = _mm_add_epi32(less_v, occ_l_v);
                    let nr_v = _mm_add_epi32(less_v, occ_r_v);

                    // alive: nl < nr_exc (non-empty child intervals)
                    let alive_v = _mm_cmpgt_epi32(nr_v, nl_v);
                    let alive_mask = _mm_movemask_ps(_mm_castsi128_ps(alive_v));
                    if alive_mask == 0 { continue; }

                    // Store nl, nr for dispatch.
                    let mut nl_arr = [0u32; 4];
                    let mut nr_arr = [0u32; 4];
                    _mm_storeu_si128(nl_arr.as_mut_ptr() as *mut __m128i, nl_v);
                    _mm_storeu_si128(nr_arr.as_mut_ptr() as *mut __m128i, nr_v);

                    for buf in base_survivors.iter_mut() { buf.clear(); }

                    // ── Vertical SIMD Automaton (AVX2) / SSE2 fallback ───
                    if has_avx2 {
                        vertical_automaton_avx2(
                            &work_items, alive_mask,
                            &mut base_survivors,
                        );
                    } else {
                        for item in &work_items {
                            let tidx = base_char_to_idx(item.target_char) as i32;
                            let target_v = _mm_set1_epi32(tidx);
                            let match_v = _mm_cmpeq_epi32(target_v, bases_idx_v);
                            let cost_v = _mm_andnot_si128(match_v, ones_v);

                            let budget_v = _mm_set1_epi32(item.budget as i32);
                            let remaining_v = _mm_sub_epi32(budget_v, cost_v);

                            let over_budget =
                                _mm_cmpgt_epi32(_mm_setzero_si128(), remaining_v);
                            let viable_v = _mm_andnot_si128(over_budget, alive_v);

                            let mask = _mm_movemask_ps(_mm_castsi128_ps(viable_v)) as u32;
                            if mask == 0 { continue; }

                            let mut rem = [0i32; 4];
                            _mm_storeu_si128(rem.as_mut_ptr() as *mut __m128i, remaining_v);
                            let mut costs = [0i32; 4];
                            _mm_storeu_si128(costs.as_mut_ptr() as *mut __m128i, cost_v);

                            let mut m = mask;
                            while m != 0 {
                                let bi = m.trailing_zeros() as usize;
                                base_survivors[bi].push((item.vqid, rem[bi] as u8, costs[bi] as u8));
                                m &= m - 1;
                            }
                        }
                    }

                    // ── Dispatch ─────────────────────────────────────────
                    dispatch_survivors(
                        &base_survivors, &nl_arr, &nr_arr, is_terminal,
                        query_len, max_mm, index, chroms, hits, next, lineage,
                        &mut bkt, clog_vqids,
                    );
                }
            }

            if !is_terminal && next.len() > 1 {
                dedup_frontier(next, lineage);
            }
            // Return workspace vectors (capacity preserved for next depth step).
            ws.order = order;
            ws.group_starts = group_starts;
            return;
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    //  Scalar fallback (non-x86_64 or no InterleavedRank)
    // ══════════════════════════════════════════════════════════════════════
    for g in 0..n_groups {
        let gi = group_starts[g];
        let gj = if g + 1 < n_groups { group_starts[g + 1] } else { n };
        let l = frontier.l[order[gi] as usize] as usize;
        let r = frontier.r[order[gi] as usize] as usize;

        // Prefetch next group.
        if g + 1 < n_groups {
            let ngi = group_starts[g + 1];
            let next_l = frontier.l[order[ngi] as usize] as usize;
            let next_r = frontier.r[order[ngi] as usize] as usize;
            index.prefetch_lf(next_l, next_r);
        }

        // Build work_items.
        work_items.clear();
        for gk in gi..gj {
            let k = order[gk] as usize;
            let budget = frontier.budget[k];
            let lin_id = frontier.query_id[k];
            for &vqid in lineage.get(lin_id) {
                let target_char = seqs[(vqid & 1) as usize][(vqid >> 1) as usize][col];
                work_items.push(WorkItem { vqid, budget, target_char });
            }
        }

        // Scalar LF mapping + per-base survivor scan.
        let lf = index.lf_map(l, r);
        for buf in base_survivors.iter_mut() { buf.clear(); }
        for (bi, &base) in BASES.iter().enumerate() {
            let (nl, nr_exc) = lf[bi];
            if nl >= nr_exc { continue; }
            for item in &work_items {
                let mm_cost = (base != item.target_char) as u8;
                if item.budget >= mm_cost {
                    base_survivors[bi].push((item.vqid, item.budget - mm_cost, mm_cost));
                }
            }
        }
        let nl_arr = [lf[0].0, lf[1].0, lf[2].0, lf[3].0];
        let nr_arr = [lf[0].1, lf[1].1, lf[2].1, lf[3].1];
        dispatch_survivors(
            &base_survivors, &nl_arr, &nr_arr, is_terminal,
            query_len, max_mm, index, chroms, hits, next, lineage, &mut bkt,
            clog_vqids,
        );
    }

    if !is_terminal && next.len() > 1 {
        dedup_frontier(next, lineage);
    }

    // Return buffers to workspace (capacity preserved for next depth step).
    ws.order = order;
    ws.group_starts = group_starts;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Wavefront Processor — block-address sorted linear streaming + 8-wide AVX2 rank
// ═══════════════════════════════════════════════════════════════════════════════

/// Drop-in replacement for `step_depth_dedup` using the wavefront strategy:
///
/// 1. **Sort by block_idx(l-1)** — memory access sweeps linearly through the
///    RankBlock array so the hardware prefetcher pulls data into L1 automatically.
/// 2. **Detect groups of identical (l,r)** — work-item sharing is preserved.
/// 3. **Batch 8 groups' rank computations** via `rank_8way_avx2` — two gather
///    batches (l-1 positions, then r positions) amortise latency across 8 groups.
/// 4. Per-group: existing `vertical_automaton_avx2` + `dispatch_survivors`.
/// 5. `dedup_frontier` on next (unchanged).
///
/// Requires AVX2 and `rank_blocks() → Some`. Falls back to `step_depth_dedup`
/// at the call site when unavailable.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn step_depth_wavefront<I: FmOcc>(
    index: &I,
    frontier: &SearchFrontier,
    depth: usize,
    query_len: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
    next: &mut SearchFrontier,
    lineage: &mut LineageTable,
    clog_vqids: &mut Vec<u32>,
    ws: &mut StepWorkspace,
) {
    use std::arch::x86_64::*;

    let n = frontier.len();
    if n == 0 {
        return;
    }
    next.clear();

    let (blocks, less_tbl) = index.rank_blocks().unwrap();
    let blocks_ptr = blocks.as_ptr() as *const i32;

    let col = query_len - 1 - depth;
    let is_terminal = depth + 1 == query_len;
    let seqs: [&[Vec<u8>]; 2] = [queries_fwd, queries_rc];

    // Precompute less values for the 4 bases (A, C, G, T).
    let less_arr: [u32; 4] = [
        less_tbl[b'A' as usize] as u32,
        less_tbl[b'C' as usize] as u32,
        less_tbl[b'G' as usize] as u32,
        less_tbl[b'T' as usize] as u32,
    ];

    // ── Sort frontier by block_idx(l-1) for linear streaming ────────────
    let mut order = std::mem::take(&mut ws.order);
    order.clear();
    order.extend(0..n as u32);
    order.sort_unstable_by(|&a, &b| {
        let (a, b) = (a as usize, b as usize);
        let blk_a = if frontier.l[a] > 0 { (frontier.l[a] - 1) >> 6 } else { u32::MAX };
        let blk_b = if frontier.l[b] > 0 { (frontier.l[b] - 1) >> 6 } else { u32::MAX };
        blk_a.cmp(&blk_b)
            .then(frontier.l[a].cmp(&frontier.l[b]))
            .then(frontier.r[a].cmp(&frontier.r[b]))
    });

    // ── Detect groups of identical (l, r) ───────────────────────────────
    let mut group_starts = std::mem::take(&mut ws.group_starts);
    group_starts.clear();
    group_starts.push(0);
    for i in 1..n {
        let prev = order[i - 1] as usize;
        let curr = order[i] as usize;
        if frontier.l[curr] != frontier.l[prev] || frontier.r[curr] != frontier.r[prev] {
            group_starts.push(i);
        }
    }
    let n_groups = group_starts.len();

    // Reusable buffers.
    let mut work_items: Vec<WorkItem> = Vec::new();
    let mut bkt: [Vec<u32>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    let mut base_survivors: [Vec<(u32, u8, u8)>; 4] =
        [Vec::new(), Vec::new(), Vec::new(), Vec::new()];

    // ── Process groups in batches of 8 ──────────────────────────────────
    let mut g = 0usize;
    while g < n_groups {
        let batch_end = (g + 8).min(n_groups);
        let batch_len = batch_end - g;

        // Collect l-1 and r positions for rank_8way.
        let mut pos_l = [u32::MAX; 8]; // u32::MAX = sentinel for l=0
        let mut pos_r = [0u32; 8];
        for k in 0..batch_len {
            let gi = group_starts[g + k];
            let idx = order[gi] as usize;
            let l = frontier.l[idx];
            let r = frontier.r[idx];
            pos_l[k] = if l > 0 { l - 1 } else { u32::MAX };
            pos_r[k] = r;
        }
        // Pad unused lanes with safe values (will be ignored).
        for k in batch_len..8 {
            pos_l[k] = u32::MAX;
            pos_r[k] = 0;
        }

        // 16-wide interleaved rank computation — overlaps l/r gather latency.
        let (ranks_l, ranks_r) = rank_16way_interleaved(blocks_ptr, pos_l, pos_r);

        // Process each group in this batch.
        for k in 0..batch_len {
            let group_idx = g + k;
            let gi = group_starts[group_idx];
            let gj = if group_idx + 1 < n_groups { group_starts[group_idx + 1] } else { n };
            let first = order[gi] as usize;
            let _l = frontier.l[first] as usize;
            let _r = frontier.r[first] as usize;

            // Compute nl[4] and nr[4] from precomputed ranks + less table.
            let mut nl_arr = [0u32; 4];
            let mut nr_arr = [0u32; 4];
            for bi in 0..4 {
                nl_arr[bi] = less_arr[bi] + ranks_l[k][bi];
                nr_arr[bi] = less_arr[bi] + ranks_r[k][bi];
            }

            // alive: nl < nr_exc (non-empty child intervals)
            let nl_v = _mm_loadu_si128(nl_arr.as_ptr() as *const __m128i);
            let nr_v = _mm_loadu_si128(nr_arr.as_ptr() as *const __m128i);
            let alive_v = _mm_cmpgt_epi32(nr_v, nl_v);
            let alive_mask = _mm_movemask_ps(_mm_castsi128_ps(alive_v));
            if alive_mask == 0 { continue; }

            // Build work_items from lineage.
            work_items.clear();
            for gk in gi..gj {
                let idx = order[gk] as usize;
                let budget = frontier.budget[idx];
                let lin_id = frontier.query_id[idx];
                for &vqid in lineage.get(lin_id) {
                    let target_char =
                        seqs[(vqid & 1) as usize][(vqid >> 1) as usize][col];
                    work_items.push(WorkItem { vqid, budget, target_char });
                }
            }

            // Vertical SIMD Automaton.
            for buf in base_survivors.iter_mut() { buf.clear(); }
            vertical_automaton_avx2(
                &work_items, alive_mask,
                &mut base_survivors,
            );

            // Dispatch survivors.
            dispatch_survivors(
                &base_survivors, &nl_arr, &nr_arr, is_terminal,
                query_len, max_mm, index, chroms, hits, next, lineage,
                &mut bkt, clog_vqids,
            );
        }

        g = batch_end;
    }

    if !is_terminal && next.len() > 1 {
        dedup_frontier(next, lineage);
    }

    ws.order = order;
    ws.group_starts = group_starts;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Hit Deduplication
// ═══════════════════════════════════════════════════════════════════════════════

/// Remove duplicate (query_id, position) pairs, keeping the best score.
fn dedup_hits(hits: &mut HitAccumulator) {
    let n = hits.query_id.len();
    if n <= 1 { return; }

    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_unstable_by(|&a, &b| {
        hits.query_id[a].cmp(&hits.query_id[b])
            .then(hits.position[a].cmp(&hits.position[b]))
            .then(hits.strand[a].cmp(&hits.strand[b]))
    });

    let mut new_qi = Vec::with_capacity(n);
    let mut new_pos = Vec::with_capacity(n);
    let mut new_strand = Vec::with_capacity(n);
    let mut new_score = Vec::with_capacity(n);

    let mut prev_qi = u32::MAX;
    let mut prev_pos = u32::MAX;
    let mut prev_strand = false;

    for &i in &idx {
        let qi = hits.query_id[i];
        let pos = hits.position[i];
        let strand = hits.strand[i];
        if qi == prev_qi && pos == prev_pos && strand == prev_strand {
            // Duplicate — keep best score (highest = fewest mismatches).
            let last = new_score.len() - 1;
            if hits.score[i] > new_score[last] {
                new_score[last] = hits.score[i];
            }
            continue;
        }
        new_qi.push(qi);
        new_pos.push(pos);
        new_strand.push(strand);
        new_score.push(hits.score[i]);
        prev_qi = qi;
        prev_pos = pos;
        prev_strand = strand;
    }

    hits.query_id = new_qi;
    hits.position = new_pos;
    hits.strand = new_strand;
    hits.score = new_score;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Driver: width-first search with k-mer seeding + mmap frontier
// ═══════════════════════════════════════════════════════════════════════════════

/// Width-first BWT search with tunable mismatch-branch width cap.
///
/// Expands the full query length through the BWT. Mismatch branches whose
/// interval width exceeds `max_width` are pruned — these correspond to
/// repetitive genomic regions that are biologically uninformative.
/// Exact-match branches are never pruned regardless of width.
pub fn search_width_first_seeded<I: FmOcc>(
    index: &I,
    _seed_table: &KmerSeedTable,
    _pos_table: &PosTable,
    _genome_text: &[u8],
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    max_width: u32,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) -> io::Result<()> {
    search_width_first(
        index, queries_fwd, queries_rc, query_len,
        max_mm, max_width, chroms, hits,
    );
    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Legacy driver (no k-mer seeding) — kept for API compatibility
// ═══════════════════════════════════════════════════════════════════════════════

pub fn search_width_first<I: FmOcc>(
    index: &I,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    max_width: u32,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) {
    debug_assert_eq!(queries_fwd.len(), queries_rc.len());
    debug_assert!(max_mm <= 3);

    let sa_len = index.sa_len() as u32;
    let mut a = SearchFrontier::seed(queries_fwd.len(), sa_len, max_mm);
    let mut b = SearchFrontier::with_capacity(a.len() * 2);
    let mut ws = StepWorkspace::new();

    for depth in 0..query_len {
        step_depth(
            index, &a, depth, query_len, queries_fwd, queries_rc,
            max_mm, max_width, chroms, hits, &mut b, &mut ws,
        );
        std::mem::swap(&mut a, &mut b);
    }
}
