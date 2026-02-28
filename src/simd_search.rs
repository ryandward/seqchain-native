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

use crate::fm_index::{RankBlock, BLOCK_SIZE};
use crate::kmer_index::{self, KmerSeed, KmerSeedTable, PosTable};

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
            let valid = (nl < nr_exc) & (frontier.budget[i] >= mm_cost);
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
//  Verify-from-Leading-Seeds (Pigeonhole Supplement)
// ═══════════════════════════════════════════════════════════════════════════════

/// Verify pass: search the LEADING K-mer of each query with up to `seed_mm`
/// mismatches using the PosTable, then verify the full query against the
/// stored genome text character-by-character.
///
/// This supplements the BWT extension pass (which seeds from the TRAILING K-mer).
/// Together they satisfy the pigeonhole principle: for max_mm=k mismatches across
/// 2 non-overlapping segments, at least one segment has ≤ floor(k/2) mismatches.
///
/// Parallelized with rayon: queries are processed in chunks across threads.
fn verify_leading_seeds(
    genome_text: &[u8],
    pos_table: &PosTable,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) {
    use rayon::prelude::*;

    let k = pos_table.k();
    // Dynamic seed_mm based on K vs query_len overlap:
    //   Non-overlapping (2K ≤ L): seed_mm=1 suffices by pigeonhole
    //   Overlapping (2K > L):     seed_mm=2 needed for coverage
    let seed_mm = if k * 2 <= query_len { max_mm.min(1) } else { max_mm.min(2) };
    let text_len = genome_text.len();
    let n_queries = queries_fwd.len();

    // Process queries in parallel; each thread collects (qid, pos, strand, score).
    let chunk_results: Vec<Vec<(u32, u32, bool, f32)>> = (0..n_queries)
        .into_par_iter()
        .map(|qid| {
            let mut local_hits: Vec<(u32, u32, bool, f32)> = Vec::new();
            let mut variants: Vec<(usize, u8)> = Vec::with_capacity(64);

            // Forward strand: leading K-mer = query_fwd[0..K]
            let fwd_leading = &queries_fwd[qid][0..k];
            variants.clear();
            kmer_index::enumerate_kmer_variants(fwd_leading, seed_mm, &mut variants);
            for &(rank, _mm_used) in &variants {
                for &gpos in pos_table.positions_for_rank(rank) {
                    let start = gpos as usize;
                    if start + query_len > text_len { continue; }
                    if !chroms.is_valid(start, query_len) { continue; }

                    let genome = &genome_text[start..start + query_len];
                    let query = &queries_fwd[qid];
                    let mut mm = 0u8;
                    let mut ok = true;
                    for j in 0..query_len {
                        if genome[j] != query[j] {
                            mm += 1;
                            if mm > max_mm { ok = false; break; }
                        }
                    }
                    if ok {
                        local_hits.push((qid as u32, start as u32, true, SCORE_LUT[mm as usize]));
                    }
                }
            }

            // Reverse complement strand: leading K-mer = query_rc[0..K]
            let rc_leading = &queries_rc[qid][0..k];
            variants.clear();
            kmer_index::enumerate_kmer_variants(rc_leading, seed_mm, &mut variants);
            for &(rank, _mm_used) in &variants {
                for &gpos in pos_table.positions_for_rank(rank) {
                    let start = gpos as usize;
                    if start + query_len > text_len { continue; }
                    if !chroms.is_valid(start, query_len) { continue; }

                    let genome = &genome_text[start..start + query_len];
                    let query = &queries_rc[qid];
                    let mut mm = 0u8;
                    let mut ok = true;
                    for j in 0..query_len {
                        if genome[j] != query[j] {
                            mm += 1;
                            if mm > max_mm { ok = false; break; }
                        }
                    }
                    if ok {
                        local_hits.push((qid as u32, start as u32, false, SCORE_LUT[mm as usize]));
                    }
                }
            }

            local_hits
        })
        .collect();

    // Merge all thread-local results into the accumulator.
    for local in chunk_results {
        for (qi, pos, strand, score) in local {
            hits.push(qi, pos, strand, score);
        }
    }
}

/// Targeted verify for clog-affected queries only: uses seed_mm = max_mm to search
/// the leading K-mer with all possible mismatches. Only called for the small subset
/// of queries that had clog-pruned BWT paths.
fn verify_leading_seeds_targeted(
    genome_text: &[u8],
    pos_table: &PosTable,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    clog_qids: &[u32],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) {
    use rayon::prelude::*;

    let k = pos_table.k();
    let seed_mm = max_mm; // full mismatch budget for clog recovery
    let text_len = genome_text.len();

    let chunk_results: Vec<Vec<(u32, u32, bool, f32)>> = clog_qids
        .par_iter()
        .map(|&qid| {
            let qid = qid as usize;
            let mut local_hits: Vec<(u32, u32, bool, f32)> = Vec::new();
            let mut variants: Vec<(usize, u8)> = Vec::with_capacity(4096);

            // Forward strand
            let fwd_leading = &queries_fwd[qid][0..k];
            variants.clear();
            kmer_index::enumerate_kmer_variants(fwd_leading, seed_mm, &mut variants);
            for &(rank, _) in &variants {
                for &gpos in pos_table.positions_for_rank(rank) {
                    let start = gpos as usize;
                    if start + query_len > text_len { continue; }
                    if !chroms.is_valid(start, query_len) { continue; }
                    let genome = &genome_text[start..start + query_len];
                    let query = &queries_fwd[qid];
                    let mut mm = 0u8;
                    let mut ok = true;
                    for j in 0..query_len {
                        if genome[j] != query[j] {
                            mm += 1;
                            if mm > max_mm { ok = false; break; }
                        }
                    }
                    if ok {
                        local_hits.push((qid as u32, start as u32, true, SCORE_LUT[mm as usize]));
                    }
                }
            }

            // Reverse complement strand
            let rc_leading = &queries_rc[qid][0..k];
            variants.clear();
            kmer_index::enumerate_kmer_variants(rc_leading, seed_mm, &mut variants);
            for &(rank, _) in &variants {
                for &gpos in pos_table.positions_for_rank(rank) {
                    let start = gpos as usize;
                    if start + query_len > text_len { continue; }
                    if !chroms.is_valid(start, query_len) { continue; }
                    let genome = &genome_text[start..start + query_len];
                    let query = &queries_rc[qid];
                    let mut mm = 0u8;
                    let mut ok = true;
                    for j in 0..query_len {
                        if genome[j] != query[j] {
                            mm += 1;
                            if mm > max_mm { ok = false; break; }
                        }
                    }
                    if ok {
                        local_hits.push((qid as u32, start as u32, false, SCORE_LUT[mm as usize]));
                    }
                }
            }

            local_hits
        })
        .collect();

    for local in chunk_results {
        for (qi, pos, strand, score) in local {
            hits.push(qi, pos, strand, score);
        }
    }
}

/// Verify pass for TRAILING K-mer: complements the BWT pass (which seeds trailing
/// K-mer with 0mm only). Enumerates trailing K-mer variants with up to `seed_mm`
/// mismatches, looks up positions in PosTable, verifies full query against genome.
///
/// The alignment start = genome_trailing_pos - (query_len - K).
fn verify_trailing_seeds(
    genome_text: &[u8],
    pos_table: &PosTable,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) {
    let k = pos_table.k();
    // Dynamic seed_mm based on K vs query_len overlap:
    //   Non-overlapping (2K ≤ L): seed_mm=1 suffices by pigeonhole
    //   Overlapping (2K > L):     seed_mm=2 needed for coverage
    let seed_mm = if k * 2 <= query_len { max_mm.min(1) } else { max_mm.min(2) };
    let text_len = genome_text.len();
    let prefix_len = query_len - k; // offset from alignment start to trailing K-mer

    let mut variants: Vec<(usize, u8)> = Vec::with_capacity(64);

    for qid in 0..queries_fwd.len() {
        // Forward strand: trailing K-mer = query_fwd[query_len-K..]
        let fwd_trailing = &queries_fwd[qid][prefix_len..];
        variants.clear();
        kmer_index::enumerate_kmer_variants(fwd_trailing, seed_mm, &mut variants);
        for &(rank, _mm_used) in &variants {
            for &gpos in pos_table.positions_for_rank(rank) {
                let kmer_start = gpos as usize;
                // Alignment starts at kmer_start - prefix_len
                if kmer_start < prefix_len { continue; }
                let start = kmer_start - prefix_len;
                if start + query_len > text_len { continue; }
                if !chroms.is_valid(start, query_len) { continue; }

                let genome = &genome_text[start..start + query_len];
                let query = &queries_fwd[qid];
                let mut mm = 0u8;
                let mut ok = true;
                for j in 0..query_len {
                    if genome[j] != query[j] {
                        mm += 1;
                        if mm > max_mm { ok = false; break; }
                    }
                }
                if ok {
                    hits.push(qid as u32, start as u32, true, SCORE_LUT[mm as usize]);
                }
            }
        }

        // Reverse complement strand: trailing K-mer = query_rc[query_len-K..]
        let rc_trailing = &queries_rc[qid][prefix_len..];
        variants.clear();
        kmer_index::enumerate_kmer_variants(rc_trailing, seed_mm, &mut variants);
        for &(rank, _mm_used) in &variants {
            for &gpos in pos_table.positions_for_rank(rank) {
                let kmer_start = gpos as usize;
                if kmer_start < prefix_len { continue; }
                let start = kmer_start - prefix_len;
                if start + query_len > text_len { continue; }
                if !chroms.is_valid(start, query_len) { continue; }

                let genome = &genome_text[start..start + query_len];
                let query = &queries_rc[qid];
                let mut mm = 0u8;
                let mut ok = true;
                for j in 0..query_len {
                    if genome[j] != query[j] {
                        mm += 1;
                        if mm > max_mm { ok = false; break; }
                    }
                }
                if ok {
                    hits.push(qid as u32, start as u32, false, SCORE_LUT[mm as usize]);
                }
            }
        }
    }
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

/// Threshold: if initial frontier exceeds this, use mmap-backed frontiers.
const MMAP_THRESHOLD: usize = 500_000;

/// Width-first search with on-disk k-mer seeding.
///
/// 1. Load (or build) the 10-mer seed table.
/// 2. For each query × strand, generate mismatch neighborhood seeds at depth K.
/// 3. If frontier is small, run in-memory. If large, use mmap ping/pong.
/// 4. Step depths K..query_len with branchless 4-way expansion.
pub fn search_width_first_seeded<I: FmOcc>(
    index: &I,
    seed_table: &KmerSeedTable,
    pos_table: &PosTable,
    genome_text: &[u8],
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) -> io::Result<()> {
    debug_assert_eq!(queries_fwd.len(), queries_rc.len());
    debug_assert!(max_mm <= 3);

    let k = seed_table.k();

    if query_len <= k {
        // Query is shorter than seed length — fall back to direct table lookup.
        return search_short_queries(index, seed_table, queries_fwd, queries_rc,
                                     query_len, max_mm, chroms, hits);
    }

    // ── Phase 0: BWT extension from trailing seed ──────────────────────
    // Seeds the trailing K-mer with up to `seed_mm` mismatches, then gives
    // the remaining budget to width-first BWT extension through the leading
    // positions. This covers all alignments where mm_trailing ≤ seed_mm.
    //
    // The verify pass (Phase 1) covers mm_leading ≤ seed_mm via PosTable,
    // which together with this phase satisfies the pigeonhole principle:
    // for max_mm mismatches across 2 segments, at least one has ≤ floor(max_mm/2).

    // Dynamic seed_mm based on K vs query_len overlap:
    //   Non-overlapping (2K ≤ L): seed_mm=1 suffices by pigeonhole
    //   Overlapping (2K > L):     seed_mm=2 needed for coverage
    let seed_mm = if k * 2 <= query_len { max_mm.min(1) } else { max_mm.min(2) };
    let n_queries = queries_fwd.len();
    let mut seed_l = Vec::new();
    let mut seed_r = Vec::new();
    let mut seed_bud = Vec::new();
    let mut seed_qid = Vec::new();
    let mut seed_str = Vec::new();

    let mut seed_buf: Vec<KmerSeed> = Vec::with_capacity(64);

    for qid in 0..n_queries {
        // Forward strand: trailing K-mer = query[query_len-K..]
        let fwd_kmer = &queries_fwd[qid][query_len - k..];
        seed_buf.clear();
        kmer_index::seed_with_mismatches(seed_table, fwd_kmer, seed_mm, &mut seed_buf);
        for s in &seed_buf {
            seed_l.push(s.l);
            seed_r.push(s.r);
            seed_bud.push(max_mm - s.mm_used); // remaining budget for extension
            seed_qid.push(qid as u32);
            seed_str.push(1u8);
        }

        // Reverse complement strand: trailing K-mer = rc_query[query_len-K..]
        let rc_kmer = &queries_rc[qid][query_len - k..];
        seed_buf.clear();
        kmer_index::seed_with_mismatches(seed_table, rc_kmer, seed_mm, &mut seed_buf);
        for s in &seed_buf {
            seed_l.push(s.l);
            seed_r.push(s.r);
            seed_bud.push(max_mm - s.mm_used); // remaining budget for extension
            seed_qid.push(qid as u32);
            seed_str.push(0u8);
        }
    }

    let initial_size = seed_l.len();
    let remaining_depths = query_len - k; // depths K..query_len-1

    if remaining_depths == 0 {
        // K == query_len: seeds are terminal — resolve SA positions directly.
        for i in 0..initial_size {
            let mm_used = max_mm - seed_bud[i];
            let score = SCORE_LUT[mm_used as usize];
            let fwd = seed_str[i] != 0;
            let lo = seed_l[i] as usize;
            let hi = seed_r[i] as usize + 1;
            for sa_idx in lo..hi {
                let sa_pos = index.sa(sa_idx);
                if chroms.is_valid(sa_pos, query_len) {
                    hits.push(seed_qid[i], sa_pos as u32, fwd, score);
                }
            }
        }
        return Ok(());
    }

    // ── Dispatch: dedup vs simple in-memory ─────────────────────────────
    // Dedup pays for itself when there are enough seeds to share intervals.
    // At 0mm with small batches, the overhead of sort-merge isn't worth it.
    let use_dedup = initial_size > 1000 || max_mm > 0;

    let clog_vqids = if use_dedup {
        search_dedup_from_seeds(
            index, &seed_l, &seed_r, &seed_bud, &seed_qid, &seed_str,
            k, queries_fwd, queries_rc, query_len, max_mm, chroms, hits,
        )
    } else {
        search_inmemory_from_seeds(
            index, &seed_l, &seed_r, &seed_bud, &seed_qid, &seed_str,
            k, queries_fwd, queries_rc, query_len, max_mm, chroms, hits,
        );
        Vec::new()
    };

    // ── Phase 1: Verify-from-leading-seeds (pigeonhole supplement) ─────
    // BWT Phase 0 covers mm_T ≤ seed_mm (trailing K-mer seeded with mismatches).
    // Leading verify covers mm_L ≤ seed_mm (leading K-mer searched via PosTable).
    // By pigeonhole: for any alignment with max_mm mismatches across 2 segments,
    // at least one segment has ≤ floor(max_mm/2) ≤ seed_mm. Combined coverage
    // is complete.
    //
    // Additionally, for queries clog-pruned during BWT (mm_L > seed_mm in clogged
    // intervals), we run a targeted verify with seed_mm = max_mm to recover those hits.
    if max_mm >= 1 {
        verify_leading_seeds(
            genome_text, pos_table,
            queries_fwd, queries_rc,
            query_len, max_mm, chroms, hits,
        );

        // Targeted verify for clog-affected queries: recover hits with mm_L > seed_mm
        // that the BWT pass pruned in wide intervals.
        if !clog_vqids.is_empty() && max_mm >= 2 {
            // Deduplicate clog-affected query IDs.
            let mut clog_qids: Vec<u32> = clog_vqids.iter().map(|v| v >> 1).collect();
            clog_qids.sort_unstable();
            clog_qids.dedup();

            verify_leading_seeds_targeted(
                genome_text, pos_table,
                queries_fwd, queries_rc,
                &clog_qids,
                query_len, max_mm, chroms, hits,
            );
        }

        dedup_hits(hits);
    }

    Ok(())
}

/// In-memory path for small frontiers without dedup (< DEDUP_THRESHOLD rows, 0mm).
fn search_inmemory_from_seeds<I: FmOcc>(
    index: &I,
    sl: &[u32], sr: &[u32], sb: &[u8], sq: &[u32], ss: &[u8],
    start_depth: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) {
    let n = sl.len();
    let mut a = SearchFrontier::with_capacity(n);
    for i in 0..n {
        a.push(sl[i], sr[i], sb[i], sq[i], ss[i]);
    }
    let mut b = SearchFrontier::with_capacity(n * 2);
    let mut ws = StepWorkspace::new();

    for depth in start_depth..query_len {
        step_depth(
            index, &a, depth, query_len, queries_fwd, queries_rc,
            max_mm, chroms, hits, &mut b, &mut ws,
        );
        std::mem::swap(&mut a, &mut b);
    }
}

/// Deduplicated in-memory path: lineage table + sort-merge after each depth.
/// Returns the set of vqids that were clog-pruned during BWT extension.
fn search_dedup_from_seeds<I: FmOcc>(
    index: &I,
    sl: &[u32], sr: &[u32], sb: &[u8], sq: &[u32], ss: &[u8],
    start_depth: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) -> Vec<u32> {
    let n = sl.len();
    let mut lineage = LineageTable::with_capacity(n * 2, n);

    // Seed initial lineages: each seed gets a singleton lineage.
    // vqid = query_id * 2 + strand_bit (0=fwd, 1=rc)
    let mut a = SearchFrontier::with_capacity(n);
    for i in 0..n {
        let vqid = sq[i] * 2 + (1 - ss[i] as u32);
        let lin_id = lineage.alloc(&[vqid]);
        a.push(sl[i], sr[i], sb[i], lin_id, 0);
    }

    // Dedup the initial seed frontier.
    dedup_frontier(&mut a, &mut lineage);

    let mut b = SearchFrontier::with_capacity(a.len() * 2);
    let mut clog_vqids: Vec<u32> = Vec::new();
    let mut ws = StepWorkspace::new();

    for depth in start_depth..query_len {
        step_depth_dedup(
            index, &a, depth, query_len, queries_fwd, queries_rc,
            max_mm, chroms, hits, &mut b, &mut lineage, &mut clog_vqids,
            &mut ws,
        );
        std::mem::swap(&mut a, &mut b);
    }

    clog_vqids
}

/// Mmap path for large frontiers (>= MMAP_THRESHOLD rows).
fn search_mmap_from_seeds<I: FmOcc>(
    index: &I,
    sl: &[u32], sr: &[u32], sb: &[u8], sq: &[u32], ss: &[u8],
    start_depth: usize,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) -> io::Result<()> {
    // Write seeds into the first mmap frontier.
    let mut writer_a = MmapFrontier::new()?;
    for i in 0..sl.len() {
        writer_a.push(sl[i], sr[i], sb[i], sq[i], ss[i])?;
    }

    let mut writer_b = MmapFrontier::new()?;

    for depth in start_depth..query_len {
        // Freeze the current writer into a read-only snapshot.
        let frozen = writer_a.freeze()?;

        // Reset the other writer for this depth's output.
        writer_b.reset()?;

        step_depth_mmap(
            index, &frozen, depth, query_len, queries_fwd, queries_rc,
            max_mm, chroms, hits, &mut writer_b,
        )?;

        // Ping/pong: swap writers. writer_a (now drained) becomes the next output.
        std::mem::swap(&mut writer_a, &mut writer_b);
    }

    Ok(())
}

/// Handle queries shorter than or equal to SEED_K.
fn search_short_queries<I: FmOcc>(
    index: &I,
    _seed_table: &KmerSeedTable,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    max_mm: u8,
    chroms: &ChromGeometry,
    hits: &mut HitAccumulator,
) -> io::Result<()> {
    // For short queries, we can't use the K-mer table directly since it's
    // built for K-length strings. Fall back to the original depth-0 seeding.
    let sa_len = index.sa_len() as u32;
    let mut a = SearchFrontier::seed(queries_fwd.len(), sa_len, max_mm);
    let mut b = SearchFrontier::with_capacity(a.len() * 2);
    let mut ws = StepWorkspace::new();

    for depth in 0..query_len {
        step_depth(
            index, &a, depth, query_len, queries_fwd, queries_rc,
            max_mm, chroms, hits, &mut b, &mut ws,
        );
        std::mem::swap(&mut a, &mut b);
    }
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
            max_mm, chroms, hits, &mut b, &mut ws,
        );
        std::mem::swap(&mut a, &mut b);
    }
}
