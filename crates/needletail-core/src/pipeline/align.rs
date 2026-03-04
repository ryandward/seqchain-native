//! RNA-seq alignment pipeline.
//!
//! `align_fastq` routes the Water from disk to disk:
//!
//! ```text
//! ChunkedFastq  →  BWT search  →  group by query  →
//! extend_batch (8-lane SIMD)  →  RegionSink (SAM / Parquet)
//! ```
//!
//! # O(C) RAM invariant
//!
//! At every point in the loop the live allocation is:
//!
//! | Structure             | Size               |
//! |-----------------------|--------------------|
//! | `chunk`               | C × L bytes        |
//! | `fwd + rc`            | C × 2L bytes       |
//! | `HitAccumulator`      | C × h̄ × 13 bytes  |
//! | Parquet/SAM row buf   | row_group_size × row|
//!
//! where C = `chunk_size`, L = `read_len`, h̄ = average hits per read.
//! All are O(C), independent of N (total reads).
//!
//! # Import hierarchy
//! Orchestration layer: may import all lower layers.

use std::io;
use std::path::Path;

use bio::alphabets::dna;

use crate::engine::affine::{extend_batch, SpliceParams};
use crate::engine::simd_search::{ChromGeometry, FmOcc, HitAccumulator};
use crate::engine::kmer_index::{KmerSeedTable, PosTable};
use crate::io::fastq::{ChunkedFastq, FastqRecord};
use crate::io::sam::score_to_mapq;
use crate::io::RegionSink;
use crate::models::region::{Region, Strand, TagValue};
use crate::{run_search_seeded, SearchError};

// ═══════════════════════════════════════════════════════════════════════════
//  Configuration
// ═══════════════════════════════════════════════════════════════════════════

/// Parameters governing a single alignment run.
pub struct AlignConfig {
    /// Maximum mismatches allowed in the BWT seed search (0–3).
    pub max_mm: u8,
    /// Maximum extension length in bases for affine verification.
    /// Set to 0 to skip affine extension and use BWT scores directly.
    pub max_ext: usize,
    /// Reads per processing cycle. Governs peak RAM: O(chunk_size × read_len).
    pub chunk_size: usize,
    /// Maximum BWT search width (limits multi-mapper explosion for repeat regions).
    pub max_width: u32,
    /// Affine gap / splice parameters.
    pub splice: SpliceParams,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            max_mm:     2,
            max_ext:    150,
            chunk_size: 100_000,
            max_width:  u32::MAX,
            splice:     SpliceParams::default(),
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Statistics
// ═══════════════════════════════════════════════════════════════════════════

/// Summary counters returned by `align_fastq`.
#[derive(Debug, Default)]
pub struct AlignStats {
    /// Total reads processed.
    pub total_reads:  u64,
    /// Reads with at least one alignment.
    pub mapped_reads: u64,
    /// Reads with > 1 alignment (multi-mappers).
    pub multi_mapped: u64,
    /// Reads with no alignment.
    pub unmapped:     u64,
    /// Total alignment records emitted (primary + secondary).
    pub total_hits:   u64,
}

// ═══════════════════════════════════════════════════════════════════════════
//  Main entry point
// ═══════════════════════════════════════════════════════════════════════════

/// Align all reads in `fastq_path` against the FM-index and write results
/// to `sink`.
///
/// # Memory
/// Peak RSS ≈ FM-Index + seed tables + O(chunk_size × read_len × 3).
/// Guides never accumulate across chunks.
///
/// # Arguments
/// * `index`      — Any `FmOcc + Send + Sync` (built or mmap-loaded)
/// * `seed_table` / `pos_table` — K-mer seed acceleration tables
/// * `genome_text` — Raw reference bytes for affine extension
/// * `chroms`     — Chromosome geometry for hit validation
/// * `fastq_path` — Input reads (uncompressed FASTQ)
/// * `sink`       — Terminal output (`SamSink`, `ParquetFileSink`, etc.)
/// * `cfg`        — Tuning parameters
pub fn align_fastq<I>(
    index:       &I,
    seed_table:  &KmerSeedTable,
    pos_table:   &PosTable,
    genome_text: &[u8],
    chroms:      &ChromGeometry,
    fastq_path:  &Path,
    sink:        &mut dyn RegionSink,
    cfg:         &AlignConfig,
) -> Result<AlignStats, AlignError>
where
    I: FmOcc + Send + Sync,
{
    use std::fs::File;
    use std::io::BufReader;

    let reader = BufReader::new(File::open(fastq_path)?);
    let parser = ChunkedFastq::new(reader, cfg.chunk_size);

    let mut stats = AlignStats::default();

    for chunk_result in parser {
        let chunk: Vec<FastqRecord> = chunk_result?;
        // ── chunk lives here: O(C × L) ────────────────────────────────────

        if chunk.is_empty() {
            continue;
        }

        let read_len = chunk[0].seq.len();
        stats.total_reads += chunk.len() as u64;

        // ── Forward + RC sequences: O(C × 2L) ────────────────────────────
        // Allocated once per chunk; freed at the end of this block.
        let (fwd, rc): (Vec<Vec<u8>>, Vec<Vec<u8>>) = chunk
            .iter()
            .map(|r| {
                let f = r.seq.clone();  // already uppercased by ChunkedFastq
                let rc = dna::revcomp(&f);
                (f, rc)
            })
            .unzip();

        // ── BWT search: O(C × h̄ × 13B) ──────────────────────────────────
        // MmapFrontier handles the BWT frontier out-of-core.
        // HitAccumulator is bounded by chunk × avg_hits.
        let hits = run_search_seeded(
            index,
            seed_table,
            pos_table,
            genome_text,
            &fwd,
            &rc,
            read_len,
            cfg.max_mm,
            cfg.max_width,
            chroms,
        )?;

        drop(fwd);  // ← O(C × L) freed
        drop(rc);   // ← O(C × L) freed

        // ── Group hits by query_id: O(H) ──────────────────────────────────
        // Index-sort to avoid allocating a separate Vec<HitEntry> per read.
        let mut idx: Vec<usize> = (0..hits.query_id.len()).collect();
        idx.sort_unstable_by_key(|&i| hits.query_id[i]);

        // ── Emit: extend + write, O(row_group_size) at a time in the sink ─
        process_hits(
            &chunk,
            &hits,
            &idx,
            genome_text,
            read_len,
            cfg,
            sink,
            &mut stats,
        )?;
        // hits, idx, chunk all dropped here.
    }

    Ok(stats)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Inner processing — kept separate to keep `align_fastq` readable
// ═══════════════════════════════════════════════════════════════════════════

/// Extend BWT hits with the affine splice kernel and drain to the sink.
///
/// Reads are batched in groups of 8 for SIMD extension. Within each group
/// the tensor pivot is applied once; the extension loop is then fed by
/// contiguous array loads with no per-step gather.
#[allow(clippy::too_many_arguments)]
fn process_hits(
    chunk:       &[FastqRecord],
    hits:        &HitAccumulator,
    sorted_idx:  &[usize],        // hit indices sorted by query_id
    genome_text: &[u8],
    read_len:    usize,
    cfg:         &AlignConfig,
    sink:        &mut dyn RegionSink,
    stats:       &mut AlignStats,
) -> io::Result<()> {
    let n_reads = chunk.len();

    // ── Build per-read best-hit map ────────────────────────────────────────
    // For each query in [0, n_reads): record its best hit (highest score)
    // and total hit count.
    //
    // best_hit[qid] = Some((hit_index, score))
    // nhits[qid]    = number of BWT hits
    let mut best_hit: Vec<Option<(usize, f32)>> = vec![None; n_reads];
    let mut nhits:    Vec<u32>                   = vec![0;    n_reads];

    // sorted_idx is sorted by query_id → linear group scan.
    let mut si = 0;
    while si < sorted_idx.len() {
        let qid = hits.query_id[sorted_idx[si]] as usize;
        // Gather all indices for this query_id.
        let group_end = sorted_idx[si..]
            .partition_point(|&i| hits.query_id[i] as usize == qid)
            + si;
        let group = &sorted_idx[si..group_end];

        nhits[qid] = group.len() as u32;
        // Best hit = highest BWT score.
        let best_i = *group.iter().max_by(|&&a, &&b| {
            hits.score[a].partial_cmp(&hits.score[b]).unwrap()
        }).unwrap();
        best_hit[qid] = Some((best_i, hits.score[best_i]));

        si = group_end;
    }

    // ── Update stats ──────────────────────────────────────────────────────
    for qid in 0..n_reads {
        if best_hit[qid].is_some() {
            stats.mapped_reads += 1;
            if nhits[qid] > 1 {
                stats.multi_mapped += 1;
            }
            stats.total_hits += nhits[qid] as u64;
        } else {
            stats.unmapped += 1;
        }
    }

    // ── Emit unmapped reads ────────────────────────────────────────────────
    for qid in 0..n_reads {
        if best_hit[qid].is_none() {
            let rec = &chunk[qid];
            let region = unmapped_region(rec);
            sink.consume(region)?;
        }
    }

    // ── Extend mapped reads in SIMD batches of 8 ──────────────────────────
    // Collect all (qid, hit_index) pairs, then iterate in groups of 8.
    let mapped: Vec<usize> = (0..n_reads)
        .filter(|&qid| best_hit[qid].is_some())
        .collect();

    for batch_start in (0..mapped.len()).step_by(8) {
        let batch_end = (batch_start + 8).min(mapped.len());
        let batch_size = batch_end - batch_start;

        // ── Collect reads and anchors for this SIMD batch ─────────────────
        let mut reads_refs: Vec<&[u8]> = Vec::with_capacity(batch_size);
        let mut anchors:    Vec<u32>   = Vec::with_capacity(batch_size);
        let mut bwt_scores: Vec<f32>   = Vec::with_capacity(batch_size);
        let mut strands:    Vec<bool>  = Vec::with_capacity(batch_size); // true = fwd

        for &qid in &mapped[batch_start..batch_end] {
            let (hi, bwt_score) = best_hit[qid].unwrap();
            reads_refs.push(&chunk[qid].seq);
            anchors.push(hits.position[hi]);
            bwt_scores.push(bwt_score);
            strands.push(!hits.strand[hi]); // strand==false in HitAccumulator is forward
        }

        // ── Affine extension (tensor pivot happens inside extend_batch) ────
        let ext_len = cfg.max_ext.min(read_len);
        let extended: Vec<(i32, u32)> = if ext_len > 0 {
            extend_batch(genome_text, &reads_refs, &anchors, &cfg.splice, ext_len)
        } else {
            // Extension disabled: use anchor as endpoint, BWT score as score.
            anchors.iter().map(|&a| (1, a + read_len as u32)).collect()
        };

        // ── Emit primary alignment for each read in the batch ─────────────
        for (lane, &qid) in mapped[batch_start..batch_end].iter().enumerate() {
            let (hi, bwt_score) = best_hit[qid].unwrap();
            let (ext_score, ext_end) = extended[lane];

            // Use extended score if available, else fall back to BWT score.
            let final_score = if ext_len > 0 && ext_score > 0 {
                ext_score as f64 / (cfg.splice.match_score as f64 * read_len as f64)
            } else {
                bwt_score as f64
            };

            let ref_start = anchors[lane] as i64;
            let ref_end   = if ext_len > 0 && ext_score > 0 {
                ext_end as i64 + 1
            } else {
                ref_start + read_len as i64
            };

            let strand = if strands[lane] { Strand::Forward } else { Strand::Reverse };
            let mapq   = score_to_mapq(final_score);
            let rec    = &chunk[qid];
            let mut region = Region::new(
                hits_chrom(hits.position[hi], genome_text.len()),
                ref_start,
                ref_end,
                )
                .with_strand(strand)
                .with_score(final_score)
                .with_name(String::from_utf8_lossy(&rec.id).into_owned());

            region = region
                .with_tag("seq",    TagValue::Str(
                    String::from_utf8_lossy(&rec.seq).into_owned()
                ))
                .with_tag("qual",   TagValue::Str(
                    String::from_utf8_lossy(&rec.qual).into_owned()
                ))
                .with_tag("nm",     TagValue::Int(mismatch_count(bwt_score) as i64))
                .with_tag("mapq",   TagValue::Int(mapq as i64))
                .with_tag("nhits",  TagValue::Int(nhits[qid] as i64));

            // MAPQ = 0 for multi-mappers (SAM convention).
            if nhits[qid] > 1 {
                region = region.with_tag("mapq", TagValue::Int(0));
            }

            sink.consume(region)?;
        }
    }

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════
//  Helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Build an unmapped Region for a read with no BWT hits.
fn unmapped_region(rec: &FastqRecord) -> Region {
    Region::new("*", 0, 0)
        .with_score(0.0)
        .with_name(String::from_utf8_lossy(&rec.id).into_owned())
        .with_tag("seq",  TagValue::Str(String::from_utf8_lossy(&rec.seq).into_owned()))
        .with_tag("qual", TagValue::Str(String::from_utf8_lossy(&rec.qual).into_owned()))
        .with_tag("nhits", TagValue::Int(0))
}

/// Resolve a reference position to a chromosome name.
///
/// In the master-buffer architecture, chromosome boundaries are stored in
/// `ChromGeometry`. Here we use a simplified placeholder that returns the
/// raw offset as a string; callers with a proper `FmIndexSearcher` should
/// resolve via `chrom_name_at(pos)`.
///
/// TODO: accept `&ChromGeometry` + chrom name table for full resolution.
#[inline]
fn hits_chrom(_pos: u32, _genome_len: usize) -> String {
    // Placeholder: real resolution requires the chrom name table passed
    // from the IndexHandle. See `FmIndexSearcher::chrom_at(pos)`.
    "unknown".to_string()
}

/// Kept for forward-compatibility with a future chrom-name lookup.
#[inline]
#[allow(dead_code)]
fn chrom_for_pos(_pos: u32, _chrom_boundaries: &[(&str, usize)]) -> &'static str {
    "unknown"
}

/// Recover approximate mismatch count from BWT SCORE_LUT value.
/// SCORE_LUT[mm] = 1 / (1 + mm) → mm = round(1/score - 1).
#[inline]
fn mismatch_count(score: f32) -> u8 {
    let mm = (1.0 / score - 1.0).round() as u8;
    mm.min(3)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Error type
// ═══════════════════════════════════════════════════════════════════════════

/// Errors that can occur during alignment.
#[derive(Debug)]
pub enum AlignError {
    Io(io::Error),
    Search(SearchError),
}

impl std::fmt::Display for AlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignError::Io(e)     => write!(f, "I/O error: {e}"),
            AlignError::Search(e) => write!(f, "search error: {e}"),
        }
    }
}

impl std::error::Error for AlignError {}

impl From<io::Error> for AlignError {
    fn from(e: io::Error) -> Self { AlignError::Io(e) }
}

impl From<SearchError> for AlignError {
    fn from(e: SearchError) -> Self { AlignError::Search(e) }
}
