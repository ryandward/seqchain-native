//! needletail-core — Pure Rust CRISPR guide design engine.
//!
//! No Python, no FFI. Exposes five foundation pillars plus orchestration
//! helpers for building frontends (PyO3, HTTP, CLI, etc.).
//!
//! ```text
//!   geometry    ← pure 1D coordinate math (leaf, no deps)
//!   chemistry   ← molecular rules: PAM masks, guide IDs (leaf, no deps)
//!   engine/     ← BWT, SIMD search, k-mer seeds (depends on geometry)
//!   io/         ← mmap persistence, Parquet sink (depends on engine)
//!   operations/ ← PAM scanning, off-target validation (depends on all above)
//! ```

#![deny(clippy::all)]

// ─── The five pillars ────────────────────────────────────────────────────────
pub mod annotation;
pub mod chemistry;
pub mod engine;
pub mod error;
pub mod geometry;
pub mod io;
pub mod models;
pub mod operations;
pub mod pipeline;

// ─── Re-exports for convenience ──────────────────────────────────────────────
pub use chemistry::{CompiledPam, PamDirection};
pub use engine::fm_index::FmIndexSearcher;
pub use engine::kmer_index::{sa_sweep_build, KmerSeedTable, PosTable, SEED_K_LARGE, SEED_K_SMALL};
pub use engine::simd_search::{
    search_width_first, search_width_first_seeded, ChromGeometry, FmOcc, HitAccumulator,
};
pub use error::SearchError;
pub use io::persist::MappedIndex;
pub use io::{RegionSink, CountingSink};
pub use io::json::FileSink;
pub use operations::pam_scanner::{
    enrich_hits, filter_hits_by_pam, find_pam_sites, GuideHits,
};

use std::path::Path;
use std::sync::Arc;

use bio::alphabets::dna;

// ─── IndexHandle ────────────────────────────────────────────────────────────

/// Unified handle for both in-memory-built and mmap-loaded FM-Indexes.
pub enum IndexHandle {
    Built(Arc<FmIndexSearcher>),
    Loaded(Arc<MappedIndex>),
}

impl IndexHandle {
    pub fn chrom_geometry(&self) -> ChromGeometry {
        match self {
            IndexHandle::Built(s) => s.chrom_geometry(),
            IndexHandle::Loaded(m) => m.chrom_geometry(),
        }
    }

    pub fn chrom_names(&self) -> Vec<String> {
        match self {
            IndexHandle::Built(s) => s.chrom_names().into_iter().map(|s| s.to_string()).collect(),
            IndexHandle::Loaded(m) => m.chrom_names(),
        }
    }

    pub fn text(&self) -> &[u8] {
        match self {
            IndexHandle::Built(s) => s.text(),
            IndexHandle::Loaded(m) => m.text(),
        }
    }
}

// ─── Seed tier: paired KmerSeedTable + PosTable for a given K ────────────────

pub struct SeedTier {
    pub seed_table: Arc<KmerSeedTable>,
    pub pos_table: Arc<PosTable>,
}

// ─── Generic search helpers ─────────────────────────────────────────────────

pub fn run_search_unseeded<I: FmOcc + Send + Sync>(
    index: &I,
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    mismatches: u8,
    max_width: u32,
    chroms: &ChromGeometry,
) -> HitAccumulator {
    let mut acc = HitAccumulator::new();
    search_width_first(
        index,
        queries_fwd,
        queries_rc,
        query_len,
        mismatches,
        max_width,
        chroms,
        &mut acc,
    );
    acc
}

pub fn run_search_seeded<I: FmOcc + Send + Sync>(
    index: &I,
    seed_table: &KmerSeedTable,
    pos_table: &PosTable,
    genome_text: &[u8],
    queries_fwd: &[Vec<u8>],
    queries_rc: &[Vec<u8>],
    query_len: usize,
    mismatches: u8,
    max_width: u32,
    chroms: &ChromGeometry,
) -> Result<HitAccumulator, SearchError> {
    use rayon::prelude::*;

    let n = queries_fwd.len();
    if n == 0 {
        return Ok(HitAccumulator::new());
    }

    let n_threads = rayon::current_num_threads();
    let chunk_size = ((n + n_threads - 1) / n_threads).max(256);

    let chunks: Vec<(usize, usize)> = (0..n)
        .step_by(chunk_size)
        .map(|s| (s, (s + chunk_size).min(n)))
        .collect();

    let results: Result<Vec<HitAccumulator>, SearchError> = chunks
        .par_iter()
        .map(|&(start, end)| {
            let mut acc = HitAccumulator::new();
            search_width_first_seeded(
                index,
                seed_table,
                pos_table,
                genome_text,
                &queries_fwd[start..end],
                &queries_rc[start..end],
                query_len,
                mismatches,
                max_width,
                chroms,
                &mut acc,
            )?;
            for qi in &mut acc.query_id {
                *qi += start as u32;
            }
            Ok(acc)
        })
        .collect();

    let accs = results?;

    let total: usize = accs.iter().map(|a| a.query_id.len()).sum();
    let mut merged = HitAccumulator::new();
    merged.query_id.reserve(total);
    merged.position.reserve(total);
    merged.strand.reserve(total);
    merged.score.reserve(total);
    for acc in accs {
        merged.query_id.extend(acc.query_id);
        merged.position.extend(acc.position);
        merged.strand.extend(acc.strand);
        merged.score.extend(acc.score);
    }

    Ok(merged)
}

// ─── Query preparation (no PyO3) ───────────────────────────────────────────

/// Validate and prepare forward + revcomp query sequences.
///
/// All queries must be the same length and non-empty.
/// Returns `(query_len, queries_fwd, queries_rc)`.
pub fn prepare_queries(
    queries: &[String],
    mismatches: u8,
) -> Result<(usize, Vec<Vec<u8>>, Vec<Vec<u8>>), String> {
    if mismatches > 3 {
        return Err("mismatches must be 0, 1, 2, or 3".into());
    }
    if queries.is_empty() {
        return Err("queries must not be empty".into());
    }

    let query_len = queries[0].len();
    if query_len == 0 {
        return Err("query length must be > 0".into());
    }
    for (i, q) in queries.iter().enumerate() {
        if q.len() != query_len {
            return Err(format!(
                "query {} has length {} but expected {}",
                i,
                q.len(),
                query_len,
            ));
        }
    }

    let queries_fwd: Vec<Vec<u8>> = queries
        .iter()
        .map(|q| q.bytes().map(|b| b.to_ascii_uppercase()).collect())
        .collect();
    let queries_rc: Vec<Vec<u8>> = queries_fwd.iter().map(|fwd| dna::revcomp(fwd)).collect();

    Ok((query_len, queries_fwd, queries_rc))
}

// ─── Seed tier construction ─────────────────────────────────────────────────

/// Build both seed tiers (K=10 and K=14) from an FM-Index.
///
/// When the suffix array is available as a contiguous slice (i.e. built from
/// FASTA, not loaded from mmap), uses a single O(N) SA sweep per K value.
/// This eliminates all backward searches and HashMap allocations.
pub fn build_seed_tiers<I: FmOcc>(
    index: &I,
    text: &[u8],
    fasta_path: &str,
) -> Result<(SeedTier, SeedTier), String> {
    use engine::kmer_index::sa_sweep_build;
    use std::time::Instant;

    if let Some(sa) = index.sa_slice() {
        // ── SA Sweep path: O(N) per K, no backward searches, no HashMap ──
        let t0 = Instant::now();
        let (seed_small, pos_small) = sa_sweep_build(sa, text, SEED_K_SMALL);
        eprintln!(
            "[seed] K={} SA sweep (KmerSeedTable + PosTable): {:.3}s",
            SEED_K_SMALL,
            t0.elapsed().as_secs_f64()
        );

        let t1 = Instant::now();
        let (seed_large, pos_large) = sa_sweep_build(sa, text, SEED_K_LARGE);
        eprintln!(
            "[seed] K={} SA sweep (KmerSeedTable + PosTable): {:.3}s",
            SEED_K_LARGE,
            t1.elapsed().as_secs_f64()
        );

        Ok((
            SeedTier { seed_table: Arc::new(seed_small), pos_table: Arc::new(pos_small) },
            SeedTier { seed_table: Arc::new(seed_large), pos_table: Arc::new(pos_large) },
        ))
    } else {
        // ── Fallback: backward search + text scan (mmap-loaded indexes) ──
        let t0 = Instant::now();
        let seed_small = KmerSeedTable::open_or_build_from_text(index, text, fasta_path, SEED_K_SMALL)
            .map_err(|e| format!("Failed to build/load K={} seed table: {}", SEED_K_SMALL, e))?;
        eprintln!("[seed] K={} KmerSeedTable: {:.3}s", SEED_K_SMALL, t0.elapsed().as_secs_f64());

        let t1 = Instant::now();
        let pos_small = PosTable::build(text, SEED_K_SMALL);
        eprintln!("[seed] K={} PosTable: {:.3}s", SEED_K_SMALL, t1.elapsed().as_secs_f64());

        let t2 = Instant::now();
        let seed_large = KmerSeedTable::open_or_build_from_text(index, text, fasta_path, SEED_K_LARGE)
            .map_err(|e| format!("Failed to build/load K={} seed table: {}", SEED_K_LARGE, e))?;
        eprintln!("[seed] K={} KmerSeedTable: {:.3}s", SEED_K_LARGE, t2.elapsed().as_secs_f64());

        let t3 = Instant::now();
        let pos_large = PosTable::build(text, SEED_K_LARGE);
        eprintln!("[seed] K={} PosTable: {:.3}s", SEED_K_LARGE, t3.elapsed().as_secs_f64());

        Ok((
            SeedTier { seed_table: Arc::new(seed_small), pos_table: Arc::new(pos_small) },
            SeedTier { seed_table: Arc::new(seed_large), pos_table: Arc::new(pos_large) },
        ))
    }
}

/// Select the appropriate seed tier based on mismatch level and query length.
///
/// Dynamic routing:
///   mm <= 2 -> K=10 (non-overlapping, seed_mm=1, 43 variants, 10 BWT steps)
///   mm = 3 -> K=14 ONLY if 2K <= L (non-overlapping), else K=10
pub fn select_tier<'a>(
    tier_small: Option<&'a SeedTier>,
    tier_large: Option<&'a SeedTier>,
    mismatches: u8,
    query_len: usize,
) -> Option<(&'a SeedTier, usize)> {
    if mismatches <= 2 {
        if let Some(tier) = tier_small {
            if query_len > SEED_K_SMALL {
                return Some((tier, SEED_K_SMALL));
            }
        }
    } else {
        if let Some(tier) = tier_large {
            if query_len >= 2 * SEED_K_LARGE {
                return Some((tier, SEED_K_LARGE));
            }
        }
    }
    if let Some(tier) = tier_small {
        if query_len > SEED_K_SMALL {
            return Some((tier, SEED_K_SMALL));
        }
    }
    None
}

/// Build a single seed tier (K=10 or K=14) from an index handle.
///
/// Uses SA sweep for `Built` handles (suffix array available in memory).
/// Falls back to backward search + text scan for `Loaded` handles.
pub fn build_seed_tier_for_handle(
    handle: &IndexHandle,
    text: &[u8],
    index_path: &str,
    k: usize,
) -> Option<SeedTier> {
    match handle {
        IndexHandle::Built(s) => {
            // SA sweep: O(N) construction from suffix array
            if let Some(sa) = s.sa_slice() {
                let (seed_table, pos_table) = engine::kmer_index::sa_sweep_build(sa, text, k);
                return Some(SeedTier {
                    seed_table: Arc::new(seed_table),
                    pos_table: Arc::new(pos_table),
                });
            }
            // Fallback (shouldn't happen for FmIndexSearcher)
            let table = KmerSeedTable::open_or_build_from_text(&**s, text, "<mmap>", k).ok()?;
            Some(SeedTier {
                seed_table: Arc::new(table),
                pos_table: Arc::new(PosTable::build(text, k)),
            })
        }
        IndexHandle::Loaded(m) => {
            let idx_path = Path::new(index_path);
            let stem = idx_path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy();
            let parent = idx_path.parent().unwrap_or_else(|| Path::new("."));
            let kmer_path = parent.join(format!("{}.{}mer.idx", stem, k));

            let table = if kmer_path.exists() {
                KmerSeedTable::open(&kmer_path, k).ok()?
            } else {
                engine::kmer_index::build_kmer_index_from_text(&**m, text, k, &kmer_path)
                    .ok()?;
                KmerSeedTable::open(&kmer_path, k).ok()?
            };
            Some(SeedTier {
                seed_table: Arc::new(table),
                pos_table: Arc::new(PosTable::build(text, k)),
            })
        }
    }
}
