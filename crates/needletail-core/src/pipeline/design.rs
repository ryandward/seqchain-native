//! Pipeline orchestrator — wires all stages into `design_library()`.
//!
//! This replaces SeqChain's `library.py::run_job()`. It is the only module
//! that imports from all other layers.
//!
//! Stages:
//! 1. Feature tiling (annotation/feature)
//! 2. PAM scanning (operations/pam_scanner)
//! 3. Guide enrichment (operations/pam_scanner)
//! 4. Off-target scoring (engine/)
//! 5. Sweep-line annotation → terminal sink (zero additional RAM)

use std::collections::HashMap;

use crate::annotation::feature::annotate_features;
use crate::annotation::sweep::SweepAnnotator;
use crate::chemistry::{CompiledPam, PamDirection};
use crate::geometry::signed_distance;
use crate::io::RegionSink;
use crate::models::genome::Genome;
use crate::models::preset::{CRISPRPreset, FeatureConfig};
use crate::models::region::{Region, Strand, TagValue};
use crate::operations::pam_scanner::{enrich_hits, find_pam_sites};
use crate::{
    filter_hits_by_pam, prepare_queries, run_search_seeded, run_search_unseeded,
    select_tier, IndexHandle, SeedTier,
};

/// Progress callback for reporting pipeline stage advancement.
///
/// The step/stage/items model:
/// - **step**: UI step key (matches method schema `steps[].key`)
/// - **stage**: free-text sublabel within the current step
/// - **items**: discrete n-of-x counter for countable work
///
/// `pct_complete` is per-step: it resets when step changes.
pub trait ProgressSink: Send + Sync {
    fn report(&self, stage: &str, current: usize, total: usize);

    /// Set the current UI step key (e.g., "scanning", "scoring").
    fn set_step(&self, _step: &str) {}

    /// Set a free-text sublabel within the current step.
    fn set_stage(&self, _stage: &str) {}

    /// Report item-level progress within a step.
    fn set_items(&self, _complete: usize, _total: usize) {}

    /// Clear item counter (when entering a step with no countable items).
    fn clear_items(&self) {}

    fn is_cancelled(&self) -> bool { false }
}

/// Null progress sink — does nothing.
pub struct NullProgress;
impl ProgressSink for NullProgress {
    fn report(&self, _stage: &str, _current: usize, _total: usize) {}
}

/// Result of the design_library pipeline.
///
/// Guides are drained directly into the `RegionSink` — they are NOT held
/// in memory.  This struct carries only metadata.
pub struct LibraryResult {
    /// Total guides scored (before annotation).
    pub total_guides_scored: usize,
    /// Total guides written to the sink.
    pub guides_written: usize,
}

/// Design a whole-genome CRISPR guide library: the full Rust pipeline.
///
/// This is the equivalent of SeqChain's `run_job()`:
/// 1. Tile genome into feature regions
/// 2. Scan for PAM sites → guides
/// 3. Score off-targets
/// 4. Sweep-line annotate against features → drain into sink
///
/// Guides are streamed through the sink one at a time.  The annotation
/// phase runs with O(k) auxiliary memory, regardless of genome size.
pub fn design_library(
    genome: &Genome,
    index: &IndexHandle,
    tier_small: Option<&SeedTier>,
    tier_large: Option<&SeedTier>,
    preset: &CRISPRPreset,
    feature_config: &FeatureConfig,
    progress: &dyn ProgressSink,
    sink: &mut dyn RegionSink,
) -> Result<LibraryResult, String> {
    fn peak_mb() -> f64 {
        std::fs::read_to_string("/proc/self/status").ok()
            .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse::<f64>().ok()))
            .unwrap_or(0.0) / 1024.0
    }
    fn rss_mb() -> f64 {
        std::fs::read_to_string("/proc/self/status").ok()
            .and_then(|s| s.lines().find(|l| l.starts_with("VmRSS:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse::<f64>().ok()))
            .unwrap_or(0.0) / 1024.0
    }
    eprintln!("[mem] pipeline start: RSS={:.0}MB peak={:.0}MB", rss_mb(), peak_mb());

    // ── Step: indexing (brief — index was built at upload time) ─────────────
    progress.set_step("indexing");
    progress.set_stage("Preparing genome data");
    progress.clear_items();
    progress.report("indexing", 0, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    let genes: Vec<Region> = genome.genes().into_iter().cloned().collect();
    let chrom_sizes = genome.chrom_length_map();
    // ── Step: scanning ──────────────────────────────────────────────────────
    progress.set_step("scanning");
    progress.set_stage("Tiling genome features");
    progress.clear_items();
    progress.report("tiling", 1, 7);

    let feature_tiles = annotate_features(&genes, feature_config, &chrom_sizes);

    progress.set_stage("Scanning PAM sites");
    progress.report("scanning", 2, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    let direction = match preset.pam_direction.as_str() {
        "upstream" => PamDirection::Upstream,
        _ => PamDirection::Downstream,
    };

    let chroms = index.chrom_geometry();
    let chrom_names = index.chrom_names();
    let topologies = genome.topology_vec();
    let topo_ref = Some(topologies.as_slice());

    let pam_hits = find_pam_sites(
        index.text(),
        &chroms,
        &preset.pam,
        preset.spacer_len,
        direction,
        topo_ref,
    )?;
    eprintln!("[mem] after PAM scan: RSS={:.0}MB peak={:.0}MB", rss_mb(), peak_mb());

    progress.set_stage("Enriching guide hits");
    progress.report("enriching", 3, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    let chrom_ranges = chroms.ranges.clone();
    let guide_hits = enrich_hits(pam_hits, &chrom_names, direction, topo_ref, &chrom_ranges);
    eprintln!("[mem] after enrich ({} guides): RSS={:.0}MB peak={:.0}MB", guide_hits.count, rss_mb(), peak_mb());

    // ── Step: scoring ───────────────────────────────────────────────────────
    progress.set_step("scoring");
    progress.set_stage("Scoring off-targets");
    progress.clear_items();
    progress.report("scoring", 4, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    // Collect unique spacers for deduplication
    let sl = guide_hits.spacer_len;
    let n_guides = guide_hits.count;

    let mut unique_spacers: Vec<String> = Vec::new();
    let mut spacer_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut guide_spacer_map: Vec<usize> = Vec::with_capacity(n_guides);

    for i in 0..n_guides {
        let spacer = &guide_hits.spacers[i * sl..(i + 1) * sl];
        let spacer_idx = if let Some(&idx) = spacer_to_idx.get(spacer) {
            idx
        } else {
            let idx = unique_spacers.len();
            unique_spacers.push(
                std::str::from_utf8(spacer)
                    .unwrap_or("")
                    .to_string(),
            );
            spacer_to_idx.insert(spacer.to_vec(), idx);
            idx
        };
        guide_spacer_map.push(spacer_idx);
    }

    // Search all unique spacers
    let n_unique = unique_spacers.len();
    progress.set_items(0, n_unique);
    progress.set_stage(&format!("Aligning {} unique spacers", n_unique));
    let t_score = std::time::Instant::now();
    let hit_counts = if !unique_spacers.is_empty() {
        score_spacers(
            &unique_spacers,
            index,
            tier_small,
            tier_large,
            &preset.pam,
            direction,
            preset.mismatches,
            topo_ref,
        )?
    } else {
        vec![]
    };
    progress.set_items(n_unique, n_unique);
    eprintln!("[pipeline] score_spacers: {:.3}s ({} unique, mm={})",
        t_score.elapsed().as_secs_f64(), n_unique, preset.mismatches);

    // Drop dedup structures — no longer needed
    drop(unique_spacers);
    drop(spacer_to_idx);
    eprintln!("[mem] after scoring: RSS={:.0}MB peak={:.0}MB", rss_mb(), peak_mb());

    // ── Build scored guide Regions ────────────────────────────────────────────
    progress.set_stage("Building scored regions");
    progress.set_items(0, n_guides);
    progress.report("annotating", 5, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    let pl = guide_hits.pam_len;
    let gl = sl + pl;
    let report_interval = (n_guides / 100).max(1000);

    let mut scored_guides: Vec<Region> = Vec::with_capacity(n_guides);
    for i in 0..n_guides {
        let cid = guide_hits.chrom_ids[i] as usize;
        let chrom_name = &chrom_names[cid];
        let gs = guide_hits.guide_starts[i] as i64;
        let ge = guide_hits.guide_ends[i] as i64;
        let strand = if guide_hits.strands[i] > 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        };

        let spacer = std::str::from_utf8(&guide_hits.spacers[i * sl..(i + 1) * sl])
            .unwrap_or("")
            .to_string();
        let pam_seq = std::str::from_utf8(&guide_hits.pam_seqs[i * pl..(i + 1) * pl])
            .unwrap_or("")
            .to_string();
        let guide_seq = std::str::from_utf8(&guide_hits.guide_seqs[i * gl..(i + 1) * gl])
            .unwrap_or("")
            .to_string();
        let guide_id = std::str::from_utf8(&guide_hits.guide_ids[i * 8..(i + 1) * 8])
            .unwrap_or("")
            .to_string();

        let spacer_idx = guide_spacer_map[i];
        let total_hits = hit_counts.get(spacer_idx).copied().unwrap_or(0);
        let off_targets = if total_hits > 0 { total_hits - 1 } else { 0 };

        let r = Region::new(chrom_name, gs, ge)
            .with_strand(strand)
            .with_name(&guide_id)
            .with_score(off_targets as f64)
            .with_tag("spacer", spacer)
            .with_tag("pam_seq", pam_seq)
            .with_tag("guide_seq", guide_seq)
            .with_tag("guide_id", guide_id)
            .with_tag("off_targets", off_targets as i64)
            .with_tag("total_hits", total_hits as i64);

        scored_guides.push(r);

        if (i + 1) % report_interval == 0 {
            progress.set_items(i + 1, n_guides);
        }
    }
    progress.set_items(n_guides, n_guides);

    let total_guides_scored = scored_guides.len();
    eprintln!("[mem] after {} guide regions: RSS={:.0}MB peak={:.0}MB", total_guides_scored, rss_mb(), peak_mb());

    // Sort guides by (chrom, start) for sweep-line
    scored_guides.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    // ── Sweep-line annotation → terminal drain (O(k) auxiliary RAM) ──────────
    progress.set_step("annotation");
    progress.set_stage("Annotating & streaming to sink");
    progress.clear_items();
    progress.report("sweep", 6, 7);
    if progress.is_cancelled() {
        return Err("cancelled".into());
    }

    let chrom_len_map = genome.chrom_length_map();

    // Ownership of scored_guides is transferred via into_iter().
    // The Vec's memory is freed as the iterator advances — annotation
    // footprint is O(k) where k is the active overlap window.
    let annotator = SweepAnnotator::new(
        scored_guides.into_iter(),
        &feature_tiles,
        |chrom| chrom_len_map.get(chrom).map(|&len| len as i64),
    );

    // Terminal drain: sweep → TSS distance → sink.  No collect().
    // Report progress every ~1% so GenomeHub's polling sees a smooth bar.
    let drain_report_interval = (total_guides_scored / 100).max(1000);
    let mut guides_written: usize = 0;
    progress.set_items(0, total_guides_scored);
    for mut region in annotator {
        // Compute TSS distance for regions that have a landmark
        if let (Some(landmark), Some(feat_strand)) = (
            region.tags.get("landmark").and_then(|v| v.as_i64()),
            region
                .tags
                .get("feature_strand")
                .and_then(|v| v.as_str())
                .or_else(|| region.tags.get("gene_strand").and_then(|v| v.as_str())),
        ) {
            let fwd = feat_strand == "+";
            let dist = signed_distance(region.start, landmark, fwd);
            region
                .tags
                .insert("signed_distance".into(), TagValue::Int(dist));
        }

        sink.consume(region).map_err(|e| format!("sink error: {}", e))?;
        guides_written += 1;

        if guides_written % drain_report_interval == 0 {
            progress.set_items(guides_written, total_guides_scored);
        }
    }

    eprintln!(
        "[mem] after annotation+sink ({} guides written): RSS={:.0}MB peak={:.0}MB",
        guides_written,
        rss_mb(),
        peak_mb()
    );

    progress.set_items(guides_written, guides_written);
    progress.set_stage(&format!("Streamed {} annotated guides", guides_written));
    progress.report("complete", 7, 7);

    Ok(LibraryResult {
        total_guides_scored,
        guides_written,
    })
}

/// Score unique spacers: search + PAM filter → hit counts per spacer.
#[allow(clippy::too_many_arguments)]
fn score_spacers(
    spacers: &[String],
    index: &IndexHandle,
    tier_small: Option<&SeedTier>,
    tier_large: Option<&SeedTier>,
    pam: &str,
    direction: PamDirection,
    mismatches: u8,
    topologies: Option<&[bool]>,
) -> Result<Vec<u64>, String> {
    let (query_len, queries_fwd, queries_rc) = prepare_queries(spacers, mismatches)?;
    let chroms = index.chrom_geometry();
    let max_width = 8u32;

    let mut acc = if let Some((tier, _k)) =
        select_tier(tier_small, tier_large, mismatches, query_len)
    {
        let seed_table = tier.seed_table.clone();
        let pos_table = tier.pos_table.clone();

        match index {
            IndexHandle::Built(idx) => run_search_seeded(
                &**idx,
                &seed_table,
                &pos_table,
                idx.text(),
                &queries_fwd,
                &queries_rc,
                query_len,
                mismatches,
                max_width,
                &chroms,
            )
            .map_err(|e| e.to_string())?,
            IndexHandle::Loaded(idx) => run_search_seeded(
                &**idx,
                &seed_table,
                &pos_table,
                idx.text(),
                &queries_fwd,
                &queries_rc,
                query_len,
                mismatches,
                max_width,
                &chroms,
            )
            .map_err(|e| e.to_string())?,
        }
    } else {
        match index {
            IndexHandle::Built(idx) => run_search_unseeded(
                &**idx,
                &queries_fwd,
                &queries_rc,
                query_len,
                mismatches,
                max_width,
                &chroms,
            ),
            IndexHandle::Loaded(idx) => run_search_unseeded(
                &**idx,
                &queries_fwd,
                &queries_rc,
                query_len,
                mismatches,
                max_width,
                &chroms,
            ),
        }
    };

    // PAM filter
    let compiled = CompiledPam::compile(pam)?;
    let text = index.text();
    let filter_chroms = index.chrom_geometry();
    filter_hits_by_pam(
        &mut acc,
        text,
        &filter_chroms,
        &compiled,
        direction,
        query_len,
        topologies,
    );

    // Count hits per query
    let n = spacers.len();
    let mut counts = vec![0u64; n];
    for &qid in &acc.query_id {
        if (qid as usize) < n {
            counts[qid as usize] += 1;
        }
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_null_progress() {
        let p = NullProgress;
        p.report("test", 0, 1);
        assert!(!p.is_cancelled());
    }
}
