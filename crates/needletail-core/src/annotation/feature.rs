//! Feature tiling — tile genome into promoter/gene_body/terminator/intergenic
//! regions using a sliding window over sorted genes.
//!
//! Port of: `seqchain.operations.annotate.feature.annotate_features`

use std::collections::HashMap;

use crate::geometry::resolve_landmark;
use crate::models::preset::{Anchor, FeatureConfig, FeatureDefinition};
use crate::models::region::{Region, Strand, TagValue};

/// Tile a genome into feature-type regions using gene locations and a config.
///
/// Genes must be sorted by (chrom, start). Returns feature regions covering
/// the entire genome: promoter, gene_body, terminator, and intergenic tiles.
///
/// Port of: `seqchain.operations.annotate.feature.annotate_features`
pub fn annotate_features(
    genes: &[Region],
    config: &FeatureConfig,
    chrom_sizes: &HashMap<&str, usize>,
) -> Vec<Region> {
    let mut result = Vec::new();

    // Group genes by chromosome
    let mut genes_by_chrom: HashMap<&str, Vec<&Region>> = HashMap::new();
    for gene in genes {
        genes_by_chrom
            .entry(gene.chrom.as_str())
            .or_default()
            .push(gene);
    }

    // Process each chromosome
    for (&chrom, &chrom_len) in chrom_sizes {
        let chrom_len = chrom_len as i64;

        let chrom_genes = match genes_by_chrom.get(chrom) {
            Some(g) => {
                let mut sorted = g.clone();
                sorted.sort_by_key(|r| r.start);
                sorted
            }
            None => {
                // No genes on this chromosome: entire chromosome is intergenic
                result.push(
                    Region::new(chrom, 0, chrom_len)
                        .with_name(&config.default_feature)
                        .with_tag("feature_type", config.default_feature.clone()),
                );
                continue;
            }
        };

        // Leading intergenic region (before first gene)
        if chrom_genes[0].start > 0 {
            let gap_regions = fill_gap(
                chrom,
                0,
                chrom_genes[0].start,
                None,
                Some(chrom_genes[0]),
                config,
            );
            result.extend(gap_regions);
        }

        // Sliding window over gene pairs
        for i in 0..chrom_genes.len() {
            let gene = chrom_genes[i];

            // Gene body: find overlap feature definition
            let overlap_def = config
                .features
                .iter()
                .find(|f| f.relation == "overlap");

            if let Some(def) = overlap_def {
                let mut body = Region::new(chrom, gene.start, gene.end)
                    .with_strand(gene.strand)
                    .with_name(&def.name)
                    .with_tag("feature_type", def.name.clone());

                try_enrich(&mut body, gene, def);
                result.push(body);
            }

            // Gap between this gene and the next
            if i + 1 < chrom_genes.len() {
                let next_gene = chrom_genes[i + 1];
                let gap_start = gene.end;
                let gap_end = next_gene.start;

                if gap_end > gap_start {
                    let gap_regions = fill_gap(
                        chrom,
                        gap_start,
                        gap_end,
                        Some(gene),
                        Some(next_gene),
                        config,
                    );
                    result.extend(gap_regions);
                }
            }
        }

        // Trailing region (after last gene)
        let last_gene = chrom_genes.last().unwrap();
        if last_gene.end < chrom_len {
            let gap_regions = fill_gap(
                chrom,
                last_gene.end,
                chrom_len,
                Some(last_gene),
                None,
                config,
            );
            result.extend(gap_regions);
        }
    }

    // Sort result by (chrom, start) for downstream sweep-line
    result.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    result
}

/// Fill a gap between two genes (or chromosome boundary) with regulatory
/// and intergenic regions.
///
/// Port of: `seqchain.operations.annotate.feature._fill_gap`
fn fill_gap(
    chrom: &str,
    gap_start: i64,
    gap_end: i64,
    left_gene: Option<&Region>,
    right_gene: Option<&Region>,
    config: &FeatureConfig,
) -> Vec<Region> {
    let gap_len = gap_end - gap_start;
    if gap_len <= 0 {
        return vec![];
    }

    // Collect regulatory claims: (priority, start, end, def_index, source_gene)
    let mut claims: Vec<(i32, i64, i64, usize, &Region)> = Vec::new();

    for (di, def) in config.features.iter().enumerate() {
        if def.relation == "overlap" {
            continue; // Gene bodies handled separately
        }

        // Left gene's downstream / right gene's upstream
        if let Some(lg) = left_gene {
            let is_downstream_of_left = match def.relation.as_str() {
                "downstream" => lg.strand.is_forward(),
                "upstream" => !lg.strand.is_forward(),
                _ => false,
            };
            if is_downstream_of_left && def.max_distance > 0 {
                let claim_end = (gap_start + def.max_distance).min(gap_end);
                if claim_end > gap_start {
                    claims.push((def.priority, gap_start, claim_end, di, lg));
                }
            }
        }

        if let Some(rg) = right_gene {
            let is_upstream_of_right = match def.relation.as_str() {
                "upstream" => rg.strand.is_forward(),
                "downstream" => !rg.strand.is_forward(),
                _ => false,
            };
            if is_upstream_of_right && def.max_distance > 0 {
                let claim_start = (gap_end - def.max_distance).max(gap_start);
                if claim_start < gap_end {
                    claims.push((def.priority, claim_start, gap_end, di, rg));
                }
            }
        }
    }

    if claims.is_empty() {
        // Entire gap is intergenic
        return vec![
            Region::new(chrom, gap_start, gap_end)
                .with_name(&config.default_feature)
                .with_tag("feature_type", config.default_feature.clone()),
        ];
    }

    // Sort claims by priority (lower = higher priority)
    claims.sort_by_key(|c| c.0);

    // Resolve overlaps: higher priority claims win
    // Use a coverage array approach for simplicity
    let mut assigned: Vec<Option<usize>> = vec![None; gap_len as usize];

    for (_, claim_start, claim_end, idx, _) in &claims {
        let local_start = (*claim_start - gap_start) as usize;
        let local_end = (*claim_end - gap_start) as usize;
        for pos in local_start..local_end {
            if assigned[pos].is_none() {
                assigned[pos] = Some(*idx);
            }
        }
    }

    // Build claim lookup for gene association
    let claim_genes: HashMap<(usize, i64, i64), &Region> = claims
        .iter()
        .map(|(_, s, e, idx, gene)| ((*idx, *s, *e), *gene))
        .collect();

    // Merge adjacent positions with same assignment into regions
    let mut regions = Vec::new();
    let mut run_start = 0usize;
    let mut run_assign = assigned[0];

    for i in 1..assigned.len() {
        if assigned[i] != run_assign {
            emit_gap_region(
                chrom,
                gap_start + run_start as i64,
                gap_start + i as i64,
                run_assign,
                config,
                &claims,
                &claim_genes,
                &mut regions,
            );
            run_start = i;
            run_assign = assigned[i];
        }
    }
    // Final run
    emit_gap_region(
        chrom,
        gap_start + run_start as i64,
        gap_end,
        run_assign,
        config,
        &claims,
        &claim_genes,
        &mut regions,
    );

    regions
}

/// Emit a region for a contiguous run in the gap.
#[allow(clippy::too_many_arguments)]
fn emit_gap_region(
    chrom: &str,
    start: i64,
    end: i64,
    assignment: Option<usize>,
    config: &FeatureConfig,
    claims: &[(i32, i64, i64, usize, &Region)],
    claim_genes: &HashMap<(usize, i64, i64), &Region>,
    out: &mut Vec<Region>,
) {
    match assignment {
        None => {
            out.push(
                Region::new(chrom, start, end)
                    .with_name(&config.default_feature)
                    .with_tag("feature_type", config.default_feature.clone()),
            );
        }
        Some(def_idx) => {
            let def = &config.features[def_idx];

            // Find the gene associated with this claim
            let source_gene = claims
                .iter()
                .find(|(_, cs, ce, idx, _)| {
                    *idx == def_idx && *cs <= start && *ce >= end
                })
                .or_else(|| {
                    // Fallback: any claim with this def_idx
                    claims.iter().find(|(_, _, _, idx, _)| *idx == def_idx)
                })
                .map(|(_, cs, ce, idx, gene)| {
                    claim_genes.get(&(*idx, *cs, *ce)).copied().unwrap_or(gene)
                });

            let mut region = Region::new(chrom, start, end)
                .with_name(&def.name)
                .with_tag("feature_type", def.name.clone());

            if let Some(gene) = source_gene {
                region = region.with_strand(gene.strand);
                try_enrich(&mut region, gene, def);
            }

            out.push(region);
        }
    }
}

/// Attempt to enrich a region with landmark coordinates from its source gene.
///
/// Port of: `seqchain.operations.annotate.feature._try_enrich`
fn try_enrich(region: &mut Region, gene: &Region, def: &FeatureDefinition) {
    if def.anchor == Anchor::None {
        return;
    }
    if gene.strand == Strand::Unstranded {
        return;
    }

    let fwd = gene.strand.is_forward();
    let landmark = resolve_landmark(def.anchor, gene.start, gene.end, fwd);

    region.tags.insert(
        "landmark".into(),
        TagValue::Int(landmark),
    );
    region.tags.insert(
        "gene_strand".into(),
        TagValue::Str(gene.strand.as_str().into()),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_config() -> FeatureConfig {
        FeatureConfig::by_name("saccer3").unwrap()
    }

    #[test]
    fn test_annotate_features_basic() {
        let genes = vec![
            Region::new("chr1", 1000, 2000)
                .with_strand(Strand::Forward)
                .with_name("gene1")
                .with_tag("feature_type", "gene"),
            Region::new("chr1", 3000, 4000)
                .with_strand(Strand::Reverse)
                .with_name("gene2")
                .with_tag("feature_type", "gene"),
        ];
        let config = make_config();
        let mut sizes = HashMap::new();
        sizes.insert("chr1", 5000usize);

        let tiles = annotate_features(&genes, &config, &sizes);

        // Should have tiles covering the whole chromosome
        assert!(!tiles.is_empty());

        // Check that gene bodies are present
        let gene_bodies: Vec<_> = tiles
            .iter()
            .filter(|r| {
                r.tags
                    .get("feature_type")
                    .and_then(|v| v.as_str())
                    .map_or(false, |t| t == "gene_body")
            })
            .collect();
        assert_eq!(gene_bodies.len(), 2);
        assert_eq!(gene_bodies[0].start, 1000);
        assert_eq!(gene_bodies[0].end, 2000);

        // Check promoter presence
        let promoters: Vec<_> = tiles
            .iter()
            .filter(|r| {
                r.tags
                    .get("feature_type")
                    .and_then(|v| v.as_str())
                    .map_or(false, |t| t == "promoter")
            })
            .collect();
        assert!(!promoters.is_empty());
    }

    #[test]
    fn test_annotate_features_no_genes() {
        let config = make_config();
        let mut sizes = HashMap::new();
        sizes.insert("chr1", 5000usize);

        let tiles = annotate_features(&[], &config, &sizes);
        assert_eq!(tiles.len(), 1);
        assert_eq!(
            tiles[0].tags["feature_type"].as_str(),
            Some("intergenic")
        );
        assert_eq!(tiles[0].start, 0);
        assert_eq!(tiles[0].end, 5000);
    }

    #[test]
    fn test_annotate_features_landmarks() {
        let genes = vec![
            Region::new("chr1", 500, 2000)
                .with_strand(Strand::Forward)
                .with_name("gene1")
                .with_tag("feature_type", "gene"),
        ];
        let config = make_config();
        let mut sizes = HashMap::new();
        sizes.insert("chr1", 3000usize);

        let tiles = annotate_features(&genes, &config, &sizes);

        // Gene body should have five_prime landmark = gene start (forward)
        let gene_body = tiles.iter().find(|r| {
            r.tags.get("feature_type").and_then(|v| v.as_str()) == Some("gene_body")
        }).unwrap();
        assert_eq!(gene_body.tags["landmark"].as_i64(), Some(500));

        // Promoter should have five_prime landmark = gene start
        let promoter = tiles.iter().find(|r| {
            r.tags.get("feature_type").and_then(|v| v.as_str()) == Some("promoter")
        });
        if let Some(p) = promoter {
            assert_eq!(p.tags["landmark"].as_i64(), Some(500));
        }
    }
}
