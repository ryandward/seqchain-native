//! Per-locus annotation — compute overlap, offset, and relative position
//! for a Region against a set of overlapping feature Regions.
//!
//! Port of: `seqchain.primitives.annotate.annotate_locus_from_features`

use crate::geometry::{interval_overlap, offset_in_feature, relative_position};
use crate::models::region::{Region, TagValue};

/// Annotate a single region against pre-filtered overlapping features.
///
/// For each overlapping feature, produces a new Region with annotation tags:
/// - `feature_name`, `feature_type`, `feature_start`, `feature_end`, `feature_strand`
/// - `overlap` (bp), `offset` (strand-aware), `relative_pos` (0.0–1.0)
/// - Copies optional tags: `locus_tag`, `gene`, `product`
///
/// If no features overlap, returns the region unchanged.
pub fn annotate_locus_from_features(
    region: &Region,
    features: &[&Region],
    chrom_len: Option<i64>,
) -> Vec<Region> {
    if features.is_empty() {
        return vec![region.clone()];
    }

    let mut results = Vec::with_capacity(features.len());

    for feat in features {
        let mut annotated = region.clone();
        annotate_locus_in_place(&mut annotated, feat, chrom_len);
        results.push(annotated);
    }

    results
}

/// Annotate a region in-place against a single feature.
///
/// Adds annotation tags directly to the region — no clone, no allocation
/// beyond the tag insertions themselves.
pub fn annotate_locus_in_place(
    region: &mut Region,
    feat: &Region,
    chrom_len: Option<i64>,
) {
    let overlap_bp = interval_overlap(
        region.start,
        region.end,
        feat.start,
        feat.end,
        chrom_len,
    );

    let feat_fwd = feat.strand.is_forward();
    let cl = chrom_len.unwrap_or(0);

    let offset = offset_in_feature(
        region.start,
        region.end,
        feat.start,
        feat.end,
        feat_fwd,
        cl,
    );

    let rel_pos = relative_position(region.start, feat.start, feat.end, cl);

    region.tags.insert(
        "feature_name".into(),
        TagValue::Str(feat.name.clone()),
    );
    if let Some(ft) = feat.tags.get("feature_type") {
        region.tags.insert("feature_type".into(), ft.clone());
    }
    region.tags.insert("feature_start".into(), TagValue::Int(feat.start));
    region.tags.insert("feature_end".into(), TagValue::Int(feat.end));
    region.tags.insert(
        "feature_strand".into(),
        TagValue::Str(feat.strand.as_str().into()),
    );
    region.tags.insert("overlap".into(), TagValue::Int(overlap_bp));
    region.tags.insert("offset".into(), TagValue::Int(offset));
    region.tags.insert("relative_pos".into(), TagValue::Float(rel_pos));

    // Copy optional feature tags
    for key in &["locus_tag", "gene", "product", "landmark", "gene_strand"] {
        if let Some(val) = feat.tags.get(*key) {
            region.tags.insert((*key).to_string(), val.clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_annotate_locus_basic() {
        let guide = Region::new("chr1", 150, 170)
            .with_strand(Strand::Forward)
            .with_name("guide1");

        let gene = Region::new("chr1", 100, 500)
            .with_strand(Strand::Forward)
            .with_name("YAL001C")
            .with_tag("feature_type", "promoter")
            .with_tag("locus_tag", "YAL001C");

        let results = annotate_locus_from_features(&guide, &[&gene], None);
        assert_eq!(results.len(), 1);

        let r = &results[0];
        assert_eq!(r.tags["overlap"].as_i64(), Some(20)); // 170-150
        assert_eq!(r.tags["offset"].as_i64(), Some(50)); // 150-100
        assert_eq!(r.tags["feature_name"].as_str(), Some("YAL001C"));
        assert_eq!(r.tags["locus_tag"].as_str(), Some("YAL001C"));
    }

    #[test]
    fn test_annotate_locus_no_features() {
        let guide = Region::new("chr1", 150, 170);
        let results = annotate_locus_from_features(&guide, &[], None);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 150);
    }
}
