//! Sweep-line annotation — join guides against feature tiles using a
//! sorted-merge algorithm with an active window deque.
//!
//! `SweepAnnotator` is a lazy, ownership-transferring Iterator: it pulls
//! one guide at a time from an owning iterator, annotates it against
//! feature tiles, and yields the result.  The source Vec is destructured
//! as the sweep advances, so annotation runs with O(k) auxiliary memory
//! where k is the active overlap window size.

use std::collections::VecDeque;
use std::iter::Peekable;

use crate::models::region::Region;

use super::locus::annotate_locus_in_place;

/// Lazy, ownership-transferring sweep-line annotator.
///
/// Takes ownership of guides via `IntoIterator`, pulling one at a time.
/// Features are borrowed (they're the static reference layer).
/// The full annotated set never lives in memory simultaneously.
pub struct SweepAnnotator<'a, I: Iterator<Item = Region>, F: Fn(&str) -> Option<i64>> {
    guides: Peekable<I>,
    features: &'a [Region],
    chrom_len_fn: F,

    /// Index into `features` — next feature to consider adding to the window.
    feat_cursor: usize,
    /// Active feature window (indices into `features`).
    active: VecDeque<usize>,

    /// When a guide overlaps multiple features we buffer the extra annotated
    /// copies here and drain them before advancing to the next guide.
    pending: VecDeque<Region>,
}

impl<'a, I: Iterator<Item = Region>, F: Fn(&str) -> Option<i64>> SweepAnnotator<'a, I, F> {
    pub fn new(
        guides: impl IntoIterator<IntoIter = I>,
        features: &'a [Region],
        chrom_len_fn: F,
    ) -> Self {
        SweepAnnotator {
            guides: guides.into_iter().peekable(),
            features,
            chrom_len_fn,
            feat_cursor: 0,
            active: VecDeque::new(),
            pending: VecDeque::new(),
        }
    }
}

impl<'a, I: Iterator<Item = Region>, F: Fn(&str) -> Option<i64>> Iterator
    for SweepAnnotator<'a, I, F>
{
    type Item = Region;

    fn next(&mut self) -> Option<Region> {
        // Drain any buffered multi-overlap results first.
        if let Some(r) = self.pending.pop_front() {
            return Some(r);
        }

        // Pull the next guide from the owning iterator.
        let region = self.guides.next()?;

        // Evict: remove features entirely behind the current region.
        while let Some(&front) = self.active.front() {
            let f = &self.features[front];
            if f.chrom != region.chrom || f.end <= region.start {
                self.active.pop_front();
            } else {
                break;
            }
        }

        // Advance: load features whose start < region.end.
        while self.feat_cursor < self.features.len() {
            let f = &self.features[self.feat_cursor];
            if f.chrom > region.chrom {
                break;
            }
            if f.chrom == region.chrom && f.start >= region.end {
                break;
            }
            if f.chrom == region.chrom && f.end > region.start {
                self.active.push_back(self.feat_cursor);
            }
            self.feat_cursor += 1;
        }

        // Collect overlapping feature indices.
        let overlapping: Vec<usize> = self
            .active
            .iter()
            .copied()
            .filter(|&idx| {
                let f = &self.features[idx];
                f.chrom == region.chrom && f.start < region.end && f.end > region.start
            })
            .collect();

        let chrom_len = (self.chrom_len_fn)(&region.chrom);

        if overlapping.is_empty() {
            // No overlap — pass through unchanged.
            Some(region)
        } else {
            // First overlap → returned directly.  Extras → buffered in pending.
            let mut first = region.clone();
            annotate_locus_in_place(&mut first, &self.features[overlapping[0]], chrom_len);

            for &idx in &overlapping[1..] {
                let mut extra = region.clone();
                annotate_locus_in_place(&mut extra, &self.features[idx], chrom_len);
                self.pending.push_back(extra);
            }

            Some(first)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (lower, _) = self.guides.size_hint();
        (lower + self.pending.len(), None) // upper bound unknown due to multi-overlap expansion
    }
}

/// Convenience wrapper — returns a collecting Vec for callers that need it.
/// Prefer using `SweepAnnotator` directly for streaming composition.
pub fn sorted_overlap_annotate(
    regions: Vec<Region>,
    features: &[Region],
    chrom_len_fn: impl Fn(&str) -> Option<i64>,
) -> Vec<Region> {
    SweepAnnotator::new(regions.into_iter(), features, chrom_len_fn).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_sweep_annotate() {
        let guides = vec![
            Region::new("chr1", 50, 70).with_strand(Strand::Forward),
            Region::new("chr1", 150, 170).with_strand(Strand::Forward),
            Region::new("chr1", 550, 570).with_strand(Strand::Forward),
        ];

        let features = vec![
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_name("promoter_A")
                .with_tag("feature_type", "promoter"),
            Region::new("chr1", 100, 500)
                .with_strand(Strand::Forward)
                .with_name("gene_A")
                .with_tag("feature_type", "gene_body"),
            Region::new("chr1", 500, 600)
                .with_strand(Strand::Forward)
                .with_name("term_A")
                .with_tag("feature_type", "terminator"),
        ];

        let results = sorted_overlap_annotate(guides, &features, |_| None);

        assert_eq!(results.len(), 3);
        assert_eq!(results[0].tags["feature_name"].as_str(), Some("promoter_A"));
        assert_eq!(results[1].tags["feature_name"].as_str(), Some("gene_A"));
        assert_eq!(results[2].tags["feature_name"].as_str(), Some("term_A"));
    }

    #[test]
    fn test_sweep_no_overlap() {
        let guides = vec![
            Region::new("chr1", 1000, 1020).with_strand(Strand::Forward),
        ];
        let features = vec![
            Region::new("chr1", 0, 100)
                .with_name("feat")
                .with_tag("feature_type", "gene"),
        ];

        let results = sorted_overlap_annotate(guides, &features, |_| None);
        assert_eq!(results.len(), 1);
        assert!(!results[0].tags.contains_key("feature_name"));
    }

    #[test]
    fn test_sweep_iterator_lazy() {
        let guides = vec![
            Region::new("chr1", 50, 70).with_strand(Strand::Forward),
            Region::new("chr1", 150, 170).with_strand(Strand::Forward),
        ];

        let features = vec![
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_name("promoter_A")
                .with_tag("feature_type", "promoter"),
            Region::new("chr1", 100, 500)
                .with_strand(Strand::Forward)
                .with_name("gene_A")
                .with_tag("feature_type", "gene_body"),
        ];

        // Use as iterator — filter without collecting all.
        let promoters: Vec<Region> =
            SweepAnnotator::new(guides.into_iter(), &features, |_| None)
                .filter(|r| {
                    r.tags
                        .get("feature_type")
                        .and_then(|v| v.as_str())
                        .map_or(false, |t| t == "promoter")
                })
                .collect();

        assert_eq!(promoters.len(), 1);
        assert_eq!(
            promoters[0].tags["feature_name"].as_str(),
            Some("promoter_A")
        );
    }

    #[test]
    fn test_sweep_ownership_transfer() {
        // Verify that guides Vec is consumed — no double allocation.
        let guides = vec![
            Region::new("chr1", 50, 70).with_strand(Strand::Forward),
        ];
        let features = vec![
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_name("feat_A")
                .with_tag("feature_type", "promoter"),
        ];

        // guides is moved here — cannot be used after this line
        let mut annotator = SweepAnnotator::new(guides.into_iter(), &features, |_| None);
        let first = annotator.next().unwrap();
        assert_eq!(first.tags["feature_name"].as_str(), Some("feat_A"));
        assert!(annotator.next().is_none());
    }
}
