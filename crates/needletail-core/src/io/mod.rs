//! The Disk Boundary — serialization, mmap, and columnar export.
//!
//! This layer handles bytes touching the disk: rkyv-serialized FM-Index
//! archives loaded via huge-page mmap, and zero-copy Arrow/Parquet sinks.
//! It knows how to persist and resurrect engine data structures but has
//! no concept of biology, chemistry, or Python.
//!
//! Import hierarchy: `engine`, `error` (never chemistry, operations, or lib).

pub mod fastq;
pub mod genbank;
pub mod json;
pub mod parquet_hits;
pub mod parquet_regions;
pub mod persist;
pub mod sam;

use crate::models::region::Region;

/// Tags that are internal pipeline state — never exported to output files.
///
/// Both JSON and Parquet sinks filter these out.  Add new internal tags
/// here — this is the single source of truth.
pub const HIDDEN_TAGS: &[&str] = &["landmark", "gene_strand"];

/// Canonical ordering for tag columns in output files.
///
/// Tags listed here appear in this sequence after the fixed Region fields
/// (chrom, start, end, strand, score).  Any tags NOT listed here follow
/// in alphabetical order.  This keeps output stable across pipeline
/// versions while remaining open to new tags without code changes.
pub const TAG_ORDER: &[&str] = &[
    "guide_id", "guide_seq", "spacer", "pam_seq",
    "off_targets", "total_hits",
    "feature_type", "feature_name", "feature_strand",
    "feature_start", "feature_end",
    "signed_distance", "relative_pos", "offset", "overlap",
];

/// Terminal sink that consumes annotated Regions one at a time.
///
/// This is the I/O boundary: the pipeline transfers ownership of each
/// Region into the sink, which writes it to disk (or buffers it).
/// By draining through a sink, the annotation phase runs with O(k)
/// auxiliary memory where k is the active overlap window size.
pub trait RegionSink {
    fn consume(&mut self, region: Region) -> std::io::Result<()>;
}

/// Counting sink — discards regions but counts them.
///
/// Useful for benchmarks and dry runs.
pub struct CountingSink {
    count: usize,
}

impl CountingSink {
    pub fn new() -> Self {
        CountingSink { count: 0 }
    }

    pub fn count(&self) -> usize {
        self.count
    }
}

impl RegionSink for CountingSink {
    fn consume(&mut self, _region: Region) -> std::io::Result<()> {
        self.count += 1;
        Ok(())
    }
}
