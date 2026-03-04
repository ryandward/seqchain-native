//! Region — the universal genomic interval.
//!
//! A Region is an oriented interval on a named 1D coordinate space, plus
//! a metadata map. It is the common currency of every pipeline stage:
//! genes, guides, features, peaks — all are Regions.

use std::collections::HashMap;

/// Strand orientation in a 1D coordinate space.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
    Unstranded,
}

impl Strand {
    /// Forward = true, Reverse = false, Unstranded = true (treated as forward).
    #[inline(always)]
    pub fn is_forward(self) -> bool {
        self != Strand::Reverse
    }

    /// Strand sign: +1 for Forward/Unstranded, -1 for Reverse.
    #[inline(always)]
    pub fn sign(self) -> i64 {
        if self == Strand::Reverse { -1 } else { 1 }
    }

    /// Parse from string: "+", "-", or ".".
    pub fn from_str_strand(s: &str) -> Self {
        match s {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unstranded,
        }
    }

    /// Convert to canonical string representation.
    pub fn as_str(self) -> &'static str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unstranded => ".",
        }
    }
}

/// Dynamically-typed annotation value for Region tags.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    Int(i64),
    Float(f64),
    Str(String),
    Bool(bool),
}

impl TagValue {
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            TagValue::Int(v) => Some(*v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            TagValue::Float(v) => Some(*v),
            TagValue::Int(v) => Some(*v as f64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            TagValue::Str(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            TagValue::Bool(v) => Some(*v),
            _ => None,
        }
    }
}

impl From<i64> for TagValue {
    fn from(v: i64) -> Self { TagValue::Int(v) }
}

impl From<f64> for TagValue {
    fn from(v: f64) -> Self { TagValue::Float(v) }
}

impl From<String> for TagValue {
    fn from(v: String) -> Self { TagValue::Str(v) }
}

impl From<&str> for TagValue {
    fn from(v: &str) -> Self { TagValue::Str(v.to_string()) }
}

impl From<bool> for TagValue {
    fn from(v: bool) -> Self { TagValue::Bool(v) }
}

/// A genomic interval with orientation and metadata.
///
/// Coordinates are 0-based, half-open: `[start, end)`.
/// Tags carry domain-specific metadata (feature_type, gene_name, etc.)
/// without encoding it in the type system.
#[derive(Debug, Clone)]
pub struct Region {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub score: Option<f64>,
    pub name: String,
    pub tags: HashMap<String, TagValue>,
}

impl Region {
    /// Create a new Region with minimal fields; tags default to empty.
    pub fn new(chrom: impl Into<String>, start: i64, end: i64) -> Self {
        Region {
            chrom: chrom.into(),
            start,
            end,
            strand: Strand::Unstranded,
            score: None,
            name: String::new(),
            tags: HashMap::new(),
        }
    }

    /// Builder: set strand.
    pub fn with_strand(mut self, strand: Strand) -> Self {
        self.strand = strand;
        self
    }

    /// Builder: set name.
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Builder: set score.
    pub fn with_score(mut self, score: f64) -> Self {
        self.score = Some(score);
        self
    }

    /// Builder: insert a tag.
    pub fn with_tag(mut self, key: impl Into<String>, value: impl Into<TagValue>) -> Self {
        self.tags.insert(key.into(), value.into());
        self
    }

    /// Length in bases: `end - start`.
    #[inline(always)]
    pub fn length(&self) -> i64 {
        self.end - self.start
    }

    /// Center position: `(start + end) / 2` (integer division).
    #[inline(always)]
    pub fn center(&self) -> i64 {
        (self.start + self.end) / 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_region_basics() {
        let r = Region::new("chr1", 100, 200)
            .with_strand(Strand::Forward)
            .with_name("test_gene")
            .with_score(1.5)
            .with_tag("feature_type", "gene");

        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.length(), 100);
        assert_eq!(r.center(), 150);
        assert_eq!(r.strand, Strand::Forward);
        assert_eq!(r.name, "test_gene");
        assert_eq!(r.score, Some(1.5));
        assert_eq!(r.tags["feature_type"].as_str(), Some("gene"));
    }

    #[test]
    fn test_strand_properties() {
        assert!(Strand::Forward.is_forward());
        assert!(!Strand::Reverse.is_forward());
        assert!(Strand::Unstranded.is_forward());
        assert_eq!(Strand::Forward.sign(), 1);
        assert_eq!(Strand::Reverse.sign(), -1);
    }

    #[test]
    fn test_tag_value_conversions() {
        let v: TagValue = 42i64.into();
        assert_eq!(v.as_i64(), Some(42));
        assert_eq!(v.as_f64(), Some(42.0));

        let v: TagValue = 3.14f64.into();
        assert_eq!(v.as_f64(), Some(3.14));

        let v: TagValue = "hello".into();
        assert_eq!(v.as_str(), Some("hello"));

        let v: TagValue = true.into();
        assert_eq!(v.as_bool(), Some(true));
    }
}
