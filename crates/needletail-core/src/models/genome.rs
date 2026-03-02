//! Genome — named collection of sequences, features, and topologies.

use std::collections::HashMap;

use super::region::Region;

/// Chromosome topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    Linear,
    Circular,
}

impl Topology {
    pub fn is_circular(self) -> bool {
        self == Topology::Circular
    }

    pub fn from_str_topo(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "circular" => Topology::Circular,
            _ => Topology::Linear,
        }
    }
}

/// Metadata for a single chromosome in the master buffer.
#[derive(Debug, Clone)]
pub struct ChromMeta {
    pub name: String,
    /// Byte offset into `Genome.text` where this chromosome starts.
    pub start: usize,
    /// Sequence length (excluding the '$' sentinel).
    pub len: usize,
}

/// A genome: master sequence buffer, features, and metadata for one or more chromosomes.
///
/// The `text` buffer contains all chromosome sequences concatenated with '$' sentinels:
/// `chr1_bytes $ chr2_bytes $ ...` — uppercase, FM-Index-ready.
#[derive(Debug, Clone)]
pub struct Genome {
    pub name: String,
    /// Master sequence buffer: chr1_bytes $ chr2_bytes $ ...
    /// Uppercase, sentinel-delimited, FM-Index-ready.
    pub text: Vec<u8>,
    /// Chromosome metadata in parse order.
    pub chromosomes: Vec<ChromMeta>,
    /// All features across all chromosomes.
    pub features: Vec<Region>,
    /// Chromosome name → topology.
    pub topologies: HashMap<String, Topology>,
}

impl Genome {
    pub fn new(name: impl Into<String>) -> Self {
        Genome {
            name: name.into(),
            text: Vec::new(),
            chromosomes: Vec::new(),
            features: Vec::new(),
            topologies: HashMap::new(),
        }
    }

    /// Append a chromosome's sequence directly into the master buffer.
    ///
    /// Uppercases the sequence, appends a '$' sentinel, and records metadata.
    pub fn push_sequence(&mut self, name: String, seq: &[u8], topology: Topology) {
        let start = self.text.len();
        self.text.reserve(seq.len() + 1);
        for &b in seq {
            self.text.push(b.to_ascii_uppercase());
        }
        let len = self.text.len() - start;
        self.text.push(b'$');
        self.chromosomes.push(ChromMeta { name: name.clone(), start, len });
        self.topologies.insert(name, topology);
    }

    /// Chromosome lengths in insertion order.
    pub fn chrom_lengths(&self) -> Vec<(&str, usize)> {
        self.chromosomes
            .iter()
            .map(|cm| (cm.name.as_str(), cm.len))
            .collect()
    }

    /// Chromosome lengths as a HashMap.
    pub fn chrom_length_map(&self) -> HashMap<&str, usize> {
        self.chromosomes
            .iter()
            .map(|cm| (cm.name.as_str(), cm.len))
            .collect()
    }

    /// Sorted list of chromosome names.
    pub fn chroms(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.chromosomes.iter().map(|cm| cm.name.as_str()).collect();
        names.sort();
        names
    }

    /// Filter features to genes only (tags["feature_type"] == "gene").
    pub fn genes(&self) -> Vec<&Region> {
        self.features
            .iter()
            .filter(|r| {
                r.tags
                    .get("feature_type")
                    .and_then(|v| v.as_str())
                    .map_or(false, |t| t == "gene")
            })
            .collect()
    }

    /// Filter features to a single chromosome.
    pub fn features_on(&self, chrom: &str) -> Vec<&Region> {
        self.features.iter().filter(|r| r.chrom == chrom).collect()
    }

    /// Check if a chromosome has circular topology.
    pub fn is_circular(&self, chrom: &str) -> bool {
        self.topologies
            .get(chrom)
            .map_or(false, |t| t.is_circular())
    }

    /// Get sequence for a chromosome (slice into the master buffer).
    pub fn sequence(&self, chrom: &str) -> Option<&[u8]> {
        self.chromosomes
            .iter()
            .find(|cm| cm.name == chrom)
            .map(|cm| &self.text[cm.start..cm.start + cm.len])
    }

    /// Boolean topology vector in insertion order (for compatibility with engine).
    pub fn topology_vec(&self) -> Vec<bool> {
        self.chromosomes
            .iter()
            .map(|cm| self.is_circular(&cm.name))
            .collect()
    }

    /// Number of chromosomes.
    pub fn chrom_count(&self) -> usize {
        self.chromosomes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::{Region, Strand};

    #[test]
    fn test_genome_basics() {
        let mut g = Genome::new("test");
        g.push_sequence("chr1".into(), b"ATCGATCG", Topology::Linear);
        g.push_sequence("chr2".into(), b"GGCC", Topology::Circular);

        assert_eq!(g.chrom_lengths(), vec![("chr1", 8), ("chr2", 4)]);
        assert!(!g.is_circular("chr1"));
        assert!(g.is_circular("chr2"));
        assert_eq!(g.sequence("chr1"), Some(b"ATCGATCG".as_slice()));
        assert_eq!(g.sequence("chr2"), Some(b"GGCC".as_slice()));
        assert_eq!(g.chrom_count(), 2);

        // Master buffer: ATCGATCG$GGCC$
        assert_eq!(g.text.len(), 8 + 1 + 4 + 1);
        assert_eq!(&g.text[..9], b"ATCGATCG$");
    }

    #[test]
    fn test_genome_genes() {
        let mut g = Genome::new("test");
        g.push_sequence("chr1".into(), b"ATCGATCG", Topology::Linear);
        g.features.push(
            Region::new("chr1", 0, 100)
                .with_strand(Strand::Forward)
                .with_tag("feature_type", "gene"),
        );
        g.features.push(
            Region::new("chr1", 200, 300)
                .with_strand(Strand::Forward)
                .with_tag("feature_type", "CDS"),
        );

        let genes = g.genes();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].start, 0);
    }

    #[test]
    fn test_push_sequence_uppercases() {
        let mut g = Genome::new("test");
        g.push_sequence("chr1".into(), b"acgtACGT", Topology::Linear);
        assert_eq!(g.sequence("chr1"), Some(b"ACGTACGT".as_slice()));
    }
}
