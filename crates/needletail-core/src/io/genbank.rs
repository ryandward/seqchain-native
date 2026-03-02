//! GenBank parser — `.gb`/`.gb.gz` files to `Genome`.
//!
//! Uses the `gb-io` crate for parsing. Extracts sequences, gene features,
//! topologies, and handles compound locations and origin-wrapping genes.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::models::genome::{Genome, Topology};
use crate::models::region::{Region, Strand, TagValue};

/// Load a GenBank file (`.gb` or `.gb.gz`) into a `Genome`.
pub fn load_genbank(path: &Path) -> Result<Genome, String> {
    if !path.exists() {
        return Err(format!("Genome file not found: {}", path.display()));
    }

    let name = path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    // Strip .gb from .gb.gz stems
    let name = name.strip_suffix(".gb").unwrap_or(&name).to_string();

    let mut genome = Genome::new(name);

    let records = parse_genbank_file(path)?;

    for seq in records {
        let chrom = seq.name.clone().unwrap_or_default();
        if chrom.is_empty() {
            continue;
        }

        let seq_len = seq.seq.len();

        // Topology
        let topology = if seq.is_circular() {
            Topology::Circular
        } else {
            Topology::Linear
        };

        // Single-pass: uppercase + push into master buffer + sentinel
        genome.push_sequence(chrom.clone(), &seq.seq, topology);

        // Features
        for feature in &seq.features {
            // Skip "source" features
            if feature.kind.as_ref() == "source" {
                continue;
            }

            let regions = extract_feature(feature, &chrom, seq_len, topology.is_circular());
            genome.features.extend(regions);
        }
    }

    Ok(genome)
}

/// Parse a GenBank file, handling gzip transparently.
fn parse_genbank_file(path: &Path) -> Result<Vec<gb_io::seq::Seq>, String> {
    let path_str = path.to_string_lossy();
    let is_gz = path_str.ends_with(".gz");

    let file = File::open(path)
        .map_err(|e| format!("Failed to open {}: {}", path.display(), e))?;

    let reader: Box<dyn Read> = if is_gz {
        Box::new(GzDecoder::new(BufReader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let seqs: Result<Vec<_>, _> = gb_io::reader::SeqReader::new(reader).collect();
    seqs.map_err(|e| format!("Failed to parse GenBank file {}: {}", path.display(), e))
}

/// Extract a Region from a gb-io Feature.
///
/// Handles simple locations, complement (reverse strand), and compound
/// locations (join). Origin-wrapping genes on circular chromosomes get
/// virtual coordinates beyond `seq_len`.
fn extract_feature(
    feature: &gb_io::seq::Feature,
    chrom: &str,
    seq_len: usize,
    circular: bool,
) -> Vec<Region> {
    let feature_type = feature.kind.to_string();

    // Extract qualifiers
    let locus_tag = get_qualifier(feature, "locus_tag");
    let gene_name = get_qualifier(feature, "gene");
    let product = get_qualifier(feature, "product");

    // Determine display name: prefer locus_tag, then gene
    let name = locus_tag
        .as_deref()
        .or(gene_name.as_deref())
        .unwrap_or("")
        .to_string();

    // Parse location → (start, end, strand, parts, wraps_origin)
    let parsed = parse_location(&feature.location, seq_len, circular);

    parsed
        .into_iter()
        .map(|loc| {
            let mut region = Region::new(chrom, loc.start, loc.end)
                .with_strand(loc.strand)
                .with_name(name.clone())
                .with_tag("feature_type", feature_type.clone());

            if let Some(ref lt) = locus_tag {
                region.tags.insert("locus_tag".into(), TagValue::Str(lt.clone()));
            }
            if let Some(ref gn) = gene_name {
                region.tags.insert("gene".into(), TagValue::Str(gn.clone()));
            }
            if let Some(ref pr) = product {
                region.tags.insert("product".into(), TagValue::Str(pr.clone()));
            }
            if loc.wraps_origin {
                region.tags.insert("wraps_origin".into(), TagValue::Bool(true));
            }

            region
        })
        .collect()
}

/// Get the first value for a qualifier key.
fn get_qualifier(feature: &gb_io::seq::Feature, key: &str) -> Option<String> {
    for (qk, qv) in &feature.qualifiers {
        if qk.as_ref() == key {
            return qv.as_ref().map(|v| v.to_string());
        }
    }
    None
}

/// Parsed location result.
struct ParsedLocation {
    start: i64,
    end: i64,
    strand: Strand,
    wraps_origin: bool,
}

/// Parse a gb-io Location into (start, end, strand) tuples.
///
/// gb-io uses 0-based exclusive coordinates, matching our convention.
fn parse_location(
    loc: &gb_io::seq::Location,
    seq_len: usize,
    circular: bool,
) -> Vec<ParsedLocation> {
    match loc {
        gb_io::seq::Location::Range((start, _before), (end, _after)) => {
            vec![ParsedLocation {
                start: *start,
                end: *end,
                strand: Strand::Forward,
                wraps_origin: false,
            }]
        }

        gb_io::seq::Location::Between(a, b) => {
            vec![ParsedLocation {
                start: *a,
                end: *b,
                strand: Strand::Forward,
                wraps_origin: false,
            }]
        }

        gb_io::seq::Location::Complement(inner) => {
            let mut results = parse_location(inner, seq_len, circular);
            for r in &mut results {
                r.strand = Strand::Reverse;
            }
            results
        }

        gb_io::seq::Location::Join(parts) | gb_io::seq::Location::Order(parts) => {
            parse_compound_location(parts, seq_len, circular)
        }

        _ => vec![],
    }
}

/// Parse a compound (join/order) location.
///
/// For origin-wrapping genes on circular genomes:
///   join(end_part..seq_len, 0..start_part)
/// We detect this pattern and create virtual coordinates.
fn parse_compound_location(
    parts: &[gb_io::seq::Location],
    seq_len: usize,
    circular: bool,
) -> Vec<ParsedLocation> {
    // First, collect all simple ranges from parts
    let mut ranges: Vec<(i64, i64, Strand)> = Vec::new();
    let mut is_complement = false;

    for part in parts {
        match part {
            gb_io::seq::Location::Range((s, _), (e, _)) => {
                ranges.push((*s, *e, Strand::Forward));
            }
            gb_io::seq::Location::Complement(inner) => {
                is_complement = true;
                if let gb_io::seq::Location::Range((s, _), (e, _)) = inner.as_ref() {
                    ranges.push((*s, *e, Strand::Reverse));
                }
            }
            _ => {}
        }
    }

    if ranges.is_empty() {
        return vec![];
    }

    let strand = if is_complement {
        Strand::Reverse
    } else {
        Strand::Forward
    };

    // Detect origin-wrapping: 2 parts, one ends at seq_len, one starts at 0
    if circular && ranges.len() == 2 {
        let (s0, e0, _) = ranges[0];
        let (s1, e1, _) = ranges[1];
        let sl = seq_len as i64;

        // Pattern: join(X..seq_len, 0..Y) — wraps forward
        if e0 == sl && s1 == 0 {
            let real_start = s0;
            let span = (e0 - s0) + (e1 - s1);
            return vec![ParsedLocation {
                start: real_start,
                end: real_start + span,
                strand,
                wraps_origin: true,
            }];
        }

        // Pattern: join(0..Y, X..seq_len) — reversed wrap order
        if s0 == 0 && e1 == sl {
            let real_start = s1;
            let span = (e1 - s1) + (e0 - s0);
            return vec![ParsedLocation {
                start: real_start,
                end: real_start + span,
                strand,
                wraps_origin: true,
            }];
        }
    }

    // Non-wrapping compound: return full span
    let min_start = ranges.iter().map(|(s, _, _)| *s).min().unwrap();
    let max_end = ranges.iter().map(|(_, e, _)| *e).max().unwrap();

    vec![ParsedLocation {
        start: min_start,
        end: max_end,
        strand,
        wraps_origin: false,
    }]
}

// ─── FASTA loader ────────────────────────────────────────────────────────────

/// Load a FASTA file into a `Genome` (sequences only, no features).
pub fn load_fasta(path: &Path, circular_chroms: Option<&[&str]>) -> Result<Genome, String> {
    use bio::io::fasta;

    if !path.exists() {
        return Err(format!("FASTA file not found: {}", path.display()));
    }

    let name = path
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();

    let mut genome = Genome::new(name);

    let reader = fasta::Reader::from_file(path)
        .map_err(|e| format!("Failed to open FASTA {}: {}", path.display(), e))?;

    for record in reader.records() {
        let record = record
            .map_err(|e| format!("Failed to parse FASTA record: {}", e))?;
        let chrom = record.id().to_string();

        let is_circular = circular_chroms.map_or(false, |circs| {
            circs.iter().any(|&c| c == chrom)
        });
        let topology = if is_circular {
            Topology::Circular
        } else {
            Topology::Linear
        };

        genome.push_sequence(chrom, record.seq(), topology);
    }

    Ok(genome)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../test/fixtures")
            .join(name)
    }

    #[test]
    fn test_load_fasta_saccer3() {
        let path = fixture("sacCer3.fa");
        if !path.exists() {
            return;
        }
        let genome = load_fasta(&path, None).unwrap();
        assert!(!genome.chromosomes.is_empty());
        assert!(genome.features.is_empty());
        assert_eq!(genome.chromosomes.len(), 17);
    }

    #[test]
    fn test_load_fasta_phix174() {
        let path = fixture("phiX174.fa");
        if !path.exists() {
            return;
        }
        let genome = load_fasta(&path, Some(&["NC_001422.1"])).unwrap();
        assert_eq!(genome.chromosomes.len(), 1);
        assert!(genome.is_circular(&genome.chromosomes[0].name));
    }

    #[test]
    fn test_load_genbank_synthetic() {
        let path = fixture("synthetic_genome.gb");
        if !path.exists() {
            return;
        }
        let genome = load_genbank(&path).unwrap();

        // Two records: synth_chrom (circular) + synth_plasmid (linear)
        assert_eq!(genome.chromosomes.len(), 2);
        assert_eq!(genome.chromosomes[0].name, "synth_chrom");
        assert_eq!(genome.chromosomes[1].name, "synth_plasmid");
        assert_eq!(genome.chromosomes[0].len, 100);
        assert_eq!(genome.chromosomes[1].len, 50);

        // Topology
        assert!(genome.is_circular("synth_chrom"));
        assert!(!genome.is_circular("synth_plasmid"));

        // Features: source excluded, so we should have:
        // synth_chrom: geneA (gene), geneA (CDS), geneB (gene), geneC (gene)
        // synth_plasmid: geneD (gene)
        let genes = genome.genes();
        // gene features only (not CDS)
        assert_eq!(genes.len(), 4);

        // geneA: gene 10..30 → 0-based [9, 30)
        let gene_a = genome
            .features
            .iter()
            .find(|r| r.name == "SCH_0001" && r.tags.get("feature_type").and_then(|v| v.as_str()) == Some("gene"))
            .unwrap();
        assert_eq!(gene_a.start, 9);
        assert_eq!(gene_a.end, 30);
        assert_eq!(gene_a.strand, Strand::Forward);

        // geneB: complement(50..70) → 0-based [49, 70), reverse
        let gene_b = genome
            .features
            .iter()
            .find(|r| r.name == "SCH_0002")
            .unwrap();
        assert_eq!(gene_b.start, 49);
        assert_eq!(gene_b.end, 70);
        assert_eq!(gene_b.strand, Strand::Reverse);

        // geneC: join(90..100,1..15) on circular → wraps origin
        // 90..100 = [89, 100), 1..15 = [0, 15)
        // real_start = 89, span = 11 + 15 = 26, end = 115
        let gene_c = genome
            .features
            .iter()
            .find(|r| r.name == "SCH_0003")
            .unwrap();
        assert_eq!(gene_c.start, 89);
        assert_eq!(gene_c.end, 115);
        assert_eq!(gene_c.strand, Strand::Forward);
        assert_eq!(
            gene_c.tags.get("wraps_origin").and_then(|v| v.as_bool()),
            Some(true)
        );

        // geneD on plasmid: 10..40 → [9, 40)
        let gene_d = genome
            .features
            .iter()
            .find(|r| r.name == "SCP_0001")
            .unwrap();
        assert_eq!(gene_d.start, 9);
        assert_eq!(gene_d.end, 40);
        assert_eq!(gene_d.strand, Strand::Forward);
    }

    #[test]
    fn test_load_genbank_gz() {
        let path = fixture("synthetic_genome.gb.gz");
        if !path.exists() {
            return;
        }
        let genome = load_genbank(&path).unwrap();
        assert_eq!(genome.chromosomes.len(), 2);
        assert_eq!(genome.features.len(), 5); // 4 genes + 1 CDS
    }
}
