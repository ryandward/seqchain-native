//! SAM output sink: implements `RegionSink` for sequence alignment data.
//!
//! Writes a text SAM file (SAM v1.6, §1.4). Each `consume()` call emits
//! one tab-delimited alignment line. No intermediate buffering — bytes
//! flow directly from the Region into the kernel's page cache via
//! `BufWriter`. The Water leaves the mountain without pooling.
//!
//! # Region → SAM field mapping
//!
//! | SAM field | Source                              |
//! |-----------|-------------------------------------|
//! | QNAME     | `region.name`                       |
//! | FLAG      | `region.strand` + tags["is_secondary"] |
//! | RNAME     | `region.chrom`                      |
//! | POS       | `region.start + 1` (1-based)        |
//! | MAPQ      | tags["mapq"] or derived from score  |
//! | CIGAR     | tags["cigar"] or `{len}M`           |
//! | RNEXT     | `*` (single-end)                    |
//! | PNEXT     | `0`                                 |
//! | TLEN      | `0`                                 |
//! | SEQ       | tags["seq"] or `*`                  |
//! | QUAL      | tags["qual"] or `*`                 |
//!
//! Optional SAM auxiliary tags appended when present:
//! `NM:i` (edit distance), `AS:i` (alignment score), `NH:i` (hit count).
//!
//! # Import hierarchy
//! Depends on `models::region`. No engine, no biology.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::io::RegionSink;
use crate::models::region::{Region, Strand, TagValue};

// ── SAM FLAG bits ──────────────────────────────────────────────────────────
const FLAG_REVERSE:     u16 = 0x10;  // SEQ on reverse strand
const FLAG_SECONDARY:   u16 = 0x100; // not primary alignment
const FLAG_UNMAPPED:    u16 = 0x04;  // segment unmapped

// ═══════════════════════════════════════════════════════════════════════════
//  SamSink
// ═══════════════════════════════════════════════════════════════════════════

/// Terminal sink that serialises `Region` objects to a SAM text file.
///
/// The SAM header (`@HD`, `@SQ`, `@PG`) is written at construction.
/// Each subsequent `consume()` call writes one alignment line.
/// `finish()` flushes the writer and returns the output path.
pub struct SamSink {
    writer: BufWriter<File>,
    path:   std::path::PathBuf,
    count:  usize,
}

impl SamSink {
    /// Create a new SAM file at `path`.
    ///
    /// Writes the SAM header immediately:
    /// - `@HD VN:1.6 SO:unsorted` — no pre-sorting guarantee
    /// - `@SQ SN:{name} LN:{len}` — one line per chromosome
    /// - `@PG ID:needletail PN:needletail VN:0.2.0`
    ///
    /// `chroms`: `(name, length_in_bases)` pairs in index order.
    pub fn create(path: &Path, chroms: &[(String, usize)]) -> io::Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // @HD
        writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")?;

        // @SQ — one per reference sequence
        for (name, len) in chroms {
            writeln!(writer, "@SQ\tSN:{name}\tLN:{len}")?;
        }

        // @PG
        writeln!(
            writer,
            "@PG\tID:needletail\tPN:needletail\tVN:0.2.0\tCL:needletail-align"
        )?;

        Ok(Self {
            writer,
            path: path.to_path_buf(),
            count: 0,
        })
    }

    /// Number of alignment records written.
    pub fn count(&self) -> usize {
        self.count
    }

    /// Flush the writer and return the output path.
    pub fn finish(mut self) -> io::Result<std::path::PathBuf> {
        self.writer.flush()?;
        Ok(self.path)
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  RegionSink implementation
// ═══════════════════════════════════════════════════════════════════════════

impl RegionSink for SamSink {
    /// Serialise one `Region` as a SAM alignment line.
    ///
    /// Unmapped reads (score ≤ 0 or chrom == "*") receive FLAG 4
    /// and `RNAME=* POS=0 MAPQ=0 CIGAR=*`.
    fn consume(&mut self, region: Region) -> io::Result<()> {
        self.count += 1;

        let unmapped = region.chrom == "*" || region.score.map_or(true, |s| s <= 0.0);

        // ── FLAG ──────────────────────────────────────────────────────────
        let mut flag: u16 = 0;
        if unmapped {
            flag |= FLAG_UNMAPPED;
        } else {
            if region.strand == Strand::Reverse {
                flag |= FLAG_REVERSE;
            }
            if is_secondary(&region) {
                flag |= FLAG_SECONDARY;
            }
        }

        // ── RNAME / POS ───────────────────────────────────────────────────
        let (rname, pos) = if unmapped {
            ("*".to_string(), 0u64)
        } else {
            (region.chrom.clone(), (region.start + 1).max(0) as u64)
        };

        // ── MAPQ ──────────────────────────────────────────────────────────
        let mapq: u8 = if unmapped {
            0
        } else {
            tag_u8(&region, "mapq").unwrap_or_else(|| score_to_mapq(region.score.unwrap_or(0.0)))
        };

        // ── CIGAR ─────────────────────────────────────────────────────────
        let read_len = (region.end - region.start).max(0) as usize;
        let cigar = if unmapped {
            "*".to_string()
        } else {
            tag_str(&region, "cigar")
                .unwrap_or_else(|| format!("{read_len}M"))
        };

        // ── SEQ / QUAL ────────────────────────────────────────────────────
        let seq  = tag_str(&region, "seq").unwrap_or_else(|| "*".to_string());
        let qual = if seq == "*" {
            "*".to_string()
        } else {
            tag_str(&region, "qual").unwrap_or_else(|| "*".to_string())
        };

        // ── Mandatory fields ──────────────────────────────────────────────
        write!(
            self.writer,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}",
            qname = if region.name.is_empty() { "*" } else { &region.name },
            flag  = flag,
            rname = rname,
            pos   = pos,
            mapq  = mapq,
            cigar = cigar,
            seq   = seq,
            qual  = qual,
        )?;

        // ── Optional auxiliary tags ────────────────────────────────────────
        if let Some(nm) = tag_i64(&region, "nm") {
            write!(self.writer, "\tNM:i:{nm}")?;
        }
        if !unmapped {
            // AS:i — alignment score scaled to [0, 255]
            let as_score = (region.score.unwrap_or(0.0) * 100.0).round().clamp(0.0, 255.0) as u8;
            write!(self.writer, "\tAS:i:{as_score}")?;
        }
        if let Some(nh) = tag_i64(&region, "nhits") {
            write!(self.writer, "\tNH:i:{nh}")?;
        }

        writeln!(self.writer)?;
        Ok(())
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Private helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Convert a [0, 1] alignment score to SAM MAPQ [0, 60].
///
/// Uses the BWT SCORE_LUT convention: score = 1 / (1 + mm).
///   mm = 0  → score = 1.0     → MAPQ = 60
///   mm = 1  → score = 0.5     → MAPQ = 30
///   mm = 2  → score = 0.333   → MAPQ = 20
///   mm = 3  → score = 0.25    → MAPQ = 15
#[inline]
pub fn score_to_mapq(score: f64) -> u8 {
    (score * 60.0).round().clamp(0.0, 60.0) as u8
}

#[inline]
fn is_secondary(region: &Region) -> bool {
    matches!(region.tags.get("is_secondary"), Some(TagValue::Bool(true)))
}

#[inline]
fn tag_str(region: &Region, key: &str) -> Option<String> {
    match region.tags.get(key)? {
        TagValue::Str(s) => Some(s.clone()),
        _ => None,
    }
}

#[inline]
fn tag_i64(region: &Region, key: &str) -> Option<i64> {
    region.tags.get(key)?.as_i64()
}

#[inline]
fn tag_u8(region: &Region, key: &str) -> Option<u8> {
    let v = tag_i64(region, key)?;
    Some(v.clamp(0, 255) as u8)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Region;
    use tempfile::NamedTempFile;

    fn make_region(chrom: &str, start: i64, end: i64, strand: Strand, score: f64) -> Region {
        Region::new(chrom, start, end)
            .with_strand(strand)
            .with_score(score)
            .with_name("read1")
            .with_tag("seq",  TagValue::Str("ACGT".to_string()))
            .with_tag("qual", TagValue::Str("IIII".to_string()))
            .with_tag("nm",   TagValue::Int(0))
            .with_tag("nhits",TagValue::Int(1))
    }

    #[test]
    fn writes_header_and_record() {
        let tmp = NamedTempFile::new().unwrap();
        let chroms = vec![("chrI".to_string(), 230218usize)];
        let mut sink = SamSink::create(tmp.path(), &chroms).unwrap();

        let region = make_region("chrI", 0, 4, Strand::Forward, 1.0);
        sink.consume(region).unwrap();
        let path = sink.finish().unwrap();

        let contents = std::fs::read_to_string(path).unwrap();
        assert!(contents.contains("@HD\tVN:1.6"));
        assert!(contents.contains("@SQ\tSN:chrI\tLN:230218"));
        assert!(contents.contains("read1\t0\tchrI\t1\t60\t4M"));
    }

    #[test]
    fn reverse_strand_sets_flag_16() {
        let tmp = NamedTempFile::new().unwrap();
        let mut sink = SamSink::create(tmp.path(), &[]).unwrap();
        let region = make_region("chrI", 0, 4, Strand::Reverse, 1.0);
        sink.consume(region).unwrap();
        let path = sink.finish().unwrap();
        let contents = std::fs::read_to_string(path).unwrap();
        // FLAG = 16 for reverse strand.
        assert!(contents.contains("read1\t16\t"));
    }

    #[test]
    fn unmapped_record() {
        let tmp = NamedTempFile::new().unwrap();
        let mut sink = SamSink::create(tmp.path(), &[]).unwrap();
        let region = Region::new("*", 0, 0)
            .with_name("unmap1")
            .with_tag("seq",  TagValue::Str("ACGT".to_string()))
            .with_tag("qual", TagValue::Str("IIII".to_string()));
        sink.consume(region).unwrap();
        let path = sink.finish().unwrap();
        let contents = std::fs::read_to_string(path).unwrap();
        // FLAG = 4 (unmapped), RNAME = *, POS = 0.
        assert!(contents.contains("unmap1\t4\t*\t0\t0\t*"));
    }

    #[test]
    fn score_to_mapq_values() {
        assert_eq!(score_to_mapq(1.0), 60);
        assert_eq!(score_to_mapq(0.5), 30);
        assert_eq!(score_to_mapq(0.0), 0);
    }
}
