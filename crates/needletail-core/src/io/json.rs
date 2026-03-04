//! JSON serialization for Region output.
//!
//! Converts Regions to JSON objects and supports streaming JSON array output.
//! NaN/Inf values are sanitized to null.
//!
//! `FileSink` implements `RegionSink` for disk-backed JSONL output — the
//! terminal materialization point for the whole-genome pipeline.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use serde_json::{json, Map, Value};

use crate::models::region::{Region, TagValue};
use super::{RegionSink, HIDDEN_TAGS};

/// Convert a Region to a JSON value.
///
/// Fixed fields (chrom, start, end, strand, score) come first, then all
/// tags sorted alphabetically.  Internal tags (HIDDEN_TAGS) are excluded.
/// The Stone has no opinion on biological column ordering.
pub fn region_to_json(r: &Region) -> Value {
    let mut obj = Map::new();

    obj.insert("chrom".into(),  json!(r.chrom));
    obj.insert("start".into(),  json!(r.start));
    obj.insert("end".into(),    json!(r.end));
    obj.insert("strand".into(), json!(r.strand.as_str()));
    obj.insert("score".into(),  sanitize_float(r.score));

    let mut tags: Vec<(&str, Value)> = r.tags
        .iter()
        .filter(|(k, _)| !HIDDEN_TAGS.contains(&k.as_str()))
        .map(|(k, v)| (k.as_str(), tag_value_to_json(v)))
        .collect();
    tags.sort_by_key(|(k, _)| *k);

    for (k, v) in tags {
        obj.insert(k.into(), v);
    }

    Value::Object(obj)
}

/// Convert a TagValue to a JSON value.
fn tag_value_to_json(v: &TagValue) -> Value {
    match v {
        TagValue::Int(i) => json!(i),
        TagValue::Float(f) => sanitize_float(*f),
        TagValue::Str(s) => json!(s),
        TagValue::Bool(b) => json!(b),
    }
}

/// Sanitize a float: NaN and Inf become null.
fn sanitize_float(f: f64) -> Value {
    if f.is_finite() {
        json!(f)
    } else {
        Value::Null
    }
}

/// Write a JSON array of Regions to a writer, streaming one object at a time.
pub fn stream_json_array<W: Write>(regions: &[Region], writer: &mut W) -> std::io::Result<()> {
    writer.write_all(b"[")?;
    for (i, r) in regions.iter().enumerate() {
        if i > 0 {
            writer.write_all(b",")?;
        }
        let val = region_to_json(r);
        serde_json::to_writer(&mut *writer, &val)?;
    }
    writer.write_all(b"]")?;
    Ok(())
}

/// Serialize Regions to a JSON string.
pub fn regions_to_json_string(regions: &[Region]) -> String {
    let mut buf = Vec::new();
    stream_json_array(regions, &mut buf).expect("writing to Vec never fails");
    String::from_utf8(buf).expect("JSON is always valid UTF-8")
}

// ─── FileSink ─────────────────────────────────────────────────────────────────

/// Disk-backed sink that writes one JSON object per line (JSONL).
///
/// Each call to `consume()` serializes a Region and appends it as a line
/// to the file.  The sweep-line drains directly here — no intermediate
/// Vec, no RAM accumulation.  The file IS the product.
///
/// The first line is `[` and the last line (written by `finish()`) is `]`,
/// making the file a valid JSON array that can be streamed to HTTP clients
/// or loaded by downstream tools.
pub struct FileSink {
    writer: BufWriter<File>,
    path: PathBuf,
    count: usize,
}

impl FileSink {
    /// Open a new JSONL file for writing.
    pub fn create(path: &Path) -> std::io::Result<Self> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(64 * 1024, file);
        // Open the JSON array
        writer.write_all(b"[\n")?;
        Ok(FileSink {
            writer,
            path: path.to_path_buf(),
            count: 0,
        })
    }

    /// Flush and close the JSON array.  Must be called after the pipeline
    /// finishes to write the closing `]`.
    pub fn finish(mut self) -> std::io::Result<PathBuf> {
        self.writer.write_all(b"\n]")?;
        self.writer.flush()?;
        Ok(self.path)
    }

    /// Number of regions written so far.
    pub fn count(&self) -> usize {
        self.count
    }

    /// Path of the output file.
    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl RegionSink for FileSink {
    fn consume(&mut self, region: Region) -> std::io::Result<()> {
        if self.count > 0 {
            self.writer.write_all(b",\n")?;
        }
        let val = region_to_json(&region);
        serde_json::to_writer(&mut self.writer, &val)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        self.count += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_region_to_json() {
        let r = Region::new("chr1", 100, 200)
            .with_strand(Strand::Forward)
            .with_name("guide1")
            .with_score(2.0)
            .with_tag("spacer", "ATCG")
            .with_tag("off_targets", 5i64);

        let j = region_to_json(&r);
        assert_eq!(j["chrom"], "chr1");
        assert_eq!(j["start"], 100);
        assert_eq!(j["end"], 200);
        assert_eq!(j["strand"], "+");
        assert_eq!(j["score"], 2.0);
        assert_eq!(j["spacer"], "ATCG");
        assert_eq!(j["off_targets"], 5);
    }

    #[test]
    fn test_nan_sanitization() {
        let r = Region::new("chr1", 0, 10).with_score(f64::NAN);
        let j = region_to_json(&r);
        assert!(j["score"].is_null());
    }

    #[test]
    fn test_stream_json_array() {
        let regions = vec![
            Region::new("chr1", 0, 10).with_name("a"),
            Region::new("chr1", 20, 30).with_name("b"),
        ];
        let json_str = regions_to_json_string(&regions);
        let parsed: Vec<Value> = serde_json::from_str(&json_str).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0]["chrom"], "chr1");
        assert_eq!(parsed[1]["chrom"], "chr1");
    }

    #[test]
    fn test_file_sink() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("guides.json");

        let mut sink = FileSink::create(&path).unwrap();
        sink.consume(Region::new("chr1", 0, 20).with_name("g1").with_score(0.0)).unwrap();
        sink.consume(Region::new("chr1", 100, 120).with_name("g2").with_score(1.0)).unwrap();
        assert_eq!(sink.count(), 2);

        let result_path = sink.finish().unwrap();
        assert_eq!(result_path, path);

        // Read back and verify it's a valid JSON array
        let contents = std::fs::read_to_string(&path).unwrap();
        let parsed: Vec<Value> = serde_json::from_str(&contents).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0]["chrom"], "chr1");
        assert_eq!(parsed[1]["chrom"], "chr1");
    }
}
