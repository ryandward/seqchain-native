//! Parquet-native Region sink — columnar output for the design pipeline.
//!
//! Structurally identical to `FileSink`: same `create()` / `consume()` /
//! `finish()` / `count()` / `path()` contract, O(row-group) memory budget.
//!
//! Schema is discovered dynamically from the first row group: the fixed
//! Region fields (chrom, start, end, strand, score) are always present,
//! and every tag seen across the buffered regions becomes a nullable column.
//! Columns are ordered by `KNOWN_ORDER` first, then alphabetically for any
//! tags the pipeline didn't exist when this file was written.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::*;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::models::region::{Region, TagValue};
use super::RegionSink;

const ROW_GROUP_SIZE: usize = 100_000;

// Tags that are internal pipeline state — never exported.
const HIDDEN_TAGS: &[&str] = &["landmark", "gene_strand"];

// Preferred column order for known tags.  Tags listed here appear in this
// sequence after the five fixed fields; any tags NOT listed here follow in
// alphabetical order.  This keeps output stable across pipeline versions
// while remaining open to new tags without code changes.
const KNOWN_ORDER: &[&str] = &[
    "guide_id", "guide_seq", "spacer", "pam_seq",
    "off_targets", "total_hits",
    "feature_type", "feature_name", "feature_strand",
    "feature_start", "feature_end",
    "signed_distance", "relative_pos", "offset", "overlap",
];

// ─── Arrow type for a TagValue ───────────────────────────────────────────────

/// The three Arrow leaf types we map TagValue variants to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColType {
    Str,
    Int,
    Float,
}

impl ColType {
    fn from_tag(v: &TagValue) -> Self {
        match v {
            TagValue::Str(_)   => ColType::Str,
            TagValue::Int(_)   => ColType::Int,
            TagValue::Float(_) => ColType::Float,
            TagValue::Bool(_)  => ColType::Str, // booleans → "true"/"false"
        }
    }

    fn arrow_type(self) -> DataType {
        match self {
            ColType::Str   => DataType::Utf8,
            ColType::Int   => DataType::Int64,
            ColType::Float => DataType::Float64,
        }
    }
}

// ─── Dynamic column builder ──────────────────────────────────────────────────

/// A column builder that wraps a single Arrow builder variant.
enum DynBuilder {
    Str(StringBuilder),
    Int(Int64Builder),
    Float(Float64Builder),
}

impl DynBuilder {
    fn new(ct: ColType) -> Self {
        match ct {
            ColType::Str   => DynBuilder::Str(StringBuilder::new()),
            ColType::Int   => DynBuilder::Int(Int64Builder::new()),
            ColType::Float => DynBuilder::Float(Float64Builder::new()),
        }
    }

    fn append_tag(&mut self, v: Option<&TagValue>) {
        match self {
            DynBuilder::Str(b) => match v {
                Some(TagValue::Str(s))  => b.append_value(s),
                Some(TagValue::Bool(v)) => b.append_value(if *v { "true" } else { "false" }),
                _ => b.append_null(),
            },
            DynBuilder::Int(b) => match v.and_then(|v| v.as_i64()) {
                Some(i) => b.append_value(i),
                None    => b.append_null(),
            },
            DynBuilder::Float(b) => match v.and_then(|v| v.as_f64()).filter(|f| f.is_finite()) {
                Some(f) => b.append_value(f),
                None    => b.append_null(),
            },
        }
    }

    fn finish(&mut self) -> Arc<dyn Array> {
        match self {
            DynBuilder::Str(b)   => Arc::new(b.finish()),
            DynBuilder::Int(b)   => Arc::new(b.finish()),
            DynBuilder::Float(b) => Arc::new(b.finish()),
        }
    }
}

// ─── ParquetFileSink ─────────────────────────────────────────────────────────

/// Disk-backed sink that writes Parquet with Snappy compression.
///
/// Regions buffer in memory up to `ROW_GROUP_SIZE`.  The first flush
/// discovers the schema (fixed fields + all tags seen so far), creates the
/// `ArrowWriter`, and writes the first row group.  Subsequent flushes
/// populate the same schema — tags absent from a later region are null.
pub struct ParquetFileSink {
    path: PathBuf,
    count: usize,
    buffer: Vec<Region>,
    /// `None` until the first flush, when the schema is discovered.
    writer: Option<ArrowWriter<File>>,
    schema: Option<Arc<Schema>>,
    /// Ordered tag columns: (tag_name, ColType).  Set on first flush.
    tag_cols: Vec<(String, ColType)>,
}

impl ParquetFileSink {
    /// Open a new Parquet file for writing.
    ///
    /// The file is created immediately (so path errors fail fast), but the
    /// Arrow schema and writer are deferred until the first row-group flush.
    pub fn create(path: &Path) -> std::io::Result<Self> {
        // Validate the path is writable by creating (and keeping) the file.
        // The File handle is dropped here; we re-create it on first flush
        // when we know the schema.
        File::create(path)?;

        Ok(ParquetFileSink {
            path: path.to_path_buf(),
            count: 0,
            buffer: Vec::with_capacity(ROW_GROUP_SIZE),
            writer: None,
            schema: None,
            tag_cols: Vec::new(),
        })
    }

    // ── Schema discovery ─────────────────────────────────────────────────

    /// Scan buffered regions and build the Arrow schema + tag column list.
    fn discover_schema(&mut self) {
        // Collect every tag name → first-seen ColType, excluding hidden.
        let mut seen: HashMap<String, ColType> = HashMap::new();
        for r in &self.buffer {
            for (k, v) in &r.tags {
                if HIDDEN_TAGS.contains(&k.as_str()) {
                    continue;
                }
                seen.entry(k.clone()).or_insert_with(|| ColType::from_tag(v));
            }
        }

        // Order: KNOWN_ORDER tags first (if present), then remaining sorted.
        let mut ordered: Vec<(String, ColType)> = Vec::with_capacity(seen.len());
        for &name in KNOWN_ORDER {
            if let Some(ct) = seen.remove(name) {
                ordered.push((name.to_string(), ct));
            }
        }
        let mut extras: Vec<(String, ColType)> = seen.into_iter().collect();
        extras.sort_by(|a, b| a.0.cmp(&b.0));
        ordered.extend(extras);

        // Build Arrow schema: 5 fixed fields + N tag columns.
        let mut fields = vec![
            Field::new("chrom",  DataType::Utf8,    false),
            Field::new("start",  DataType::Int64,   false),
            Field::new("end",    DataType::Int64,   false),
            Field::new("strand", DataType::Utf8,    false),
            Field::new("score",  DataType::Float64, true),
        ];
        for (name, ct) in &ordered {
            fields.push(Field::new(name, ct.arrow_type(), true));
        }

        self.schema = Some(Arc::new(Schema::new(fields)));
        self.tag_cols = ordered;
    }

    // ── Flush ────────────────────────────────────────────────────────────

    fn flush(&mut self) -> std::io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        // First flush: discover schema and create writer.
        if self.writer.is_none() {
            self.discover_schema();
            let file = File::create(&self.path)?;
            let props = WriterProperties::builder()
                .set_compression(Compression::SNAPPY)
                .build();
            let schema = self.schema.clone().unwrap();
            let writer = ArrowWriter::try_new(file, schema, Some(props))
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            self.writer = Some(writer);
        }

        let n = self.buffer.len();
        let schema = self.schema.clone().unwrap();

        // Fixed-field builders
        let mut chrom  = StringBuilder::with_capacity(n, n * 6);
        let mut start  = Int64Builder::with_capacity(n);
        let mut end    = Int64Builder::with_capacity(n);
        let mut strand = StringBuilder::with_capacity(n, n);
        let mut score  = Float64Builder::with_capacity(n);

        // Tag builders — one per discovered column
        let mut tag_builders: Vec<DynBuilder> = self
            .tag_cols
            .iter()
            .map(|(_, ct)| DynBuilder::new(*ct))
            .collect();

        for region in &self.buffer {
            chrom.append_value(&region.chrom);
            start.append_value(region.start);
            end.append_value(region.end);
            strand.append_value(region.strand.as_str());
            if region.score.is_finite() {
                score.append_value(region.score);
            } else {
                score.append_null();
            }

            for (i, (name, _)) in self.tag_cols.iter().enumerate() {
                tag_builders[i].append_tag(region.tags.get(name));
            }
        }

        // Assemble arrays
        let mut arrays: Vec<Arc<dyn Array>> = Vec::with_capacity(5 + tag_builders.len());
        arrays.push(Arc::new(chrom.finish()));
        arrays.push(Arc::new(start.finish()));
        arrays.push(Arc::new(end.finish()));
        arrays.push(Arc::new(strand.finish()));
        arrays.push(Arc::new(score.finish()));
        for b in &mut tag_builders {
            arrays.push(b.finish());
        }

        let batch = RecordBatch::try_new(schema, arrays)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        self.writer
            .as_mut()
            .unwrap()
            .write(&batch)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        self.buffer.clear();
        Ok(())
    }

    /// Flush remaining rows and close the Parquet writer.
    pub fn finish(mut self) -> std::io::Result<PathBuf> {
        self.flush()?;
        if let Some(writer) = self.writer {
            writer
                .close()
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        }
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

// ─── RegionSink ──────────────────────────────────────────────────────────────

impl RegionSink for ParquetFileSink {
    fn consume(&mut self, region: Region) -> std::io::Result<()> {
        self.count += 1;
        self.buffer.push(region);
        if self.buffer.len() >= ROW_GROUP_SIZE {
            self.flush()?;
        }
        Ok(())
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::region::Strand;

    #[test]
    fn test_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("guides.parquet");

        let mut sink = ParquetFileSink::create(&path).unwrap();
        sink.consume(
            Region::new("chr1", 100, 120)
                .with_strand(Strand::Forward)
                .with_score(0.85)
                .with_tag("guide_id", "g1")
                .with_tag("spacer", "ATCGATCG")
                .with_tag("off_targets", 3i64)
                .with_tag("total_hits", 10i64),
        ).unwrap();
        sink.consume(
            Region::new("chr1", 200, 220)
                .with_strand(Strand::Reverse)
                .with_score(f64::NAN)
                .with_tag("guide_id", "g2")
                .with_tag("feature_type", "exon")
                .with_tag("feature_name", "BRCA1"),
        ).unwrap();
        assert_eq!(sink.count(), 2);

        let result_path = sink.finish().unwrap();
        assert_eq!(result_path, path);

        // Read back
        let file = File::open(&path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap().build().unwrap();
        let batches: Vec<RecordBatch> = reader.collect::<Result<Vec<_>, _>>().unwrap();
        let total: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total, 2);

        let batch = &batches[0];
        let schema = batch.schema();

        // Fixed fields always present
        assert_eq!(schema.field(0).name(), "chrom");
        assert_eq!(schema.field(1).name(), "start");
        assert_eq!(schema.field(2).name(), "end");
        assert_eq!(schema.field(3).name(), "strand");
        assert_eq!(schema.field(4).name(), "score");

        // Dynamic tag columns discovered — KNOWN_ORDER tags first
        let tag_names: Vec<&str> = (5..schema.fields().len())
            .map(|i| schema.field(i).name().as_str())
            .collect();
        // guide_id should come before feature_type per KNOWN_ORDER
        let gid_pos = tag_names.iter().position(|&n| n == "guide_id").unwrap();
        let ft_pos  = tag_names.iter().position(|&n| n == "feature_type").unwrap();
        assert!(gid_pos < ft_pos);

        // Verify values
        let chrom_col = batch.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(chrom_col.value(0), "chr1");

        let start_col = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(start_col.value(0), 100);
        assert_eq!(start_col.value(1), 200);

        // Score: first row finite, second row NaN → null
        let score_idx = 4;
        let score_col = batch.column(score_idx).as_any().downcast_ref::<Float64Array>().unwrap();
        assert_eq!(score_col.value(0), 0.85);
        assert!(score_col.is_null(1));
    }

    #[test]
    fn test_row_group_flushing() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("large.parquet");

        let mut sink = ParquetFileSink::create(&path).unwrap();
        for i in 0..(ROW_GROUP_SIZE + 50) {
            sink.consume(
                Region::new("chr1", i as i64, (i + 20) as i64)
                    .with_strand(Strand::Forward)
                    .with_score(0.5)
                    .with_tag("guide_id", format!("g{}", i)),
            ).unwrap();
        }
        assert_eq!(sink.count(), ROW_GROUP_SIZE + 50);

        sink.finish().unwrap();

        let file = File::open(&path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap().build().unwrap();
        let total: usize = reader.collect::<Result<Vec<_>, _>>().unwrap()
            .iter().map(|b| b.num_rows()).sum();
        assert_eq!(total, ROW_GROUP_SIZE + 50);
    }

    #[test]
    fn test_hidden_tags_excluded() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("hidden.parquet");

        let mut sink = ParquetFileSink::create(&path).unwrap();
        sink.consume(
            Region::new("chr1", 0, 20)
                .with_tag("landmark", "TSS")
                .with_tag("gene_strand", "+")
                .with_tag("guide_id", "g1"),
        ).unwrap();
        sink.finish().unwrap();

        let file = File::open(&path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap().build().unwrap();
        let batches: Vec<RecordBatch> = reader.collect::<Result<Vec<_>, _>>().unwrap();
        let schema = batches[0].schema();
        for field in schema.fields() {
            assert!(!HIDDEN_TAGS.contains(&field.name().as_str()),
                "hidden tag '{}' should not appear in schema", field.name());
        }
        // guide_id should still be present
        assert!(schema.field_with_name("guide_id").is_ok());
    }

    #[test]
    fn test_novel_tags_appear_as_columns() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("novel.parquet");

        let mut sink = ParquetFileSink::create(&path).unwrap();
        sink.consume(
            Region::new("chr1", 0, 20)
                .with_tag("guide_id", "g1")
                .with_tag("gc_content", 0.55f64)
                .with_tag("custom_label", "high"),
        ).unwrap();
        sink.consume(
            Region::new("chr1", 30, 50)
                .with_tag("guide_id", "g2")
                .with_tag("gc_content", 0.42f64),
        ).unwrap();
        sink.finish().unwrap();

        let file = File::open(&path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap().build().unwrap();
        let batches: Vec<RecordBatch> = reader.collect::<Result<Vec<_>, _>>().unwrap();
        let schema = batches[0].schema();

        // Novel tags should appear as columns
        assert!(schema.field_with_name("gc_content").is_ok());
        assert!(schema.field_with_name("custom_label").is_ok());

        // gc_content should be Float64
        let gc = schema.field_with_name("gc_content").unwrap();
        assert_eq!(gc.data_type(), &DataType::Float64);

        // custom_label should be Utf8
        let cl = schema.field_with_name("custom_label").unwrap();
        assert_eq!(cl.data_type(), &DataType::Utf8);

        // guide_id (KNOWN_ORDER) should come before gc_content/custom_label
        let gid_idx = schema.index_of("guide_id").unwrap();
        let gc_idx  = schema.index_of("gc_content").unwrap();
        let cl_idx  = schema.index_of("custom_label").unwrap();
        assert!(gid_idx < gc_idx);
        assert!(gid_idx < cl_idx);
        // Extras are alphabetical: custom_label < gc_content
        assert!(cl_idx < gc_idx);
    }

    #[test]
    fn test_empty_sink() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.parquet");

        let sink = ParquetFileSink::create(&path).unwrap();
        assert_eq!(sink.count(), 0);
        // finish() with zero regions — file exists but is effectively empty
        let result = sink.finish().unwrap();
        assert_eq!(result, path);
    }
}
