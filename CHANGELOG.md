# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`engine/affine.rs` — Tensor pivot + splice-aware affine extension.**
  `pivot_reads` transposes K row-major reads into column-major `[i32; 8]`
  SIMD lanes, paying the O(K×L) gather cost once before the inner loop.
  `extend_batch` advances up to 8 (read, anchor) pairs simultaneously
  through the affine splice automaton; all state transitions (substitution,
  gap open/extend, GT-AG intron detection) are zero-branch `blendv` /
  `max_epi32` operations. Scalar fallback for non-AVX2 targets.
  Exports: `pivot_reads`, `extend_batch`, `SpliceParams`.

- **`io/fastq.rs` — Streaming FASTQ parser enforcing O(C) invariant.**
  `ChunkedFastq<R: BufRead>` is an iterator yielding `Vec<FastqRecord>`
  of bounded size `chunk_size`. A single line buffer is reused across all
  records; sequences are ASCII-uppercased on read. QNAME is truncated at
  the first whitespace (SAM §1.4). Handles gzipped or plain FASTQ via any
  `BufRead` source. Exports: `ChunkedFastq`, `FastqRecord`.

- **`io/sam.rs` — Streaming SAM output sink.**
  `SamSink` implements `RegionSink` for SAM v1.6 text output. The SAM
  header (`@HD`, `@SQ` per chromosome, `@PG`) is written at construction.
  Each `consume()` call serialises one alignment line with optional
  auxiliary tags `NM:i`, `AS:i`, `NH:i`. Unmapped reads receive FLAG 4.
  `score_to_mapq` converts BWT SCORE_LUT scores to SAM MAPQ [0, 60].
  Exports: `SamSink`, `score_to_mapq`.

- **`pipeline/align.rs` — General-purpose short-read aligner.**
  `align_fastq` is the O(C) orchestration loop: `ChunkedFastq` → BWT
  seeded search → SIMD affine extension → `RegionSink`. Peak RAM is
  bounded by `chunk_size × read_len × 3`, independent of total read count.
  Multi-mapper detection sets MAPQ = 0 per SAM convention. Exports:
  `align_fastq`, `AlignConfig`, `AlignStats`, `AlignError`.

- **README: engine recharacterised as general-purpose.** Tagline, intro,
  architecture tree, and new "Rust API: General-purpose alignment" section
  document `align_fastq`, the O(C) invariant, and the tensor pivot.

### Changed

- **Single-pass master buffer architecture** — `Genome` struct replaced
  `sequences: Vec<(String, Vec<u8>)>` with a single contiguous `text: Vec<u8>`
  master buffer plus `chromosomes: Vec<ChromMeta>` metadata. Genome text now
  exists exactly once in RAM (owned by `FmIndexSearcher` after move).
  Eliminates the GenBank → temp FASTA → re-parse → FM-Index roundtrip that
  kept the genome in memory 3×.
- **`FmIndexSearcher::from_text()`** — new primary constructor takes ownership
  of a pre-built master buffer. Builds SA, BWT, BlockRank directly with zero
  file I/O, zero uppercasing, zero concatenation. `from_fasta()` kept for
  standalone FASTA usage and PyO3 bindings; refactored to delegate to
  `from_text()` internally.
- **`Genome::push_sequence()`** — builder method appends a chromosome's
  sequence directly into the master buffer (uppercase + sentinel) in a single
  pass. Used by both `load_genbank()` and `load_fasta()`.
- **Server upload path simplified** — `POST /api/files/upload` no longer
  generates an intermediate FASTA file for GenBank uploads. Parses GenBank
  once, moves the master buffer to `FmIndexSearcher::from_text()`. Removed
  `fasta_path` parameter from `GenomeStore::upload()` and `UploadRequest`.
- **Reverted accession version stripping** — the `.9` suffix hack in
  `from_fasta()` and `load_fasta()` is removed. Chrom names now pass through
  unchanged from the source file, matching GenBank accessions directly via the
  shared `ChromMeta` objects.

## [0.2.0] - 2026-03-01

### Added

- **Rust-native PAM validation in `search_batch`** — new optional parameters
  `pam`, `pam_direction`, and `topologies` enable in-Rust IUPAC bitmask PAM
  filtering. Eliminates Python-side PAM regex validation entirely.
  28x scoring speedup (111 s → 3.95 s for 893K spacers on SacCer3).
- **`score_off_targets_fast()`** — streaming drop-in replacement for
  `score_off_targets_native()` and `score_off_targets()`. Macro-chunked
  (100K regions via `islice`), deduplicates spacers per chunk, issues one
  `search_batch` call per chunk with in-Rust PAM validation. O(chunk)
  memory, Axis 8 (Generator Purity) compliant. Accepts `**kwargs` for
  `ScorerFn` signature parity with existing SeqChain call sites.
- **`scan_guides()`** — Rust-native PAM scanning + guide enrichment via
  `FmIndex.scan_guides()`. Replaces `regex_map() + interpret_guides()`.
  All coordinate math, sequence assembly, and SHA-256 guide ID hashing
  done in Rust; Python only constructs `Region` objects.
- **`GuideBuffer` PyO3 class** — SoA container returned by `scan_guides()`,
  exposes `__len__` + `__getitem__` for zero-copy lazy iteration from Python.
- **`geometry.rs`** — pure 1D coordinate algebra (`normalize`,
  `interval_envelope`, `is_low_side`, `complement_base`, `fetch_sequence`).
  All `#[inline(always)]`, zero dependencies, branchless modular arithmetic.
- **`chemistry.rs`** — molecular rules leaf module (`PamDirection`,
  `CompiledPam`, `BASE_MASK`, `generate_guide_id`). No genome coordinates,
  no search structures, no Python.
- **`python/needletail_guides.py`** — thin Python bridge providing
  `scan_guides()` and `score_off_targets_fast()` as SeqChain-compatible
  generators.
- `sha2` dependency for deterministic guide ID hashing.

### Changed

- **5-pillar codebase reorganization** aligned with the SeqChain import
  hierarchy:
  - `src/geometry.rs` — pure 1D math (leaf, no deps)
  - `src/chemistry.rs` — molecular rules (leaf, no deps)
  - `src/engine/` — `fm_index.rs`, `kmer_index.rs`, `simd_search.rs`
  - `src/io/` — `persist.rs`, `parquet_hits.rs`, `parquet_regions.rs`, `json.rs`
  - `src/operations/` — `pam_scanner.rs`
  - `src/lib.rs` — orchestration & FFI only
- `search_batch` signature extended with `max_width`, `pam`,
  `pam_direction`, `topologies` parameters (all optional, backward
  compatible).
- `bench_annotate_saccer3.py` updated to use `score_off_targets_fast`
  and `--max-width` parameter.

### Removed

- `src/primitives.rs` — split into `geometry.rs` (coordinates) and
  `chemistry.rs` (molecular rules).
- Flat `src/` layout — replaced by pillar-based module hierarchy.

## [0.1.0] - 2026-02-28

### Added

- Initial FM-Index implementation with BlockRank, AVX2 vertical automaton.
- Width-capped BWT search with tunable mismatch-branch pruning.
- Software prefetch pipeline for depth-2 pipelining.
- K=10 and K=14 seed tables for two-tier acceleration.
- rkyv zero-copy persistence with `MAP_HUGETLB` support.
- PyO3 bindings: `FmIndex.build()`, `.load()`, `.search_batch()`.
- Parquet output via Arrow zero-copy pipeline.
- Verified bit-identical results against Bowtie 1 at all mismatch levels.
