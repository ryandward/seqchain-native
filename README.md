<p align="center">
  <img src="logo.png" alt="Needletail" width="500">
</p>

<h1 align="center">Needletail</h1>

<p align="center"><strong>A SIMD-accelerated FM-Index aligner for CRISPR guide design.</strong></p>

Needletail is a Rust + PyO3 short-read aligner purpose-built for the
off-target scoring bottleneck in CRISPR guide design pipelines.  It replaces
Bowtie 1's subprocess-and-disk workflow with an in-process, memory-mapped
FM-Index searched by an AVX2 vertical automaton — delivering 3–28x faster
mismatch search while producing bit-identical results.

Needletail ships as a native Python extension (`needletail`) and
integrates directly into the
[SeqChain](https://github.com/ryandward/SeqChain) generator pipeline as a
drop-in `ScorerFn` operation.

```python
from needletail import FmIndex

idx = FmIndex.load("sacCer3.seqchain")      # 51 µs  (mmap, zero-copy)
idx.warm_seeds("sacCer3.seqchain")           # 5.6 s  (K=10 + K=14 tables)

# Search with in-Rust PAM validation — no Python-side regex needed
qi, pos, strand, scores = idx.search_batch(
    spacers, mismatches=2, max_width=8,
    pam="NGG", pam_direction="downstream",   # IUPAC bitmask filter
)
```

---

## Why not just use Bowtie?

Bowtie 1 is a masterpiece of BWT engineering.  It has been the standard
short-read aligner for over fifteen years.  But its architecture forces
three compromises that hurt CRISPR pipelines:

| Constraint | Bowtie 1 | Needletail |
|---|---|---|
| **Process model** | Subprocess per invocation | In-process FFI, GIL released |
| **Index loading** | Read + deserialize `.ebwt` (10 s) | `mmap()` + pointer cast (51 µs) |
| **Mismatch search** | Backtracking DFS, scalar | AVX2 vertical automaton, 256-wide |
| **PAM validation** | Shell-side regex (Python) | In-Rust IUPAC bitmask engine |
| **Output format** | Write SAM/bowtie to disk, re-parse | Return NumPy-compatible arrays |
| **Pipeline integration** | Shell pipes, temp files | Python generator (`ScorerFn`) |

The result is identical hit counts — verified at every mismatch level against
Bowtie 1 on the full *S. cerevisiae* genome — with no subprocess overhead, no
temp files, and no SAM parsing.

---

## Benchmarks

All benchmarks on the *Saccharomyces cerevisiae* S288C reference genome
(SacCer3, 12.16 Mbp, 17 chromosomes).  Machine: AMD Zen 4, 16 cores.

### Off-target scoring: 893,268 unique spacers

NGG PAM, 20 bp spacer, mm=2, `max_width=8`, PAM-validated in Rust.

| Stage | Needletail | Bowtie 1 pipeline |
|-------|------------|-------------------|
| search_batch (BWT + PAM filter) | **3.95 s** | 111 s |
| Scoring speedup | | **28x** |

967,398 PAM-valid hits.  70,425 promoter-targeting guides, 5,443 genes,
93.8% perfectly specific (score = 0.0).

### End-to-end CRISPRi pipeline

Full pipeline: scan + score + annotate + filter + TSS distance + JSON export.

```
70,425 promoter-targeting guides
 5,443 unique genes targeted
66,046 perfectly specific (93.8%)

Genome load  : 1.63s
mmap load    : 169602 µs
warm_seeds   : 6.29s
Pipeline drain: 20.93s
Wall-clock    : 29.05s
```

### Chromosome 1 (230 kb): scaling across mismatch levels

48,458 unique spacers, TTN PAM, report all hits.

| Mismatch | Hits | Needletail | Bowtie 1 | Speedup |
|----------|------|------------|----------|---------|
| mm=0 | 48,458 | 0.044 s | 0.209 s | **4.7x** |
| mm=1 | 51,969 | 0.212 s | 0.530 s | **2.5x** |
| mm=2 | 56,883 | 0.118 s | 0.707 s | **6.0x** |
| mm=3 | 69,946 | 0.064 s | 1.555 s | **24.2x** |

The advantage grows with mismatch budget because the vertical automaton
processes all 4 mismatch tiers in a single pass, while Bowtie's backtracking
DFS explores an exponentially larger tree.

---

## Architecture

Needletail follows the **SeqChain 5-pillar import hierarchy**: strict
layering from pure math to FFI, where each layer depends only on layers
below it.

```
src/
├── geometry.rs         Pure 1D coordinate math (leaf, no deps)
│                       normalize · interval_envelope · is_low_side
│                       complement_base · fetch_sequence
│
├── chemistry.rs        Molecular rules (leaf, no deps)
│                       PamDirection · CompiledPam · BASE_MASK
│                       generate_guide_id
│
├── engine/             BWT search structures
│   ├── fm_index.rs     BlockRank, rank/occ queries
│   ├── kmer_index.rs   K=10 / K=14 seed tables, PosTable
│   └── simd_search.rs  AVX2 vertical automaton, width-first BFS
│
├── io/                 Disk boundary
│   ├── persist.rs      rkyv zero-copy serialization, mmap
│   └── parquet_sink.rs Arrow zero-copy Parquet export
│
├── operations/         Hot loops (composition layer)
│   └── pam_scanner.rs  PAM site scanning, guide enrichment,
│                       PAM validation, off-target filtering
│
└── lib.rs              Orchestration & FFI (only file that imports pyo3)
    FmIndex · GuideBuffer · search_batch · scan_guides
```

### The FM-Index: BlockRank

The suffix array and BWT are built by the
[`bio`](https://docs.rs/bio) crate.  Rank queries use a cache-line-aligned
BlockRank structure: every 512 symbols, a superblock stores the cumulative
rank; every 128 symbols within, a miniblock stores the local delta.  The
inner loop resolves ranks with `popcnt` on 64-bit words — no table lookups,
no branch mispredictions.

### Persistence: Zero-Copy mmap

The index serializes via [`rkyv`](https://rkyv.org) (zero-copy
deserialization) to a `.seqchain` file.  Loading is a single `mmap()` +
pointer cast — **51 µs** for the 290 MB SacCer3 index.  The kernel
demand-pages only the blocks actually touched during search.

### Search: The Vertical SIMD Automaton

The mismatch search is a width-first BFS over the FM-Index, processing
frontier nodes in groups.  For groups of 256+ items, an **AVX2 vertical
automaton** replaces the per-item scalar loop:

**Match mask precomputation** — For each BWT symbol `c ∈ {A, C, G, T}`,
a 256-bit register `M_c` encodes which of the 256 work items have
`bwt_char == c` (bit set = match):

```
M_A[lane] = bitmask of items in lane whose BWT character is A
```

**Budget-remaining recurrence** — Four 256-bit registers `R_0..R_3`
track the mismatch budget state.  Bit `i` of `R_j` is set iff work item `i`
has budget ≥ `j`.  The single-step recurrence:

```
R'_3 = R_3 ∧ M_c                        // budget=3 survives only on match
R'_2 = (R_2 ∧ M_c) ∨ (R_3 ∧ ¬M_c)      // match keeps 2, mismatch demotes 3→2
R'_1 = (R_1 ∧ M_c) ∨ (R_2 ∧ ¬M_c)      // match keeps 1, mismatch demotes 2→1
R'_0 = (R_0 ∧ M_c) ∨ (R_1 ∧ ¬M_c)      // match keeps 0, mismatch demotes 1→0
```

**Branchless pruning** — Items with budget < 0 are killed:
`alive = R_0 ≠ 0` via `_mm256_testz_si256`.

**13-register map** — The inner loop fits entirely in YMM registers:
ymm0–3 (match masks), ymm4–7 (budget state), ymm8 (all-ones), ymm9 (M_c),
ymm10 (¬M_c), ymm11–12 (scratch).  Zero spills to memory.

For groups smaller than 256, an SSE2 fallback uses 128-bit registers with
the same logic.  Runtime detection (`is_x86_feature_detected!("avx2")`)
selects the path.

### IUPAC Bitmask PAM Validation

PAM validation runs in Rust after the BWT search completes. Each IUPAC
character in the PAM pattern is compiled to a 4-bit mask (A=0001, C=0010,
G=0100, T=1000, N=1111, etc.). At each hit position, the engine checks
`BASE_MASK[genome[pos]] & pam_mask[j] != 0` — O(1) per position with no
regex overhead. Both forward and reverse-complement PAMs are checked via
precomputed mask arrays.

### Seed Tables: Two-Tier Acceleration

K-mer seed tables pre-index the genome at K=10 (8 MiB, for mm ≤ 2) and K=14
(2 GiB, for mm = 3).  For each query, the engine extracts a seed k-mer,
looks up candidate SA intervals, and restricts the BFS to only those
intervals.  The `select_tier()` router dynamically picks K=10 or K=14 based
on the mismatch budget and query length.

### Python API

```python
class FmIndex:
    def __new__(fasta_path: str) -> FmIndex: ...
    def build(fasta_path: str, index_path: str) -> FmIndex: ...
    def load(index_path: str) -> FmIndex: ...
    def warm_seeds(index_path: str) -> None: ...
    def search_batch(
        queries: list[str],
        mismatches: int = 0,
        max_width: int = 2**32,
        pam: str | None = None,
        pam_direction: str = "downstream",
        topologies: list[bool] | None = None,
    ) -> tuple[list[int], list[int], list[bool], list[float]]: ...
    def search_batch_exhaustive(
        queries: list[str], mismatches: int = 0,
    ) -> tuple[list[int], list[int], list[bool], list[float]]: ...
    def scan_guides(
        pam: str, spacer_len: int,
        pam_direction: str = "downstream",
        topologies: list[bool] | None = None,
    ) -> GuideBuffer: ...
    def chrom_names() -> list[str]: ...
    def chrom_ranges() -> list[tuple[int, int]]: ...

class GuideBuffer:
    count: int
    spacer_len: int
    pam_len: int
    def __len__() -> int: ...
    def __getitem__(i: int) -> tuple: ...
```

The GIL is released during `search_batch`, `scan_guides`, `build`, and
`warm_seeds`. Rayon parallelizes across query groups internally.

---

## SeqChain Integration

Needletail plugs into the SeqChain generator architecture with two
high-level bridges in `python/needletail_guides.py`:

### `scan_guides()` — Rust-native PAM scanning

Replaces `regex_map() + interpret_guides()`. All coordinate math, sequence
assembly, and guide ID hashing done in Rust. Python only constructs `Region`
objects.

```python
from needletail_guides import scan_guides

guides = scan_guides(idx, preset, genome)  # lazy generator
```

### `score_off_targets_fast()` — Streaming Rust-native scoring

Drop-in replacement for `score_off_targets()` and `score_off_targets_native()`.
Pulls 100K regions at a time from the upstream generator via `islice`,
deduplicates spacers within each chunk, issues a single `search_batch` call
with in-Rust PAM validation, yields scored Regions one by one, then frees
the chunk. O(100K) memory — Axis 8 (Generator Purity) compliant.

Accepts `**kwargs` for `ScorerFn` signature parity: extra parameters like
`sequences=` are silently swallowed so existing call sites work unchanged.

```python
from needletail_guides import score_off_targets_fast

scored = score_off_targets_fast(
    guides, index, preset.pam,
    spacer_len=20, pam_direction="downstream",
    mismatches=2, max_width=8,
    topologies=genome.topologies,
    sequences=genome.sequences,  # accepted, ignored — PAM validated in Rust
)

# Continue the generator chain
tiles     = anchor_features(genome.genes(), config, genome.chrom_lengths)
annotated = annotate_tracks(scored, tiles)
promoters = (g for g in annotated if g.tags.get("feature_type") == "promoter")
final     = append_tss_distance(promoters)

stream_json_export("output.json", final)
```

---

## Building

```bash
# Create venv and install
python -m venv .venv
source .venv/bin/activate
pip install maturin
maturin develop --release

# Verify
python -c "from needletail import FmIndex; print('OK')"
```

### Build an index

```python
from needletail import FmIndex

# Build from FASTA and persist
idx = FmIndex.build("genome.fa", "genome.seqchain")

# Later: mmap load (51 µs)
idx = FmIndex.load("genome.seqchain")
idx.warm_seeds("genome.seqchain")
```

### Run a search

```python
spacers = ["ATCGATCGATCGATCGATCG", "TTTTAAAACCCCGGGGAAAA"]

# Basic search
qi, pos, strand, scores = idx.search_batch(spacers, mismatches=2)

# With PAM validation (no Python post-processing needed)
qi, pos, strand, scores = idx.search_batch(
    spacers, mismatches=2, max_width=8,
    pam="NGG", pam_direction="downstream",
)

# qi[i]     = index into spacers list for hit i
# pos[i]    = global genome position
# strand[i] = True (forward) / False (reverse)
# scores[i] = mismatch count
```

---

## File format

| File | Contents | Size (SacCer3) |
|------|----------|----------------|
| `*.seqchain` | rkyv-serialized FM-Index (BWT + SA + BlockRank + metadata) | 290 MiB |
| `*.10mer.idx` | K=10 seed table (4^10 × 8 bytes) | 8 MiB |
| `*.14mer.idx` | K=14 seed table (4^14 × 8 bytes) | 2 GiB |

The `.seqchain` file is architecture-specific (pointer width, endianness).
The seed tables are flat binary arrays, regenerated by `warm_seeds()` if
missing.

---

## Correctness

Every mismatch level (0–3) has been validated against Bowtie 1 on the full
SacCer3 genome:

- **Raw hit counts match exactly** (e.g., 1,790,098 at mm=2)
- **PAM-validated counts match exactly** (1,076,797 at mm=2)
- **Unique guide counts match exactly** (852,930 at mm=2)
- **Per-spacer off-target scores match exactly**

The `search_batch_exhaustive` method provides a ground-truth unseeded path
for verifying the seeded search against.

---

## Dependencies

| Crate | Purpose |
|-------|---------|
| `pyo3` | Python extension module (abi3-py39+) |
| `bio` | Suffix array + BWT construction |
| `rkyv` | Zero-copy serialization |
| `memmap2` | Memory-mapped I/O |
| `rayon` | Data-parallel search |
| `sha2` | Guide ID hashing (SHA-256) |
| `anyhow` | Error handling |
| `tempfile` | Seed table build scratch |
| `arrow` + `parquet` | Zero-copy Parquet export |

Release builds use LTO for cross-crate inlining of the SIMD inner loops.

---

## License

MIT
