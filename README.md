# Needletail

**A SIMD-accelerated FM-Index aligner for CRISPR guide design.**

Needletail is a Rust + PyO3 short-read aligner purpose-built for the
off-target scoring bottleneck in CRISPR guide design pipelines.  It replaces
Bowtie 1's subprocess-and-disk workflow with an in-process, memory-mapped
FM-Index searched by an AVX2 vertical automaton — delivering 3–24x faster
mismatch search while producing bit-identical results.

Needletail ships as a native Python extension (`needletail`) and
integrates directly into the
[SeqChain](https://github.com/ryandward/SeqChain) generator pipeline as a
drop-in `ScorerFn` operation.

```python
from needletail import FmIndex

idx = FmIndex.load("sacCer3.seqchain")      # 51 µs  (mmap, zero-copy)
idx.warm_seeds("sacCer3.seqchain")           # 5.6 s  (K=10 + K=14 tables)

qi, pos, strand, scores = idx.search_batch(
    spacers, mismatches=2,                   # 33 µs/query, 893K queries
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
| **Output format** | Write SAM/bowtie to disk, re-parse | Return NumPy-compatible arrays |
| **Pipeline integration** | Shell pipes, temp files | Python generator (`ScorerFn`) |

The result is identical hit counts — verified at every mismatch level against
Bowtie 1 on the full *S. cerevisiae* genome — with no subprocess overhead, no
temp files, and no SAM parsing.

---

## Benchmarks

All benchmarks on the *Saccharomyces cerevisiae* S288C reference genome
(SacCer3, 12.16 Mbp, 17 chromosomes).  Machine: AMD Zen 4, 16 cores.

### Genome-wide SpCas9 guide design: 893,267 unique spacers

NGG PAM, 20 bp spacer, report all hits (`-a`), PAM-validated.

| Mismatch | Raw hits | PAM-valid | Needletail | Bowtie 1 | Speedup |
|----------|----------|-----------|------------|----------|---------|
| mm=2 | 1,790,098 | 1,076,797 | **29.8 s** | 100.8 s | **3.4x** |

Hit counts match exactly.  852,930 unique guides (score = 0) in both.

**End-to-end pipeline** (load + warm + search + PAM validate):

| | Needletail (hot) | Bowtie 1 |
|---|---|---|
| Index load | 51 µs (mmap) | 10.3 s (build) |
| Search mm=2 | 29.8 s | 100.8 s |
| PAM validate | 1.9 s | 3.1 s |
| **Total** | **31.7 s** | **114.1 s** |

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

### Idiomatic SeqChain pipeline

Full CRISPRi promoter library: scan + interpret + score + annotate + filter +
distance + stream to JSON.  Pure generator chain, O(chunk) memory.

```
70,425 promoter-targeting guides
 5,443 unique genes targeted
65,001 perfectly specific (92.3%)

Genome load  : 1.48 s
mmap load    : 51 µs
warm_seeds   : 5.62 s
Pipeline drain: 69.0 s
Wall-clock    : 76.1 s
```

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        Python (PyO3)                        │
│  FmIndex.build()  ·  FmIndex.load()  ·  search_batch()     │
├──────────────┬──────────────┬───────────────────────────────┤
│  fm_index.rs │  persist.rs  │        simd_search.rs         │
│  BWT + SA    │  rkyv + mmap │  AVX2 vertical automaton      │
│  BlockRank   │  zero-copy   │  SSE2 LF-mapping              │
├──────────────┼──────────────┤  width-first BFS              │
│ kmer_index.rs│  strand.rs   │  seeded + unseeded dispatch   │
│ K=10, K=14   │  +/- decode  │  rayon parallel groups        │
│ seed tables  │              │                               │
└──────────────┴──────────────┴───────────────────────────────┘
```

**4,181 lines of Rust.  7 dependencies.  One purpose.**

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
        queries: list[str], mismatches: int = 0,
    ) -> tuple[list[int], list[int], list[bool], list[float]]: ...
    def search_batch_exhaustive(
        queries: list[str], mismatches: int = 0,
    ) -> tuple[list[int], list[int], list[bool], list[float]]: ...
    def chrom_names() -> list[str]: ...
    def chrom_ranges() -> list[tuple[int, int]]: ...
```

The GIL is released during `search_batch`, `build`, and `warm_seeds`.
Rayon parallelizes across query groups internally.

---

## SeqChain Integration

Needletail plugs into the SeqChain generator architecture as a `ScorerFn`
operation:

```python
from seqchain.operations.score.native_off_target import score_off_targets_native

# Standard SeqChain generators — nothing computed yet
hits    = regex_map(genome.sequences, preset.pam, topologies=genome.topologies)
guides  = interpret_guides(hits, genome.sequences, preset, topologies=genome.topologies)

# Native SIMD scorer replaces bowtie — same ScorerFn contract
scored  = score_off_targets_native(
    guides, index, genome.sequences, preset.pam,
    mismatches=2, chunk_size=65536,
)

# Continue the generator chain as usual
tiles     = anchor_features(genome.genes(), config, genome.chrom_lengths)
annotated = annotate_tracks(scored, tiles)
promoters = (g for g in annotated if g.tags.get("feature_type") == "promoter")
final     = append_tss_distance(promoters)

# Drain — computation starts here
stream_json_export("output.json", final)
```

The `chunk_size=65536` parameter controls the FFI amortization boundary:
large enough to keep the Rust SIMD cores saturated, small enough that memory
stays O(chunk).  At 64K guides per batch, the PyO3 serialization cost becomes
negligible relative to the ~30 s search time.

### Axis 8 compliance

The scorer is a pure streaming generator.  No `sorted()`, no `list()`, no
materialization beyond the current chunk.  Upstream generators (`regex_map`,
`interpret_guides`) already produce Regions in genomic coordinate order.  The
downstream `annotate_tracks` sweep-line consumes them in that order.  Memory
is O(65536 Regions) regardless of genome size.

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
qi, pos, strand, scores = idx.search_batch(spacers, mismatches=2)

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
| `anyhow` | Error handling |
| `tempfile` | Seed table build scratch |

Release builds use LTO for cross-crate inlining of the SIMD inner loops.

---

## License

MIT
