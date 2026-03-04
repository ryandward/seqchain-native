# Benchmarking Needletail

## Quick Start

```bash
# Build release server
cargo build --release -p needletail-server

# Run all core tests
cargo test -p needletail-core

# Run the server (stderr shows memory profile)
./target/release/needletail-server 2>server.log &

# Upload a genome and run a job
curl -X POST http://localhost:8002/api/files/upload -F "file=@test/fixtures/sacCer3.gb"
curl -X POST http://localhost:8002/api/methods/design_library \
  -H "Content-Type: application/json" \
  -d '{"genome": "<genome_id>", "preset": "spcas9", "mismatches": 2}'

# Watch the memory profile in real time
tail -f server.log
```

## Architecture

The pipeline has two phases with fundamentally different memory profiles:

1. **Heavy Alignment** (BWT search frontiers) — the tolerated RAM budget.
   Transient: allocates during search, freed after scoring completes.
2. **Sweep-line Annotation → FileSink** — O(k) auxiliary memory where k
   is the active overlap window (typically < 10 features). Guides drain
   directly to a JSON file on disk. No intermediate Vec.

```
scored_guides.into_iter()  →  SweepAnnotator  →  FileSink (disk)
     (consumed)                  O(k)              O(buffer)
```

## Memory Profiling

The pipeline (`design_library()` in `crates/needletail-core/src/pipeline/design.rs`)
has built-in per-stage memory instrumentation using `/proc/self/status`.

### How it works

Two inline helpers read kernel-maintained process stats:

```rust
fn peak_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status").ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmHWM:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<f64>().ok()))
        .unwrap_or(0.0) / 1024.0
}
fn rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status").ok()
        .and_then(|s| s.lines().find(|l| l.starts_with("VmRSS:"))
            .and_then(|l| l.split_whitespace().nth(1))
            .and_then(|v| v.parse::<f64>().ok()))
        .unwrap_or(0.0) / 1024.0
}
```

- **VmRSS** — current resident set size (physical RAM in use right now)
- **VmHWM** — high water mark (peak RSS ever reached by this process)

Then `eprintln!` after each stage:

```rust
eprintln!("[mem] after scoring: RSS={:.0}MB peak={:.0}MB", rss_mb(), peak_mb());
```

### Reading the output

Run the server and submit a job. stderr shows the stage-by-stage profile:

```
[mem] pipeline start:                     RSS=246MB   peak=393MB
[mem] after PAM scan:                     RSS=287MB   peak=393MB
[mem] after enrich (944K guides):         RSS=351MB   peak=393MB
[mem] after scoring:                      RSS=1240MB  peak=2474MB  ← transient BWT frontiers
[mem] after 944K guide regions:           RSS=1887MB  peak=2474MB
[mem] after annotation+sink (986K):       RSS=1772MB  peak=2474MB  ← RSS drops: guides freed to disk
```

**VmHWM is monotonically increasing.** Comparing it between stages tells you
exactly which stage pushed the peak. A jump in peak but not RSS means the
stage allocated transiently and freed.

### Adding your own checkpoints

Add anywhere inside `design_library()`:

```rust
eprintln!("[mem] after <your label>: RSS={:.0}MB peak={:.0}MB", rss_mb(), peak_mb());
```

Linux-only (reads `/proc/self/status`). No dependencies. Works in Docker.

## Reference Profiles

### sacCer3 (S. cerevisiae, 12 Mbp, 17 chroms, mm=2)

#### Index Build Memory Profile (fresh process)

Index build happens once per upload. Held in RAM; no disk write.

| Stage | Time | RSS | Peak (HWM) | Notes |
|-------|------|-----|------------|-------|
| Process start | — | 38 MB | 68 MB | baseline |
| Genome parse (GenBank) | 0.07s | 93 MB | 93 MB | +55 MB master buffer |
| SA-IS suffix array | 1.11s | 186 MB | **393 MB** | ← dominant peak: u64 SA + working space |
| BWT + less table | 0.04s | 198 MB | 393 MB | SA freed after BWT built |
| BlockRank (189K blocks) | 0.07s | 209 MB | 393 MB | ~12 MB of rank structures |
| K=10 seed table | 0.29s | 279 MB | 393 MB | 8 MB seed + 50 MB pos table |
| **Upload total** | **1.57s** | **279 MB** | **393 MB** | index resident, ready |

SA-IS is the memory spike: rust-bio uses `u64` indices (8 bytes/position),
so the suffix array alone is 12M × 8 = 96 MB. SA-IS needs O(N) working space
on top of that, pushing the transient peak to 393 MB. After SA-IS completes
the working space is freed; BWT is built from the SA and the SA is dropped.

The resident index at rest: genome text 55 MB + BWT 12 MB + BlockRank 12 MB
+ K=10 seed/pos tables 58 MB ≈ **137 MB net** (RSS=279 MB includes glibc overhead).

---

#### Baseline — before O(C) recomposition

**Timing (baseline):**

| Stage | Time |
|-------|------|
| Genome parse (GenBank) | 0.14s |
| FM-Index (SA-IS + BWT + BlockRank) | 2.83s |
| Seed tier K=10 (SA sweep) | 0.76s |
| **Upload total** | **3.73s** |
| PAM scan + enrich (944K guides) | ~0.1s |
| Score spacers (893K unique, mm=2) | **2.50s** |
| Build 944K Regions | ~0.1s |
| Sweep annotation → FileSink (986K written) | ~0.1s |
| **Pipeline total** | **~3.0s** |

**Memory — baseline (single cycle):**

| Stage | RSS | Peak (HWM) | Notes |
|-------|-----|------------|-------|
| Pipeline start | 246 MB | 393 MB | FM-Index + K=10 seed tables resident |
| After PAM scan | 287 MB | 393 MB | +41 MB PAM hit positions |
| After enrich (944K) | 351 MB | 393 MB | +64 MB SoA GuideHits |
| After scoring | 1,240 MB | **2,474 MB** | ← DAM: transient BWT frontiers for all 893K |
| After 944K regions | 1,887 MB | **2,474 MB** | ← DAM: +647 MB AoS Region vec |
| After annotation+sink | 1,772 MB | **2,474 MB** | -115 MB: guides freed to disk |

**Multi-cycle — baseline (3 consecutive runs):**

| Cycle | RSS at start | RSS after drain | Peak HWM |
|-------|-------------|-----------------|----------|
| 1 | 246 MB | 1,772 MB | 2,474 MB |
| 2 | 1,861 MB | 1,867 MB | 3,254 MB |
| 3 | 1,889 MB | 1,858 MB | 3,254 MB |

---

#### After O(C) recomposition — chunked scoring + lazy Regions

Two architectural changes applied to `pipeline/design.rs`:

1. **O(C) chunked scoring** — `score_spacers` called 9 times with 100K spacers
   each instead of once with 893K. BWT frontier peak bounded by `C × h̄ × 13B`
   per cycle instead of `N × h̄ × 13B`.

2. **Lazy Region iterator** — `Vec<Region>` never materialised. A
   `Vec<usize>` sort index (7.5 MB) replaces the 647 MB Region vec. Regions
   are constructed one at a time inside the SweepAnnotator closure.

**Timing (recomposed):**

| Stage | Time | Δ |
|-------|------|---|
| PAM scan + enrich (944K guides) | ~0.1s | — |
| Score spacers (893K unique, mm=2, 9 chunks) | **0.54s** | **-1.96s  (-78%)** |
| Index sort (944K guides) | < 0.1s | replaces +0.1s Region build |
| Sweep annotation → FileSink (975K written) | ~0.1s | — |
| **Pipeline total** | **~0.8s** | **-2.2s (-73%)** |

**Memory — recomposed (Cycle 1, fresh process):**

| Stage | RSS | Peak (HWM) | Δ vs baseline |
|-------|-----|------------|---------------|
| Pipeline start | 248 MB | 393 MB | index resident from upload |
| After PAM scan | 289 MB | 393 MB | +41 MB — same as baseline |
| After enrich (944K) | 353 MB | 393 MB | +64 MB — same as baseline |
| Scoring chunk 1 (100K) | 641 MB | 641 MB | ← per-chunk peak replaces 2,474 MB DAM |
| Scoring chunk 9/9 | 677 MB | **677 MB** | **peak HWM -1,797 MB vs baseline** |
| After scoring | 645 MB | 673 MB | freed chunk; no merge spike |
| After index sort | 645 MB | 673 MB | +7.5 MB indices (not +647 MB Regions) |
| After annotation+sink | 950 MB | **950 MB** | **-822 MB vs baseline 1,772 MB** |

**Multi-cycle — recomposed (2 consecutive runs):**

| Cycle | RSS at start | RSS after drain | Peak HWM | Δ HWM vs baseline |
|-------|-------------|-----------------|----------|-------------------|
| 1 | 248 MB | 950 MB | 950 MB | **-1,524 MB (-62%)** |
| 2 | 1,711 MB | 2,088 MB | 1,908 MB | **-1,346 MB (-41%)** |

Cycle 2 RSS increase is glibc allocator freelist retention (same behaviour as
before — not a leak). The per-cycle peak HWM remains bounded by the chunk
size, not by N.

**Summary of Δ:**

| Metric | Baseline | Recomposed | Δ |
|--------|----------|------------|---|
| Peak HWM (single cycle) | 2,474 MB | **950 MB** | **-1,524 MB (-62%)** |
| Peak HWM (multi-cycle) | 3,254 MB | **1,908 MB** | **-1,346 MB (-41%)** |
| RSS after annotation | 1,772 MB | **950 MB** | **-822 MB (-46%)** |
| Region Vec allocation | 647 MB | **0 MB** | **-647 MB (-100%)** |
| Scoring time | 2.50 s | **0.61 s** | **-1.89 s (-76%)** |

The 78% scoring speedup is attributable to two effects:
- **Cache locality**: 100K spacers fit in L2/L3 better than 893K; the BWT
  frontier nodes stay hot.
- **Merge elimination**: `run_search_seeded`'s serial merge step is O(C × h̄)
  per chunk instead of O(N × h̄) for the full batch.

Cycle 2-3 RSS flatlines — this is the glibc allocator freelist retaining pages,
not a leak. The allocator avoids expensive kernel syscalls on subsequent runs.

### ZM56TN (Zymomonas, 11.9 Mbp, 45 contigs, mm=2)

| Stage | RSS | Peak (HWM) |
|-------|-----|------------|
| Pipeline start | 234 MB | 234 MB |
| After PAM scan | 272 MB | 272 MB |
| After enrich (921K) | 334 MB | 334 MB |
| After scoring | 478 MB | 2,445 MB |
| After 921K regions | 1,959 MB | 2,445 MB |
| After annotation+sink | ~320 MB | 2,445 MB |

---

### Head-to-head: SpCas9 / SacCer3 / mm=2 (same hardware, same run)

Both servers measured simultaneously. `/proc/<pid>/status` sampled at each
stage transition. Old = port 8002 (pre-recomposition); new = port 8003
(O(C) chunked scoring + lazy Region iterator).

**Stage-by-stage RSS trace (Δ from job-start RSS):**

| Stage | OLD ΔRSS | NEW ΔRSS | Note |
|-------|----------|----------|------|
| Tiling genome features | +0 MB | +0 MB | — |
| Aligning 893K spacers | +3 MB | **+0 MB** | chunked; freelist absorbs oscillation |
| Building scored regions | +48 → +124 MB | **stage eliminated** | Vec\<Region\> removed |
| Annotating & streaming | +124 MB | **+0 MB** | lazy iterator; no Region Vec |
| Streamed (done) | -129 MB | -186 MB | output freed |

**ΔHWM per job (incremental peak above pre-job waterline):**

| | OLD | NEW | Δ |
|-|-----|-----|---|
| ΔHWM | **+135 MB** | **+0 MB** | **-135 MB (-100%)** |

New code produces **zero net new peak allocation** above what is already
loaded in the process. Every allocation fits within the existing freelist.

**Wall time (warm process, index resident):**

| | OLD | NEW | Δ |
|-|-----|-----|---|
| Wall time | **4.07 s** | **4.00 s** | **-0.07 s (-2%)** |

**Fresh-process peak HWM** (from server `eprintln` telemetry, first run on
a clean process — see "After O(C) recomposition" section above):

| | OLD | NEW | Δ |
|-|-----|-----|---|
| Fresh peak HWM | **2,474 MB** | **950 MB** | **-1,524 MB (-62%)** |

---

## Seed Tier Sizing

K=10 seed table for sacCer3 (12 Mbp):

```
[sa_sweep] K=10: 1,017,260 distinct k-mers, 12,156,952 positions
  seed table: 8.0 MB    (4^10 entries × 8 bytes)
  pos table:  50.4 MB   (12M positions × ~4 bytes)
  total:      58.4 MB
```

K=14 was benchmarked and dropped: only 8% faster scoring for mm=3,
but costs +293 MB RAM and +1.3s upload time.

## Test Fixtures

- `test/fixtures/sacCer3.fa` — S. cerevisiae genome (12 Mbp)
- `test/fixtures/sacCer3.gb` — S. cerevisiae GenBank (12 Mbp)
- `test/fixtures/sacCer3.seqchain` — pre-built index
- `test/fixtures/sacCer3.10mer.idx` — pre-built K=10 seed table
- `test/fixtures/phiX174.fa` — small test genome (5 kb)
- `test/fixtures/zymomonas.gb` — Zymomonas GenBank (4.3 MB file)
- `test/fixtures/synthetic_genome.gb` — synthetic GenBank for unit tests

## Python Bindings Benchmarks

```bash
source .venv/bin/activate
maturin develop --release

# Correctness
python3 tests/test_persist.py
python3 tests/test_persist_saccer3.py

# vs Bowtie1 (hit counts + timing)
python3 tests/bench_vs_bowtie.py

# Full SeqChain pipeline
python3 bench_annotate_saccer3.py -o /tmp/bench_output.json --mm 2
```
