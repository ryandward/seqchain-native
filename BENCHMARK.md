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

**Timing:**

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

**Memory (single cycle):**

| Stage | RSS | Peak (HWM) | Notes |
|-------|-----|------------|-------|
| Pipeline start | 246 MB | 393 MB | FM-Index + K=10 seed tables resident |
| After PAM scan | 287 MB | 393 MB | +41 MB PAM hit positions |
| After enrich (944K) | 351 MB | 393 MB | +64 MB SoA GuideHits |
| After scoring | 1,240 MB | 2,474 MB | Transient BWT search frontiers |
| After 944K regions | 1,887 MB | 2,474 MB | +647 MB AoS Region vec |
| After annotation+sink | 1,772 MB | 2,474 MB | -115 MB: guides freed to disk |

**Multi-cycle accumulation (3 consecutive runs):**

| Cycle | RSS at start | RSS after drain | Peak HWM |
|-------|-------------|-----------------|----------|
| 1 | 246 MB | 1,772 MB | 2,474 MB |
| 2 | 1,861 MB | 1,867 MB | 3,254 MB |
| 3 | 1,889 MB | 1,858 MB | 3,254 MB |

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
