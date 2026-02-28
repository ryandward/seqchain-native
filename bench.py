"""Benchmark: width-first SIMD FM-Index vs bowtie1 vs bowtie2 on SacCer3."""

import time
import subprocess
import tempfile
import os
from needletail import FmIndex

FASTA = os.path.expanduser("~/Git/SeqChain/tests/data/saccer3/sacCer3.fa")

# Test queries — real 20bp sequences from SacCer3 chrI
QUERIES = [
    "CCACACCACACCCACACACC",   # unique hit
    "ATGATGATGATGATGATGAT",   # repeat-rich
    "TTAGGGTTAGGGTTAGGGTT",   # telomeric repeat
    "GCGATCGCGATCGCGATCGC",   # GC-rich
    "AATTAATTAATTAATTAATT",   # AT-rich repeat
    "ACTGGCACTGGCACTGGCAC",   # moderate complexity
    "NNNNNNNNNNNNNNNNNNNN",   # all-N (should get 0 hits)
    "CCACACCACACCCACACACC",   # duplicate of first (test dedup)
]

def run_bowtie(version, index_path, queries, mismatches):
    """Run bowtie1 or bowtie2 and return (wall_ms, hit_count)."""
    bowtie_bin = f"/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie{'2' if version == 2 else ''}"
    if not os.path.exists(bowtie_bin):
        return None, None

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for i, q in enumerate(queries):
            f.write(f">q{i}\n{q}\n")
        query_file = f.name

    try:
        if version == 2:
            cmd = [bowtie_bin, "-f", "-N", str(mismatches), "-x", index_path, query_file]
        else:
            cmd = [bowtie_bin, "-f", "-v", str(mismatches), index_path, query_file]

        start = time.perf_counter_ns()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        elapsed_ns = time.perf_counter_ns() - start

        hits = len([l for l in result.stdout.strip().split('\n') if l and not l.startswith('@')])
        return elapsed_ns / 1e6, hits
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None, None
    finally:
        os.unlink(query_file)

def build_bowtie_index(version, fasta, out_prefix):
    """Build a bowtie index if it doesn't exist."""
    bowtie_build = f"/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie{'2' if version == 2 else ''}-build"
    ext = ".1.bt2" if version == 2 else ".1.ebwt"
    if os.path.exists(out_prefix + ext):
        return True
    if not os.path.exists(bowtie_build):
        return False
    subprocess.run([bowtie_build, fasta, out_prefix], capture_output=True, timeout=120)
    return os.path.exists(out_prefix + ext)

def main():
    print(f"Genome: {FASTA}")
    print(f"Queries: {len(QUERIES)} × {len(QUERIES[0])}bp")
    print()

    # ── Build FM-Index ─────────────────────────────────────────────────────
    print("Building FM-Index...")
    t0 = time.perf_counter()
    idx = FmIndex(FASTA)
    build_s = time.perf_counter() - t0
    chroms = idx.chrom_names()
    print(f"  Built in {build_s:.2f}s — {len(chroms)} chromosomes")
    print()

    # ── Benchmark: needletail (width-first) ──────────────────────────
    print("=" * 70)
    print("needletail (width-first SIMD FM-Index)")
    print("=" * 70)

    for mm in [0, 1, 2, 3]:
        # Warmup
        idx.search_batch(QUERIES, mismatches=mm)

        # Timed run (average of 100 iterations for stable timing)
        iters = 100
        t0 = time.perf_counter_ns()
        for _ in range(iters):
            qi, pos, scores = idx.search_batch(QUERIES, mismatches=mm)
        elapsed_ns = time.perf_counter_ns() - t0
        avg_us = elapsed_ns / iters / 1000

        print(f"  {mm}mm: {len(qi):>6} hits | {avg_us:>8.1f} µs/batch | {avg_us/len(QUERIES):>6.1f} µs/query")

    # ── Scale test: 1K queries ────────────────────────────────────────────
    print()
    print("Scale test: 1000 queries × 3mm")
    queries_1k = QUERIES * 125  # 8 × 125 = 1000
    # Warmup
    idx.search_batch(queries_1k, mismatches=3)

    iters = 10
    t0 = time.perf_counter_ns()
    for _ in range(iters):
        qi, pos, scores = idx.search_batch(queries_1k, mismatches=3)
    elapsed_ns = time.perf_counter_ns() - t0
    avg_ms = elapsed_ns / iters / 1e6
    print(f"  {len(qi):>6} hits | {avg_ms:>8.2f} ms/batch | {avg_ms/len(queries_1k)*1000:>6.1f} µs/query")

    # ── Benchmark: bowtie1 ────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("bowtie1")
    print("=" * 70)

    bt1_prefix = "/tmp/saccer3_bt1"
    if build_bowtie_index(1, FASTA, bt1_prefix):
        for mm in [0, 1, 2, 3]:
            ms, hits = run_bowtie(1, bt1_prefix, QUERIES, mm)
            if ms is not None:
                print(f"  {mm}mm: {hits:>6} hits | {ms:>8.1f} ms (wall, incl. process launch + index load)")
    else:
        print("  bowtie1 not available")

    # ── Benchmark: bowtie2 ────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("bowtie2")
    print("=" * 70)

    bt2_prefix = "/tmp/saccer3_bt2"
    if build_bowtie_index(2, FASTA, bt2_prefix):
        for mm in [0, 1]:
            ms, hits = run_bowtie(2, bt2_prefix, QUERIES, mm)
            if ms is not None:
                print(f"  {mm}mm: {hits:>6} hits | {ms:>8.1f} ms (wall, incl. process launch + index load)")
    else:
        print("  bowtie2 not available")

    # ── Summary ───────────────────────────────────────────────────────────
    print()
    print("Note: bowtie times include process launch + index load (~20-40ms).")
    print("needletail times are pure in-process search (index already loaded).")

if __name__ == "__main__":
    main()
