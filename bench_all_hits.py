"""Fair comparison: all-alignments mode for bowtie1/2 vs needletail."""

import time
import subprocess
import tempfile
import os
from needletail import FmIndex

FASTA = os.path.expanduser("~/Git/SeqChain/tests/data/saccer3/sacCer3.fa")

QUERIES = [
    "CCACACCACACCCACACACC",
    "ATGATGATGATGATGATGAT",
    "TTAGGGTTAGGGTTAGGGTT",
    "GCGATCGCGATCGCGATCGC",
    "AATTAATTAATTAATTAATT",
    "ACTGGCACTGGCACTGGCAC",
]

def make_query_file(queries):
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    for i, q in enumerate(queries):
        f.write(f">q{i}\n{q}\n")
    f.close()
    return f.name

def run_bowtie1_all(index_prefix, query_file, mm):
    cmd = ["/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie",
           "-f", "-v", str(mm), "-a", index_prefix, query_file]
    start = time.perf_counter_ns()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    elapsed = (time.perf_counter_ns() - start) / 1e6
    hits = len([l for l in r.stdout.strip().split('\n') if l])
    return elapsed, hits

def run_bowtie2_all(index_prefix, query_file, mm):
    cmd = ["/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie2",
           "-f", "-N", str(mm), "-a", "--no-hd", "-x", index_prefix, query_file]
    start = time.perf_counter_ns()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    elapsed = (time.perf_counter_ns() - start) / 1e6
    hits = len([l for l in r.stdout.strip().split('\n') if l and not l.startswith('@')])
    return elapsed, hits

def main():
    print(f"SacCer3 | {len(QUERIES)} queries × {len(QUERIES[0])}bp | ALL alignments mode")
    print()

    # Build FM-Index
    idx = FmIndex(FASTA)

    # Build bowtie indexes
    bt1_prefix = "/tmp/saccer3_bt1"
    bt2_prefix = "/tmp/saccer3_bt2"
    if not os.path.exists(bt1_prefix + ".1.ebwt"):
        subprocess.run(["/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie-build",
                        FASTA, bt1_prefix], capture_output=True, timeout=120)
    if not os.path.exists(bt2_prefix + ".1.bt2"):
        subprocess.run(["/home/ryan.ward/miniforge3/envs/bowtie/bin/bowtie2-build",
                        FASTA, bt2_prefix], capture_output=True, timeout=120)

    qfile = make_query_file(QUERIES)

    # Precompute byte arrays
    queries_fwd = [q.upper() for q in QUERIES]

    print(f"{'mm':>3} | {'needletail':>30} | {'bowtie1 -a':>30} | {'bowtie2 -a':>30}")
    print(f"{'':>3} | {'hits':>8}  {'µs':>8}  {'µs/q':>8} | {'hits':>8}  {'ms':>10} | {'hits':>8}  {'ms':>10}")
    print("-" * 100)

    for mm in [0, 1, 2, 3]:
        # needletail: warmup + timed
        idx.search_batch(QUERIES, mismatches=mm)
        iters = 50
        t0 = time.perf_counter_ns()
        for _ in range(iters):
            qi, pos, scores = idx.search_batch(QUERIES, mismatches=mm)
        elapsed_ns = time.perf_counter_ns() - t0
        sc_us = elapsed_ns / iters / 1000
        sc_hits = len(qi)
        sc_per_q = sc_us / len(QUERIES)

        # bowtie1
        bt1_ms, bt1_hits = run_bowtie1_all(bt1_prefix, qfile, mm)

        # bowtie2
        bt2_ms, bt2_hits = run_bowtie2_all(bt2_prefix, qfile, mm)

        print(f"  {mm} | {sc_hits:>8}  {sc_us:>8.1f}  {sc_per_q:>8.1f} | {bt1_hits:>8}  {bt1_ms:>8.1f} ms | {bt2_hits:>8}  {bt2_ms:>8.1f} ms")

    os.unlink(qfile)

    # ── Scale comparison: 1000 queries ────────────────────────────────────
    print()
    print("Scale: 1000 queries × 3mm (needletail only — bowtie too slow for all-hits)")
    queries_1k = QUERIES * 167  # ~1000
    # Warmup
    idx.search_batch(queries_1k[:100], mismatches=3)

    t0 = time.perf_counter_ns()
    qi, pos, scores = idx.search_batch(queries_1k, mismatches=3)
    elapsed_ms = (time.perf_counter_ns() - t0) / 1e6
    print(f"  {len(qi)} hits | {elapsed_ms:.1f} ms | {elapsed_ms/len(queries_1k)*1000:.1f} µs/query")

if __name__ == "__main__":
    main()
