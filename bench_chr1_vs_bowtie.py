"""Benchmark: needletail (BlockRank + mmap) vs bowtie1 — chr1 only.

CRISPR guide design pipeline on SacCer3 chromosome 1 (230 kb):
  1. Scan for TTN PAM sites, extract 20bp spacers
  2. Search all spacers for off-targets at 0–3 mismatches
  3. Compare hit counts and wall-clock time

Usage:
    python bench_chr1_vs_bowtie.py
"""

import os
import re
import subprocess
import tempfile
import time
from pathlib import Path

# ═════════════════════════════════════════════════════════════════════════════
#  Config
# ═════════════════════════════════════════════════════════════════════════════

FASTA = str(Path(__file__).parent / "test/fixtures/sacCer3_chr1.fa")
PAM = "TT"
SPACER_LEN = 20

RC_TABLE = str.maketrans("ACGT", "TGCA")

def revcomp(seq):
    return seq.translate(RC_TABLE)[::-1]


def read_fasta(path):
    chroms = []
    name = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    chroms.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        chroms.append((name, "".join(seq_parts).upper()))
    return chroms


def scan_pam_sites(chroms):
    guides = []
    for _, seq in chroms:
        seq_len = len(seq)
        # Forward: TTN downstream of spacer
        for i in range(SPACER_LEN, seq_len - 2):
            if seq[i] == "T" and seq[i + 1] == "T":
                spacer = seq[i - SPACER_LEN : i]
                if "N" not in spacer:
                    guides.append(spacer)
        # Reverse: NAA upstream
        for i in range(0, seq_len - SPACER_LEN - 2):
            if seq[i + 1] == "A" and seq[i + 2] == "A":
                spacer_region = seq[i + 3 : i + 3 + SPACER_LEN]
                if "N" not in spacer_region:
                    guides.append(revcomp(spacer_region))
    return guides


# ═════════════════════════════════════════════════════════════════════════════
#  needletail
# ═════════════════════════════════════════════════════════════════════════════

def bench_native(unique_spacers, mm_levels):
    from needletail import FmIndex

    idx_path = "/tmp/sacCer3_chr1.seqchain"

    # Build + save
    t0 = time.perf_counter()
    idx = FmIndex.build(FASTA, idx_path)
    t_build = time.perf_counter() - t0

    # mmap load
    t0 = time.perf_counter()
    idx_mmap = FmIndex.load(idx_path)
    t_load = time.perf_counter() - t0

    # Warm seeds on loaded index
    t0 = time.perf_counter()
    idx_mmap.warm_seeds(idx_path)
    t_warm = time.perf_counter() - t0

    results = {}
    for mm in mm_levels:
        t0 = time.perf_counter()
        qi, pos, strand, scores = idx_mmap.search_batch(unique_spacers, mismatches=mm)
        t_search = time.perf_counter() - t0
        per_query_us = t_search / len(unique_spacers) * 1e6
        results[mm] = {
            "hits": len(qi),
            "time_s": t_search,
            "per_query_us": per_query_us,
        }

    # Also verify built (BlockRank) path matches mmap (flat) path
    for mm in mm_levels:
        qi_b, pos_b, str_b, _ = idx.search_batch(unique_spacers, mismatches=mm)
        qi_m, pos_m, str_m, _ = idx_mmap.search_batch(unique_spacers, mismatches=mm)
        if len(qi_b) != len(qi_m):
            print(f"  WARNING: built vs loaded mismatch at mm={mm}: {len(qi_b)} vs {len(qi_m)}")

    os.unlink(idx_path)
    # Clean up seed table files
    for ext in [".10mer.idx", ".14mer.idx"]:
        p = idx_path.rsplit(".", 1)[0] + ext
        if os.path.exists(p):
            os.unlink(p)

    return {
        "build_s": t_build,
        "load_us": t_load * 1e6,
        "warm_s": t_warm,
        "results": results,
    }


# ═════════════════════════════════════════════════════════════════════════════
#  bowtie1
# ═════════════════════════════════════════════════════════════════════════════

def bench_bowtie(unique_spacers, mm_levels):
    bt_idx = "/tmp/sacCer3_chr1_bt"

    # Build bowtie index
    t0 = time.perf_counter()
    subprocess.run(
        ["bowtie-build", "--quiet", FASTA, bt_idx], check=True, capture_output=True
    )
    t_build = time.perf_counter() - t0

    # Write all spacers to a FASTA file
    qfa = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    for i, s in enumerate(unique_spacers):
        qfa.write(f">q{i}\n{s}\n")
    qfa.close()

    results = {}
    for mm in mm_levels:
        sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
        sam.close()

        t0 = time.perf_counter()
        subprocess.run(
            ["bowtie", "-f", "-a", "-v", str(mm), bt_idx, qfa.name, sam.name],
            check=True,
            capture_output=True,
        )
        t_search = time.perf_counter() - t0

        # Count hits
        with open(sam.name) as f:
            hits = sum(1 for line in f if not line.startswith("@"))

        per_query_us = t_search / len(unique_spacers) * 1e6
        results[mm] = {
            "hits": hits,
            "time_s": t_search,
            "per_query_us": per_query_us,
        }
        os.unlink(sam.name)

    os.unlink(qfa.name)
    for ext in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]:
        try:
            os.unlink(bt_idx + ext)
        except FileNotFoundError:
            pass

    return {"build_s": t_build, "results": results}


# ═════════════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("CRISPR Guide Design Benchmark — SacCer3 chr1 (230 kb)")
    print("needletail (BlockRank + mmap) vs bowtie1")
    print("=" * 72)

    # Load genome + scan PAMs
    chroms = read_fasta(FASTA)
    total_bp = sum(len(s) for _, s in chroms)
    all_spacers = scan_pam_sites(chroms)
    unique_spacers = sorted(set(all_spacers))

    print(f"\nGenome: {total_bp:,} bp")
    print(f"PAM sites: {len(all_spacers):,} total, {len(unique_spacers):,} unique spacers")

    mm_levels = [0, 1, 2, 3]

    # ── needletail ──────────────────────────────────────────────────
    print("\n--- needletail ---")
    native = bench_native(unique_spacers, mm_levels)
    print(f"  Build:      {native['build_s']:.2f}s")
    print(f"  mmap load:  {native['load_us']:.0f} µs")
    print(f"  warm_seeds: {native['warm_s']:.2f}s")
    for mm in mm_levels:
        r = native["results"][mm]
        print(f"  mm={mm}: {r['hits']:>8,} hits | {r['time_s']:.3f}s | {r['per_query_us']:.1f} µs/query")

    # ── bowtie1 ──────────────────────────────────────────────────────────
    print("\n--- bowtie1 ---")
    bt = bench_bowtie(unique_spacers, mm_levels)
    print(f"  Build: {bt['build_s']:.2f}s")
    for mm in mm_levels:
        r = bt["results"][mm]
        print(f"  mm={mm}: {r['hits']:>8,} hits | {r['time_s']:.3f}s | {r['per_query_us']:.1f} µs/query")

    # ── Comparison ───────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("Comparison")
    print("=" * 72)
    print(f"{'':>6} {'native hits':>12} {'bowtie hits':>12} {'match':>6}  {'native':>10} {'bowtie':>10} {'speedup':>8}")
    print("-" * 72)
    for mm in mm_levels:
        rn = native["results"][mm]
        rb = bt["results"][mm]
        hit_match = "OK" if rn["hits"] == rb["hits"] else "DIFF"
        speedup = rb["per_query_us"] / rn["per_query_us"] if rn["per_query_us"] > 0 else float("inf")
        print(
            f"mm={mm}: {rn['hits']:>12,} {rb['hits']:>12,} {hit_match:>6}"
            f"  {rn['per_query_us']:>8.1f}µs {rb['per_query_us']:>8.1f}µs {speedup:>7.1f}x"
        )

    print(f"\nIndex build: native {native['build_s']:.2f}s vs bowtie {bt['build_s']:.2f}s")
    print(f"mmap load: {native['load_us']:.0f} µs")
    print("=" * 72)


if __name__ == "__main__":
    main()
