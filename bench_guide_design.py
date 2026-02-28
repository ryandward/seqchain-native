"""Real-world benchmark: genome-wide CRISPR guide design on SacCer3.

Replicates the core computational bottleneck from
SeqChain/examples/design_custom_yeast_cas12m.py:

  1. Scan genome for all TTN PAM sites (Cas12m)
  2. Extract 20bp spacers
  3. Search ALL spacers against the genome for off-targets (0–3 mismatches)
  4. PAM-validate each hit

This is the step that currently uses bowtie1 via subprocess.
"""

import re
import time
import sys
from pathlib import Path
from collections import Counter

from needletail import FmIndex

FASTA = str(Path(__file__).parent / "test/fixtures/sacCer3.fa")
PAM = "TT"          # TTN PAM — we match the TT prefix
SPACER_LEN = 20
PAM_LEN = 3         # TTN = 3 bases


# ═════════════════════════════════════════════════════════════════════════════
#  Step 1: PAM scan + spacer extraction
# ═════════════════════════════════════════════════════════════════════════════

def read_fasta(path):
    """Minimal FASTA reader — returns list of (name, sequence)."""
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


RC_TABLE = str.maketrans("ACGT", "TGCA")

def revcomp(seq):
    return seq.translate(RC_TABLE)[::-1]


def scan_pam_sites(chroms):
    """Scan all chromosomes for TTN PAM sites, extract 20bp spacers.

    Cas12m: PAM is downstream of spacer.
      + strand: ...SPACER(20bp)TTN...  → spacer = seq[pos-20:pos]
      - strand: ...NAA SPACER_RC(20bp)... → spacer = revcomp(seq[pos+3:pos+23])
    """
    guides = []
    for chrom_name, seq in chroms:
        seq_len = len(seq)

        # Forward strand: find TT at each position
        for i in range(SPACER_LEN, seq_len - PAM_LEN + 1):
            if seq[i] == 'T' and seq[i+1] == 'T':
                spacer = seq[i - SPACER_LEN : i]
                if 'N' not in spacer:
                    guides.append(spacer)

        # Reverse strand: find AA (complement of TT on reverse strand)
        # On - strand, PAM is upstream in the reference: NAA...spacer_rc
        for i in range(0, seq_len - SPACER_LEN - PAM_LEN + 1):
            if seq[i+1] == 'A' and seq[i+2] == 'A':
                spacer_region = seq[i + PAM_LEN : i + PAM_LEN + SPACER_LEN]
                if 'N' not in spacer_region:
                    guides.append(revcomp(spacer_region))

    return guides


# ═════════════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("Genome-wide CRISPR guide design benchmark — SacCer3 × Cas12m")
    print("=" * 70)
    print()

    # ── Load genome ───────────────────────────────────────────────────────
    print("Loading genome...", flush=True)
    t0 = time.perf_counter()
    chroms = read_fasta(FASTA)
    total_bp = sum(len(s) for _, s in chroms)
    t_load = time.perf_counter() - t0
    print(f"  {len(chroms)} chromosomes, {total_bp:,} bp in {t_load:.2f}s")

    # ── PAM scan ──────────────────────────────────────────────────────────
    print("\nScanning for TTN PAM sites + extracting 20bp spacers...", flush=True)
    t0 = time.perf_counter()
    all_spacers = scan_pam_sites(chroms)
    t_scan = time.perf_counter() - t0

    unique_spacers = list(set(all_spacers))
    print(f"  {len(all_spacers):,} total guides")
    print(f"  {len(unique_spacers):,} unique spacers")
    print(f"  Scanned in {t_scan:.2f}s")

    # ── Build FM-Index + persist ────────────────────────────────────────
    index_path = FASTA.rsplit(".", 1)[0] + ".seqchain"
    print(f"\nBuilding FM-Index + saving to {index_path}...", flush=True)
    t0 = time.perf_counter()
    idx = FmIndex.build(FASTA, index_path)
    t_build = time.perf_counter() - t0
    print(f"  Built + saved in {t_build:.2f}s")

    # ── Reload via mmap + warm seeds ─────────────────────────────────────
    print("\nLoading via mmap + warming seeds...", flush=True)
    t0 = time.perf_counter()
    idx = FmIndex.load(index_path)
    t_load_mmap = time.perf_counter() - t0
    t0 = time.perf_counter()
    idx.warm_seeds(index_path)
    t_warm = time.perf_counter() - t0
    print(f"  mmap load: {t_load_mmap*1e6:.0f}µs, warm_seeds: {t_warm:.2f}s")

    # ── Off-target search ─────────────────────────────────────────────────
    # Process in batches to manage memory.
    BATCH_SIZE = 10_000

    for mm in [0, 1, 2, 3]:
        print(f"\nSearching {len(unique_spacers):,} unique spacers × {mm}mm...", flush=True)
        t0 = time.perf_counter()

        total_hits = 0
        batches = 0
        for start in range(0, len(unique_spacers), BATCH_SIZE):
            batch = unique_spacers[start : start + BATCH_SIZE]
            qi, pos, strand, scores = idx.search_batch(batch, mismatches=mm)
            total_hits += len(qi)
            batches += 1

        t_search = time.perf_counter() - t0

        per_query_us = t_search / len(unique_spacers) * 1e6
        print(f"  {total_hits:,} total hits in {t_search:.2f}s")
        print(f"  {per_query_us:.1f} µs/query | {batches} batches of {BATCH_SIZE}")
        print(f"  {total_hits / len(unique_spacers):.1f} avg hits/spacer")

        # Don't continue to higher mm if this is already slow
        if t_search > 300:
            print(f"  (skipping higher mismatch levels — already {t_search:.0f}s)")
            break

    # ── Summary ───────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Genome       : SacCer3 ({total_bp:,} bp, {len(chroms)} chroms)")
    print(f"  Nuclease     : Cas12m (TTN PAM, {SPACER_LEN}bp spacer)")
    print(f"  Total guides : {len(all_spacers):,}")
    print(f"  Unique spacers: {len(unique_spacers):,}")
    print(f"  Index build  : {t_build:.2f}s")
    print(f"  mmap load    : {t_load_mmap*1e6:.0f}µs")
    print(f"  warm_seeds   : {t_warm:.2f}s")
    print(f"  PAM scan     : {t_scan:.2f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
