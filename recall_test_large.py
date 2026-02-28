"""Larger recall test: 50K queries including edge cases."""

import time
import random
from pathlib import Path
from needletail import FmIndex

FASTA = str(Path.home() / "Git/SeqChain/tests/data/saccer3/sacCer3.fa")

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

RC_TABLE = str.maketrans("ACGT", "TGCA")
def revcomp(seq):
    return seq.translate(RC_TABLE)[::-1]

def scan_pam_sites(chroms, limit=50000):
    guides = set()
    for chrom_name, seq in chroms:
        seq_len = len(seq)
        for i in range(20, seq_len - 3):
            if seq[i] == 'T' and seq[i+1] == 'T':
                spacer = seq[i - 20 : i]
                if 'N' not in spacer:
                    guides.add(spacer)
                    if len(guides) >= limit:
                        return list(guides)
        for i in range(0, seq_len - 23):
            if seq[i+1] == 'A' and seq[i+2] == 'A':
                spacer_region = seq[i + 3 : i + 23]
                if 'N' not in spacer_region:
                    guides.add(revcomp(spacer_region))
                    if len(guides) >= limit:
                        return list(guides)
    return list(guides)

print("Loading genome + building index...")
idx = FmIndex(FASTA)
chroms = read_fasta(FASTA)
queries = scan_pam_sites(chroms, limit=50000)

# Add some adversarial homopolymer queries
for base in "ACGT":
    queries.append(base * 20)
# Add AT-repeat queries
queries.append("AT" * 10)
queries.append("TA" * 10)
queries.append("ATAT" * 5)

print(f"  {len(queries)} queries (including homopolymers)")

for mm in [0, 1, 2, 3]:
    t0 = time.perf_counter()
    qi_s, pos_s, sc_s = idx.search_batch(queries, mismatches=mm)
    t_seeded = time.perf_counter() - t0

    t0 = time.perf_counter()
    qi_e, pos_e, sc_e = idx.search_batch_exhaustive(queries, mismatches=mm)
    t_exhaust = time.perf_counter() - t0

    seeded = set(zip(qi_s, pos_s))
    exhaust = set(zip(qi_e, pos_e))

    missing = exhaust - seeded
    extra = seeded - exhaust
    recall = len(seeded & exhaust) / max(len(exhaust), 1) * 100

    print(f"\n{mm}mm: seeded={len(seeded):,}  exhaustive={len(exhaust):,}  "
          f"recall={recall:.1f}%  missing={len(missing)}  extra={len(extra)}")
    print(f"  seeded: {t_seeded:.2f}s  exhaustive: {t_exhaust:.2f}s")

    if missing and len(missing) <= 20:
        for qi, pos in sorted(missing)[:20]:
            q = queries[qi]
            print(f"    MISSING: qid={qi} pos={pos} query={q}")
    elif missing:
        print(f"    (showing first 5 of {len(missing)} missing)")
        for qi, pos in sorted(missing)[:5]:
            q = queries[qi]
            print(f"    MISSING: qid={qi} pos={pos} query={q}")
