"""Test persistence: build → save → mmap load → search correctness."""

from needletail import FmIndex
import time
import tempfile
import os

genome = os.path.join(os.path.dirname(__file__), "..", "test", "fixtures", "phiX174.fa")

# Build from FASTA (baseline — has seeds)
print("Building from FASTA...")
idx1 = FmIndex(genome)
print(f"  chroms: {idx1.chrom_names()}")

# Build + save to disk (has seeds)
with tempfile.NamedTemporaryFile(suffix=".seqchain", delete=False) as f:
    idx_path = f.name

print(f"Building + saving to {idx_path}...")
idx2 = FmIndex.build(genome, idx_path)
file_size = os.path.getsize(idx_path)
print(f"  .seqchain file: {file_size:,} bytes")

# Load from mmap (no seeds — uses unseeded search path)
print("Loading from mmap...")
t = time.perf_counter()
idx3 = FmIndex.load(idx_path)
load_us = (time.perf_counter() - t) * 1e6
print(f"  mmap load: {load_us:.0f}µs")

# Correctness: all three must return identical results
# Use sequences from phiX174
queries = [
    "GAGTTTTATCGCTTCCA",  # 17bp from phiX174 start
    "ATGAAAGCAATTTTCGT",  # another phiX sequence
]

print()
for mm in [0, 1, 2]:
    r1 = idx1.search_batch(queries, mismatches=mm)
    r2 = idx2.search_batch(queries, mismatches=mm)
    r3 = idx3.search_batch(queries, mismatches=mm)

    # Sort results for comparison (hit order may differ)
    def sort_hits(r):
        return sorted(zip(r[0], r[1], r[2]))

    s1 = sort_hits(r1)
    s2 = sort_hits(r2)
    s3 = sort_hits(r3)

    assert s1 == s2, f"Built vs Build+Save differ at mm={mm}!\n  built: {s1}\n  saved: {s2}"
    assert s1 == s3, f"Built vs Loaded differ at mm={mm}!\n  built: {s1}\n  loaded: {s3}"
    print(f"  mm={mm}: {len(r1[0])} hits, all 3 paths match")

# Chrom names and ranges should match too
assert idx1.chrom_names() == idx2.chrom_names() == idx3.chrom_names()
assert idx1.chrom_ranges() == idx2.chrom_ranges() == idx3.chrom_ranges()
print(f"  chrom metadata matches")

os.unlink(idx_path)
# Clean up kmer index if created
kmer_idx = idx_path.replace(".seqchain", f".10mer.idx")
if os.path.exists(kmer_idx):
    os.unlink(kmer_idx)

print("\nPASS")
