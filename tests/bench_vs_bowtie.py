"""Benchmark: needletail (mmap) vs bowtie1 subprocess."""
import subprocess
import time
import tempfile
import os

genome = os.path.join(os.path.dirname(__file__), "..", "test", "fixtures", "sacCer3.fa")

if not os.path.exists(genome):
    print(f"SKIP: {genome} not found (run: gunzip -k test/fixtures/sacCer3.fa.gz)")
    exit(0)

# --- needletail ---
from needletail import FmIndex

# Cold start: build + save
t = time.perf_counter()
idx = FmIndex.build(genome, "/tmp/sacCer3.seqchain")
build_ms = (time.perf_counter() - t) * 1000

# Warm start: mmap load
t = time.perf_counter_ns()
idx2 = FmIndex.load("/tmp/sacCer3.seqchain")
load_us = (time.perf_counter_ns() - t) / 1000

# Search: unique hit
queries_unique = ["TATTTATACCCATTCCCTCA"]
queries_repeat = ["TGGGATTCCATTGTTGATAA"]

print("=" * 70)
print("needletail")
print("=" * 70)

for label, queries in [("unique", queries_unique), ("repeat", queries_repeat)]:
    for mm in [0, 1, 2]:
        t = time.perf_counter()
        r = idx2.search_batch(queries, mismatches=mm)
        search_us = (time.perf_counter() - t) * 1e6
        print(f"  {label} mm={mm}: {len(r[0]):>4} hits in {search_us:>8.1f}us")

print(f"  build: {build_ms:.0f}ms, mmap load: {load_us:.0f}us")

# --- bowtie1 ---
print()
print("=" * 70)
print("bowtie1")
print("=" * 70)

bt_idx = "/tmp/sacCer3_bt"

# Check for bowtie
bowtie_build = None
bowtie = None
for path in ["/usr/bin/bowtie-build", os.path.expanduser("~/miniforge3/envs/bowtie/bin/bowtie-build")]:
    if os.path.exists(path):
        bowtie_build = path
        bowtie = path.replace("-build", "")
        break

# Try PATH
if bowtie_build is None:
    try:
        subprocess.run(["bowtie-build", "--version"], capture_output=True, check=True)
        bowtie_build = "bowtie-build"
        bowtie = "bowtie"
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass

if bowtie_build:
    # Build bowtie index
    t = time.perf_counter()
    subprocess.run([bowtie_build, "--quiet", genome, bt_idx], check=True)
    bt_build_ms = (time.perf_counter() - t) * 1000

    for label, queries in [("unique", queries_unique), ("repeat", queries_repeat)]:
        for mm in [0, 1, 2]:
            qfa = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
            for i, q in enumerate(queries):
                qfa.write(f">q{i}\n{q}\n")
            qfa.close()

            sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
            sam.close()

            t = time.perf_counter()
            subprocess.run(
                [bowtie, "-f", "-a", "-v", str(mm), bt_idx, qfa.name, sam.name],
                check=True, capture_output=True,
            )
            bt_ms = (time.perf_counter() - t) * 1000

            with open(sam.name) as f:
                hits = sum(1 for line in f if not line.startswith("@"))

            print(f"  {label} mm={mm}: {hits:>4} hits in {bt_ms:>8.1f}ms")
            os.unlink(qfa.name)
            os.unlink(sam.name)

    print(f"  build: {bt_build_ms:.0f}ms")

    # Cleanup
    for ext in [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]:
        try:
            os.unlink(bt_idx + ext)
        except OSError:
            pass
else:
    print("  bowtie1 not available (install bowtie to compare)")

os.unlink("/tmp/sacCer3.seqchain")
kmer_idx = "/tmp/sacCer3.10mer.idx"
if os.path.exists(kmer_idx):
    os.unlink(kmer_idx)
