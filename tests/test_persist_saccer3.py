"""Test persistence on SacCer3: build → save → mmap load → correctness + timing."""

from needletail import FmIndex
import time
import tempfile
import os
import sys
import traceback
import threading

genome = os.path.join(os.path.dirname(__file__), "..", "test", "fixtures", "sacCer3.fa")

if not os.path.exists(genome):
    print(f"SKIP: {genome} not found (run: gunzip -k test/fixtures/sacCer3.fa.gz)")
    exit(0)


def timed(label):
    """Context manager that prints elapsed time and a spinner."""
    class Timer:
        def __init__(self):
            self.elapsed_ms = 0
            self._stop = threading.Event()
            self._thread = None

        def __enter__(self):
            self._start = time.perf_counter()
            self._stop.clear()
            self._thread = threading.Thread(target=self._spin, daemon=True)
            self._thread.start()
            return self

        def _spin(self):
            phases = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
            i = 0
            while not self._stop.wait(1.0):
                elapsed = time.perf_counter() - self._start
                print(f"\r  {phases[i % len(phases)]} {label}... {elapsed:.1f}s", end="", flush=True)
                i += 1

        def __exit__(self, *exc):
            self._stop.set()
            self._thread.join()
            self.elapsed_ms = (time.perf_counter() - self._start) * 1000
            print(f"\r  ✓ {label}: {self.elapsed_ms:.0f}ms          ")
            return False

    return Timer()


def main():
    errors = []

    # ── Step 1: Build from FASTA ──────────────────────────────────────
    print("=" * 60)
    print("STEP 1: Build FM-Index from FASTA (sacCer3)")
    print("=" * 60)
    try:
        with timed("Building FM-Index + seed tiers") as t1:
            idx1 = FmIndex(genome)
        build_ms = t1.elapsed_ms
        print(f"  chroms: {idx1.chrom_names()}")
    except Exception as e:
        traceback.print_exc()
        errors.append(f"Build failed: {e}")
        print("FAIL (cannot continue)")
        sys.exit(1)

    # ── Step 2: Build + Save ──────────────────────────────────────────
    print()
    print("=" * 60)
    print("STEP 2: Build + Save to .seqchain")
    print("=" * 60)
    with tempfile.NamedTemporaryFile(suffix=".seqchain", delete=False) as f:
        idx_path = f.name
    try:
        with timed("Build + save") as t2:
            idx2 = FmIndex.build(genome, idx_path)
        save_ms = t2.elapsed_ms
        file_size = os.path.getsize(idx_path)
        print(f"  file: {idx_path}")
        print(f"  size: {file_size:,} bytes ({file_size / 1e6:.1f} MB)")
    except Exception as e:
        traceback.print_exc()
        errors.append(f"Build+Save failed: {e}")
        print("FAIL (cannot continue)")
        sys.exit(1)

    # ── Step 3: mmap Load ─────────────────────────────────────────────
    print()
    print("=" * 60)
    print("STEP 3: Load from mmap")
    print("=" * 60)
    try:
        t = time.perf_counter_ns()
        idx3 = FmIndex.load(idx_path)
        load_us = (time.perf_counter_ns() - t) / 1000
        print(f"  ✓ mmap load: {load_us:.0f}µs")
    except Exception as e:
        traceback.print_exc()
        errors.append(f"mmap load failed: {e}")
        print("FAIL (cannot continue)")
        sys.exit(1)

    # ── Step 4: Correctness (seeded vs unseeded) ──────────────────────
    print()
    print("=" * 60)
    print("STEP 4: Correctness — seeded (idx1, idx2) vs unseeded (idx3)")
    print("=" * 60)
    queries = ["TATTTATACCCATTCCCTCA", "TGGGATTCCATTGTTGATAA"]
    print(f"  queries: {queries}")

    def sort_hits(r):
        return sorted(zip(r[0], r[1], r[2]))

    for mm in [0, 1, 2]:
        try:
            r1 = idx1.search_batch(queries, mismatches=mm)
            r2 = idx2.search_batch(queries, mismatches=mm)
            r3 = idx3.search_batch(queries, mismatches=mm)
            s1, s2, s3 = sort_hits(r1), sort_hits(r2), sort_hits(r3)
            if s1 == s2 == s3:
                print(f"  ✓ mm={mm}: {len(r1[0])} hits, all 3 paths match")
            else:
                msg = f"mm={mm}: results differ!\n    s1={s1}\n    s2={s2}\n    s3={s3}"
                print(f"  ✗ {msg}")
                errors.append(msg)
        except Exception as e:
            traceback.print_exc()
            errors.append(f"mm={mm} search failed: {e}")

    # ── Step 5: warm_seeds on loaded index ────────────────────────────
    print()
    print("=" * 60)
    print("STEP 5: warm_seeds on mmap-loaded index")
    print("=" * 60)
    try:
        with timed("warm_seeds (K=10 + K=14)") as t5:
            idx3.warm_seeds(idx_path)
        warm_ms = t5.elapsed_ms
    except Exception as e:
        traceback.print_exc()
        errors.append(f"warm_seeds failed: {e}")
        warm_ms = -1

    # ── Step 6: Re-check with seeded search ───────────────────────────
    print()
    print("=" * 60)
    print("STEP 6: Correctness — seeded search after warm_seeds")
    print("=" * 60)
    for mm in [0, 1, 2]:
        try:
            r1 = idx1.search_batch(queries, mismatches=mm)
            r3 = idx3.search_batch(queries, mismatches=mm)
            s1, s3 = sort_hits(r1), sort_hits(r3)
            if s1 == s3:
                print(f"  ✓ mm={mm}: {len(r3[0])} hits, seeded path matches")
            else:
                msg = f"mm={mm} seeded differs after warm!\n    s1={s1}\n    s3={s3}"
                print(f"  ✗ {msg}")
                errors.append(msg)
        except Exception as e:
            traceback.print_exc()
            errors.append(f"mm={mm} seeded search failed: {e}")

    # ── Cleanup ───────────────────────────────────────────────────────
    os.unlink(idx_path)
    for suffix in [".10mer.idx", ".14mer.idx"]:
        p = idx_path.replace(".seqchain", suffix)
        if os.path.exists(p):
            os.unlink(p)

    # ── Summary ───────────────────────────────────────────────────────
    print()
    print("=" * 60)
    if errors:
        print(f"FAIL — {len(errors)} error(s):")
        for e in errors:
            print(f"  • {e}")
        sys.exit(1)
    else:
        print(f"Summary: build={build_ms:.0f}ms, save={save_ms:.0f}ms, "
              f"mmap_load={load_us:.0f}µs, warm_seeds={warm_ms:.0f}ms")
        print("PASS")


if __name__ == "__main__":
    main()
