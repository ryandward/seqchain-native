# needletail — Assessment and SeqChain Integration

## What it is

needletail is a Node.js N-API addon that builds a rust-bio FM-Index from a FASTA
file and streams alignment records into a Strand SharedArrayBuffer. It exposes one
constructor (`NativeFmIndex(fasta)`) and one async method (`streamAlignments`). Results
land in a Strand ring buffer and are consumed by the JS caller through `StrandView`.

---

## Comparison with bowtie1 and bowtie2

Benchmarked on SacCer3 (11.7 MB, 17 chromosomes). bowtie1 v1.3.1, bowtie2 v2.5.4.

### Index build

| Tool            | Build time | Index size | Persists? |
|-----------------|-----------|------------|-----------|
| needletail | 1,110 ms  | ~300 MB RAM | No (rebuilt every process start) |
| bowtie1         | 4,468 ms  | 20 MB disk  | Yes (.ebwt files) |
| bowtie2         | 5,418 ms  | 25 MB disk  | Yes (.bt2 files) |

needletail builds faster than either bowtie tool, but the index lives only in RAM
and must be rebuilt on every process restart. bowtie pre-builds once and memory-maps the
index at search time.

### Search — exact match (0 mismatches)

| Query        | Hits | needletail | bowtie1 | bowtie2 |
|--------------|------|-----------------|---------|---------|
| unique 20 nt | 1    | 0.09 ms         | 22 ms   | 39 ms   |
| repeat (Ty)  | 80   | 0.42 ms         | 24 ms   | 45 ms   |

### Search — 1 mismatch

| Query        | Hits | needletail | bowtie1 | bowtie2 |
|--------------|------|-----------------|---------|---------|
| unique 20 nt | 1    | 0.14 ms         | 31 ms   | 27 ms   |
| repeat (Ty)  | 108  | 0.14 ms         | 27 ms   | 41 ms   |

Hit counts are identical across all three tools.

**Reading the bowtie numbers:** bowtie's wall time is dominated by process launch and
index load from disk (~20–25 ms), not the actual BWT search. The search itself is
sub-millisecond in all tools. needletail's advantage is that it runs in-process —
no fork, no exec, no disk I/O after the index is built.

---

## The philosophy tension

SeqChain's design principle is **streaming-only, minimal in-memory**. Data flows from
S3 to SeqChain to S3; neither side materialises the full payload. The FM-Index approach
violates this directly:

- The suffix array for SacCer3 is ~96 MB (12M positions × 8 bytes).
- The Occ table at sampling rate 3 is ~192 MB.
- Total in-memory index: ~300 MB — **25× the genome size**.
- For human (3 GB genome): ~75 GB. Cannot run on a standard instance.

This is not a fixable inefficiency; it is the nature of FM-Index. A full-text index
requires materialising the entire genome before the first query. bowtie manages this by
pre-building and memory-mapping, which amortises the I/O cost but does not change the
RAM requirement.

---

## What needletail can realistically be

**It fits the SeqChain model for small, bounded genomes** — bacteria (~5 MB),
yeast (~12 MB), small viruses. The 300 MB index fits comfortably in any modern server.
The process starts once, builds the index, and answers thousands of queries per second
with sub-millisecond latency per query. The Strand ring buffer means results stream
zero-copy to the consumer as they are found.

**It does not fit for large genomes** — human, mouse, maize. Those require a
persistent, memory-mapped index (bowtie's .bt2 files or a similar format) so the OS
page cache can hold only the used portions. That is a separate design problem.

**The right framing for SeqChain integration:** needletail is the alignment engine
for the CRISPR library design and transposon off-target workflows in SeqChain. Those
workflows already target bacterial and yeast genomes. The tool is in scope for those
use cases and out of scope for mammalian whole-genome alignment.

---

## What bowtie has that we don't

- Index persistence (no 1.1 s rebuild on every server start)
- MAPQ multi-mapping quality scores
- Gapped alignment (insertions and deletions)
- SAM/BAM output with CIGAR strings
- `--very-sensitive` / `--local` alignment modes

The most important gap for production use is **index persistence**. Serialising the
in-memory FM-Index to a file and reloading it without rebuilding would cut cold-start
time from 1,100 ms to ~50 ms (disk read). This is on the roadmap but not yet built.

---

## Strand integration gaps (current status)

Three gaps were identified in `docs/strand-integration.md`. The status after this
session's work:

| Gap | Description | Status |
|-----|-------------|--------|
| 1 | `chrom_id: u32` not resolvable to names | **Fixed** — `chrom: utf8_ref` + metadata channel |
| 2 | Output unsorted, `findFirst()` unavailable | Open — requires sorting before write |
| 3 | No `query_id` field, multi-query results unlabelled | **Fixed** — `query_id: u32` added |
| 4 | `acknowledgeRead()` never called | Open — safe for current workload; documented |

Gap 2 (sorting) requires collecting all hits for a chromosome before writing, which
temporarily breaks the streaming model. The right approach is to sort within each
chromosome's hits and write them in order, flagging `pos` with `FIELD_FLAG_SORTED_ASC`.
Left for a future change.
