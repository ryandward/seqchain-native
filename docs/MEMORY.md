# needletail Memory

## Project Overview
PyO3 Python extension for width-first SIMD FM-Index search (CRISPR guide off-target).
Located at `/home/ryan.ward/Git/needletail/`.

## Key Architecture Facts

### V8 Atomics ↔ Rust Futex Incompatibility (CRITICAL)
V8's `Atomics.notify` uses an **internal user-space futex** (pthread condition variables in a
linked-list waiter queue). It does NOT issue a `FUTEX_WAKE` kernel syscall. Therefore:
- `atomic_wait::wait` (kernel `FUTEX_WAIT`) **cannot be woken by JS `Atomics.notify`**.
- `atomic_wait::wake_one` **cannot wake JS `Atomics.waitAsync`**.
- **Solution**: Use `std::thread::sleep(Duration::from_micros(50))` in backpressure loops.
- Producer wakes JS consumer via `CTRL_STATUS` + `CTRL_COMMIT_SEQ` visible via `Atomics.load`.

### Strand Schema Byte Offsets (corrected vs plan)
buildSchema alignment rules give f32 at offset **12** (not 10 as the plan stated):
- chrom_id u32 → offset 0
- pos      u32 → offset 4
- strand   bool8 → offset 8
- mismatches u8 → offset 9
- score    f32 → offset **12** (aligned from 10, +2 pad for 4-byte align)
- record_stride = 16

### InterleavedRank (replaced bio's Occ)
- bio's `Occ` was `Vec<Vec<usize>>` — 85 inner heap allocs, 8 scattered cache accesses per LF-map
- **InterleavedRank**: flat `Vec<u32>`, layout `data[pos*4+base_idx]` (A=0,C=1,G=2,T=3)
- LF-mapping loads 2 contiguous 16-byte blocks = 2 cache lines (was 8 scattered)
- Memory: 4 × bwt_len × 4 bytes. SacCer3: ~195 MB (was ~8 GB with bio's Occ at k=1)
- `lf_map(l, r) -> [(u32, u32); 4]` — batch method on FmOcc trait
- `prefetch_lf(l, r)` — inter-group `_mm_prefetch(_MM_HINT_T0)` for next (l,r) pair
- bio's `bwt()`, `less()`, `suffix_array()` still used for construction; `FMIndex`/`Occ` dropped entirely

### bio 1.6.0 API (still used for construction only)
- `suffix_array(&text) -> Vec<usize>`
- `bwt(&text, &sa) -> Vec<u8>` and `less(&bwt, &alphabet) -> Vec<usize>`
- `fasta::Reader::from_file(path) -> anyhow::Result<Self>`
- Record iterator items: `io::Result<Record>`

### napi 2.x API
- `Undefined` (type alias for `()`) — not `JsUndefined`
- `resolve()` returns `Ok(())` for void tasks
- `Task` trait: `type JsValue = Undefined` for Promise<void>
- Exported class name: `NativeFmIndex` (napi-rs converts FMIndex → FmIndex)
- `NativeFMIndex` exported as type alias in `.d.ts`

## Dependency Pins (rustc 1.79.0)
- `napi-build` pinned to 2.1.0 (2.3.1 requires rustc 1.88)
- `indexmap` pinned to 2.7.0 (2.13.0 requires rustc 1.82)

## Build Commands
```bash
cd ~/Git/strand && npm run build     # prerequisite
cd ~/Git/needletail
npm install
npm run build    # → needletail.linux-x64-gnu.node
npm test         # vitest run — 16 tests pass
```

## Test Files
- `test/synthetic.test.ts` — pure JS StrandWriter throughput (no native addon)
- `test/basic.test.ts` — integration: NativeFmIndex + StrandView consumer
- `test/saccer3.test.ts` — real SacCer3 genome (12MB, 17 chroms); 8 tests incl. mismatch
- `test/fixtures/test_genome.fa` — synthetic 6-line FASTA with many ATGATG patterns

## Current Schema (record_stride = 20)
```typescript
buildSchema([
  { name: 'chrom',      type: 'utf8_ref' },  // offset  0 — intern table handle
  { name: 'pos',        type: 'u32'      },  // offset  4
  { name: 'query_id',   type: 'u32'      },  // offset  8
  { name: 'strand',     type: 'bool8'    },  // offset 12
  { name: 'mismatches', type: 'u8'       },  // offset 13
  { name: 'score',      type: 'f32'      },  // offset 16 (pad 2 bytes from 14)
])
```
- `chrom` is `utf8_ref` (u32 intern handle). Consumer calls `view.updateInternTable(chromNames)`.
- Chromosome names written to SAB v4 metadata: `initStrandHeader(sab, map, { intern: chromNames })`.
- `query_id` identifies which input query produced each hit (0-indexed).
- `streamAlignments(sab, queries, maxMismatches?)` — optional 3rd arg (0–3, default 0).
- `chromNames()` — returns FASTA-order NCBI accession names (e.g. NC_001133.9 = chrI).

### Inexact Backward Search Bug (FIXED)
`Occ::get(bwt, r, a)` uses **inclusive** upper bound (counts BWT[0..=r]).
Calling `occ(bwt_len, a)` panics — the checkpoint vec has exactly `floor((n-1)/k)+1` entries.
Fix: use inclusive bounds internally, initial `r = sa.len() - 1`, convert to exclusive only at output.

## Strand Version
Strand is at v4. `strand.rs` must have `STRAND_VERSION: u32 = 4`.
v3→v4 change: metadata region in header tail (unused by needletail producer; harmless).

## Backpressure in Tests
The backpressure test uses a **drain loop** that polls `CTRL_COMMIT_SEQ` every 2ms
and advances `CTRL_READ_CURSOR` via `Atomics.store + Atomics.notify`. The producer
wakes within 50µs (spin-sleep) when the cursor advances.
