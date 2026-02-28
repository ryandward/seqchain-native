# Plan: needletail — N-API + Strand FM-Index Addon

## Context

`needletail` is a new native Node-API addon that bridges Node.js directly to a
`rust-bio` FM-Index via the Strand `SharedArrayBuffer` protocol. The alternative
approaches (`subprocess.Popen` with SAM stdout parsing, or WASM) impose either process
isolation overhead or a 32-bit / 4GB memory ceiling. The native addon loads into the
V8 process at runtime, memory-maps a full-size genome index in 64-bit space, and writes
alignment records zero-copy into a Strand ring buffer. The Node.js consumer reads results
through the existing `StrandView` API with no serialization.

---

## Architecture

```
┌──────────────────────────────────────────────────────┐
│  Node.js (main thread)                               │
│                                                      │
│  const view = new StrandView(sab)                    │
│  const cursor = view.allocateCursor()                │
│  await searcher.streamAlignments(sab, queries)  ──►  │ N-API AsyncTask
│       ↑ returns Promise<void>                        │ (libuv thread pool)
│  cursor.waitForRecord()  ◄── SharedArrayBuffer ──►   │
│  cursor.getU32('pos') / cursor.getF32('score')       │
└──────────────────────────────────────────────────────┘
                                │ SAB pointer
┌──────────────────────────────▼─────────────────────┐
│  Rust (libuv worker thread)                         │
│                                                     │
│  StrandProducer  ←→  AtomicI32 ctrl words           │
│  FmIndexSearcher (rust-bio BWT + suffix array)      │
│  writes u32/bool8/u8/f32 fields directly into SAB  │
└─────────────────────────────────────────────────────┘
```

---

## Alignment Record Schema

Fixed-width, **no heap region required**. All five fields are non-nullable scalars.

```typescript
buildSchema([
  { name: 'chrom_id',   type: 'u32'   },  // offset  0, width 4
  { name: 'pos',        type: 'u32'   },  // offset  4, width 4
  { name: 'strand',     type: 'bool8' },  // offset  8, width 1
  { name: 'mismatches', type: 'u8'    },  // offset  9, width 1
  { name: 'score',      type: 'f32'   },  // offset 10, width 4
])
// null bitmap: 0 bytes (no nullable fields)
// record_stride: ceil(14/4)*4 = 16 bytes
// heap_capacity: 0
```

These byte offsets are **mirrored verbatim as Rust constants** to prevent drift.
A startup assertion in the Rust code verifies them against the `schema_fp` in the
live SAB header (FNV-1a of the Strand binary schema descriptor).

---

## File Structure

```
needletail/
├── package.json          # napi-rs CLI devDep, build + test scripts
├── Cargo.toml            # [lib] crate-type = ["cdylib"]
├── build.rs              # napi-rs build script
├── index.js              # thin JS re-export of the .node file
├── src/
│   ├── lib.rs            # #[napi] class + async task entry
│   ├── strand.rs         # StrandProducer — exact Strand protocol in Rust
│   ├── fm_index.rs       # rust-bio FMIndex wrapper
│   └── error.rs          # StrandError, SearchError → napi::Error
└── test/
    ├── basic.test.ts     # integration: real FM-Index → StrandView consumer
    ├── synthetic.test.ts # perf: stub producer at max speed (no search)
    └── fixtures/
        └── ecoli_K12.fa  # small test genome (E. coli, ~4.6 MB)
```

---

## Rust Crate Dependencies (Cargo.toml)

```toml
[lib]
crate-type = ["cdylib"]

[dependencies]
napi        = { version = "2", features = ["async", "napi4"] }
napi-derive = "2"
bio         = "1"          # rust-bio: BWT, FM-Index, FASTA reader
atomic-wait = "1"          # Linux futex wait/wake for backpressure stall
anyhow      = "1"

[build-dependencies]
napi-build = "2"
```

---

## src/strand.rs — StrandProducer

### Constants (must mirror strand/src/constants.ts exactly)

```rust
const HEADER_SIZE: usize     = 512;
const STRAND_MAGIC: u32      = 0x5354_524E;
const STRAND_VERSION: u32    = 3;

// Int32Array indices for control words
const CTRL_WRITE_SEQ:   usize = 7;
const CTRL_COMMIT_SEQ:  usize = 8;
const CTRL_READ_CURSOR: usize = 9;
const CTRL_STATUS:      usize = 10;
const CTRL_ABORT:       usize = 13;

// Status values
const STATUS_STREAMING: i32 = 1;
const STATUS_EOS:       i32 = 2;
const STATUS_ERROR:     i32 = 3;

// Schema byte offsets (must match buildSchema() output above)
const OFF_CHROM_ID:   usize = 0;
const OFF_POS:        usize = 4;
const OFF_STRAND:     usize = 8;
const OFF_MISMATCHES: usize = 9;
const OFF_SCORE:      usize = 10;
const RECORD_STRIDE:  usize = 16;
```

### Struct

```rust
pub struct StrandProducer {
    ctrl:           NonNull<[AtomicI32; 128]>,
    index_base:     NonNull<u8>,
    index_capacity: usize,
    write_seq:      i32,
}
```

### Construction

1. Read magic at offset 0 (little-endian u32) → verify `STRAND_MAGIC`
2. Read version at offset 4 → verify `STRAND_VERSION == 3`
3. Verify header CRC: FNV-1a of bytes 0–23 must equal u32 at offset 24
4. Read `index_capacity` at offset 16; assert power-of-two
5. Read `record_stride` from header; assert equals `RECORD_STRIDE` (16)
6. Read schema_fp at offset 8; assert matches compile-time expected fingerprint
7. Cast `sab_ptr` to `&[AtomicI32; 128]`; compute `index_base = sab_ptr + HEADER_SIZE`

### `begin()`
```rust
self.ctrl()[CTRL_STATUS].store(STATUS_STREAMING, Release);
```

### `write_record(chrom_id, pos, strand, mismatches, score)`
1. **Abort check**: `ctrl[CTRL_ABORT].load(Acquire) == 1` → return `Err(Aborted)`
2. **Backpressure stall**: loop until `write_seq - ctrl[CTRL_READ_CURSOR].load(Acquire) < index_capacity as i32`.
   Use `atomic_wait::wait(&ctrl[CTRL_READ_CURSOR], current)` for zero-spin Linux futex.
3. **Slot**: `let slot = write_seq as usize & (index_capacity - 1)`
4. **Pointer**: `unsafe { index_base.add(slot * RECORD_STRIDE) }`
5. **Field writes** (all `write_unaligned`):
   - `+OFF_CHROM_ID`   → `u32::to_le_bytes(chrom_id)`
   - `+OFF_POS`        → `u32::to_le_bytes(pos)`
   - `+OFF_STRAND`     → `strand as u8`
   - `+OFF_MISMATCHES` → `mismatches`
   - `+OFF_SCORE`      → `f32::to_le_bytes(score)`
6. **Advance**: `write_seq += 1`
7. **Commit**: `ctrl[CTRL_COMMIT_SEQ].store(write_seq, Release)`
8. **Notify**: `atomic_wait::wake_one(&ctrl[CTRL_COMMIT_SEQ])` → wakes JS `Atomics.waitAsync`

### `finalize()`
```rust
ctrl[CTRL_STATUS].store(STATUS_EOS, Release);
atomic_wait::wake_all(&ctrl[CTRL_COMMIT_SEQ]);
```

### `abort()`
```rust
ctrl[CTRL_STATUS].store(STATUS_ERROR, Release);
atomic_wait::wake_all(&ctrl[CTRL_COMMIT_SEQ]);
```

---

## src/fm_index.rs — FmIndexSearcher

Uses `bio::data_structures::fmindex::FMIndex` with the DNA alphabet.

```rust
pub struct Alignment {
    pub chrom_id:   u32,
    pub pos:        u32,
    pub strand:     bool,
    pub mismatches: u8,
    pub score:      f32,
}

pub struct FmIndexSearcher {
    index:  FMIndex<BWT, Less, Occ>,
    sa:     SuffixArray,
    chroms: Vec<String>,   // chromosome name → id mapping
}

impl FmIndexSearcher {
    pub fn from_fasta(path: &str) -> anyhow::Result<Self>
    pub fn search(&self, query: &[u8]) -> Vec<Alignment>
}
```

`from_fasta`: reads all records via `bio::io::fasta::Reader`, concatenates with `$`
separator, builds BWT via `bio::alphabets::dna::revcomp` + `bio::data_structures::suffix_array`,
wraps in `FMIndex`. `search` uses `fm_index.interval(query)`, resolves suffix array positions,
maps back to chromosome + offset.

Score is a simple mismatch-count-based float: `score = 1.0 / (1.0 + mismatches as f32)`.

---

## src/lib.rs — N-API Bindings

```rust
#[napi]
pub struct NativeFMIndex {
    inner: Arc<FmIndexSearcher>,
}

#[napi]
impl NativeFMIndex {
    #[napi(constructor)]
    pub fn new(fasta_path: String) -> napi::Result<Self>

    /// Returns Promise<void>. Runs in libuv thread pool via AsyncTask.
    #[napi]
    pub async fn stream_alignments(
        &self,
        sab: Buffer,          // ArrayBuffer backing the SharedArrayBuffer
        queries: Vec<String>,
    ) -> napi::Result<()>
}
```

The async task body:
1. Construct `StrandProducer::new(sab.as_ptr(), sab.len())` → validate header
2. `producer.begin()`
3. For each query: `fm_index.search(query.as_bytes())` → for each alignment: `producer.write_record(...)` → on `Err(Aborted)` break cleanly
4. `producer.finalize()` on success, `producer.abort()` on panic/error

---

## package.json (key fields)

```json
{
  "name": "needletail",
  "version": "0.1.0",
  "main": "index.js",
  "scripts": {
    "build":   "napi build --platform --release",
    "build:debug": "napi build --platform",
    "test":    "vitest run"
  },
  "devDependencies": {
    "@napi-rs/cli": "^2",
    "vitest": "^2",
    "@strand/core": "file:../strand"
  }
}
```

---

## test/basic.test.ts (integration)

```typescript
import { buildSchema, computeStrandMap, initStrandHeader, StrandView } from '@strand/core';
import { NativeFMIndex } from '../index.js';

const schema = buildSchema([
  { name: 'chrom_id',   type: 'u32'   },
  { name: 'pos',        type: 'u32'   },
  { name: 'strand',     type: 'bool8' },
  { name: 'mismatches', type: 'u8'    },
  { name: 'score',      type: 'f32'   },
]);
const map = computeStrandMap({ schema, indexCapacity: 65536, heapCapacity: 0 });
const sab = new SharedArrayBuffer(map.total_bytes);
initStrandHeader(sab, map);

const searcher = new NativeFMIndex('./test/fixtures/ecoli_K12.fa');
const view = new StrandView(sab);
const cursor = view.allocateCursor();

await searcher.streamAlignments(Buffer.from(sab), ['ATGATGATGATGATGATGATG']);

let count = 0;
while (await cursor.waitForRecord()) {
  expect(cursor.getU32('chrom_id')).toBeGreaterThanOrEqual(0);
  expect(cursor.getF32('score')).toBeGreaterThan(0);
  count++;
  cursor.advance();
}
expect(count).toBeGreaterThan(0);
```

**Cancel test**: Set `Atomics.store(ctrl, 13, 1)` mid-stream, verify `STATUS_ERROR` and
promise resolves (no hang).

**Backpressure test**: Pause consumer for 100ms, verify producer stalls (COMMIT_SEQ
stops advancing) and resumes when consumer unblocks.

---

## Critical Reference Files

| File | Why |
|---|---|
| `/home/ryan.ward/Git/strand/src/constants.ts` | Canonical byte offsets and ctrl word indices — Rust mirrors these verbatim |
| `/home/ryan.ward/Git/strand/src/writer.ts` | Claim/commit/notify cycle to replicate in `strand.rs` |
| `/home/ryan.ward/Git/strand/src/view.ts` | Consumer API used in integration tests |
| `/home/ryan.ward/Git/strand/src/header.ts` | FNV-1a CRC logic to port to Rust for header validation |
| `/home/ryan.ward/Git/strand/src/schema.ts` | Schema fingerprint computation — to verify at Rust startup |

---

## Build + Verification Steps

```bash
# 0. Prerequisite: build @strand/core so dist/ exists (no npm publish needed —
#    needletail references it via file:../strand in package.json)
cd ~/Git/strand && npm run build && cd ~/Git/needletail

# 1. Install JS tooling
npm install

# 2. Compile Rust addon
npm run build          # → needletail.linux-x64-gnu.node

# 3. Run integration tests
npm test

# 4. Verify backpressure: slow consumer doesn't OOM or deadlock
# 5. Verify abort: Atomics.store(ctrl, 13, 1) mid-stream resolves cleanly
```
