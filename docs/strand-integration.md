# Strand Integration — Current State and Gaps

This document assesses how needletail uses `@strand/core`, what parts of
the Strand feature surface are being wasted, and what is genuinely missing.

---

## What We Currently Send

The alignment record schema is five dumb scalars with no flags:

```typescript
buildSchema([
  { name: 'chrom_id',   type: 'u32'   },  // integer index, not a name
  { name: 'pos',        type: 'u32'   },  // unsorted, SA order
  { name: 'strand',     type: 'bool8' },
  { name: 'mismatches', type: 'u8'    },
  { name: 'score',      type: 'f32'   },
])
```

No heap region, no nullable fields, no field flags, no sorting. This uses
Strand as a plain ring buffer and ignores most of its feature surface.

---

## Three Concrete Gaps in Our Usage

### 1. `chrom_id: u32` is the wrong type

Strand has `utf8_ref` precisely for this situation: a `u32` handle into an
intern table that the consumer holds. The correct design is:

```typescript
buildSchema([
  { name: 'chrom',      type: 'utf8_ref' },  // handle → chromosome name
  { name: 'pos',        type: 'u32', flags: FIELD_FLAG_SORTED_ASC },
  { name: 'strand',     type: 'bool8' },
  { name: 'mismatches', type: 'u8'    },
  { name: 'score',      type: 'f32'   },
])
```

The consumer populates the intern table once before processing:

```typescript
const view = new StrandView(sab);
view.updateInternTable(['NC_001133.9', 'NC_001134.8', /* … */]);

cursor.getRef('chrom');  // → 'NC_001133.9'
```

Right now anyone consuming our stream receives `chrom_id=3` with no way to
resolve that to `NC_001136.10` (chrIV) without out-of-band communication.
The intern table mechanism exists in Strand specifically to avoid this.

### 2. Output is unsorted, so `findFirst()` is dead code for us

Strand exports `FIELD_FLAG_SORTED_ASC`. If records were emitted sorted by
chromosome then position and the field tagged accordingly, a genomic browser
could binary-search into the result set:

```typescript
// Jump directly to the first hit at or after position 100,000 — O(log n)
const seq = view.findFirst('pos', 100_000);
cursor.seek(seq);
```

FM-index search returns hits in suffix array order (essentially arbitrary).
Without sorting before writing and flagging the field, `findFirst()` produces
undefined results and is useless to consumers. For a query returning thousands
of hits across a full genome, linear scanning is the only option today.

### 3. No `query_id` field

When multiple queries go into a single `streamAlignments` call, the results
arrive interleaved with no attribution. There is no way to tell which of the
81 records in a `[UNIQUE_QUERY, REPEAT_QUERY]` call came from which query.
A `query_id: u32` field added to the schema would fix this at zero cost (4
bytes per record).

---

## What Strand Is Genuinely Missing

### No metadata channel

To populate the `utf8_ref` intern table the consumer needs the chromosome name
list before the first record arrives. Strand has no mechanism to carry
stream-level metadata alongside records — the header carries only the schema
descriptor and geometry, not application data.

The practical workaround: expose chromosome names as a separate property on the
`NativeFmIndex` JS class (populated during index construction from the FASTA
headers), and let the caller pass them to `view.updateInternTable()` before
calling `streamAlignments`. This is not a Strand fix — it is an API contract
the caller must uphold.

---

## Latent Deadlock: `acknowledgeRead()` is Never Called

Our tests pass because `index_capacity = 65536` swallows the largest result set
comfortably (80 hits for the Ty repeat query). But the ring is a fixed-size
circular buffer. A repetitive query on a larger genome returning more hits than
`index_capacity` will stall the Rust producer waiting for the ring to drain —
which it never will because the JS consumer is blocked in `await
streamAlignments()` and never calls `acknowledgeRead()`.

The correct consumer pattern for large result sets:

```typescript
const streamPromise = searcher.streamAlignments(Buffer.from(sab), queries);

const drainLoop = async () => {
  let acked = 0;
  while (true) {
    const committed = view.committedCount;
    if (committed > acked) {
      // process records acked..committed here
      view.acknowledgeRead(committed);
      acked = committed;
    }
    if (view.status === 'eos' || view.status === 'error') break;
    await new Promise(r => setTimeout(r, 2));
  }
};

await Promise.all([streamPromise, drainLoop()]);
```

For needletail's current workload (small genomes, exact matches, bounded
hit counts) this has not mattered. It will matter for any genome larger than
roughly `index_capacity × RECORD_STRIDE` bytes of results.

---

## Feature Inventory

| Strand feature | Used | Notes |
|---|---|---|
| Scalar fields (u32, f32, u8, bool8) | Yes | All five record fields |
| `utf8_ref` + `updateInternTable()` | No | Chromosome names sent as opaque integers |
| `FIELD_FLAG_SORTED_ASC` + `findFirst()` | No | Output is in SA order; binary search unavailable |
| `query_id` attribution | No | Multi-query results are unlabelled |
| `acknowledgeRead()` | No | Latent deadlock for large result sets |
| Heap fields (utf8, json, bytes, arrays) | No | Not needed for current schema |
| Nullable fields | No | All fields are always present |
| Multi-consumer (`registerConsumer`) | No | Not needed yet |
| `signalAbort()` | No | Abort path tested producer-side only |
| `waitForCommit()` | No | Caller blocks on full Promise resolution |

---

## Recommended Next Steps

1. **Add `query_id: u32` to the schema** — zero-cost, unblocks multi-query consumers.
2. **Expose chromosome names from `NativeFmIndex`** — return the ordered list of FASTA
   record IDs from the constructor so callers can populate the intern table.
3. **Switch `chrom_id` to `utf8_ref`** — consumers get names, not integers.
4. **Sort by chrom then pos before writing; tag with `FIELD_FLAG_SORTED_ASC`** —
   enables `findFirst()` and makes output compatible with coordinate-sorted
   genomic tooling.
5. **Document the `acknowledgeRead()` contract** — callers must drain the ring
   for result sets larger than `index_capacity`.
