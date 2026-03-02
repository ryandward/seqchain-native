# Needletail Engine Contract

How Needletail implements the GenomeHub engine interface. This document
describes the mapping from the generic contract (defined in GenomeHub's
`docs/engine-interface.md`) to Needletail's concrete implementation.

## Endpoint mapping

| GenomeHub contract | Needletail route | Implementation |
|--------------------|-----------------|----------------|
| `GET /api/health` | `GET /api/health` | `routes/health.rs` |
| `GET /api/methods` | `GET /api/methods` | `routes/methods.rs` — `method_catalog()` |
| `GET /api/methods/:id` | `GET /api/methods/{method_id}` | `routes/methods.rs` — `get_method()` |
| `POST /api/methods/:id` | `POST /api/methods/{method_id}` | `routes/methods.rs` — `dispatch_method()` |
| `POST /api/files/upload` | `POST /api/files/upload` | `routes/files.rs` — multipart form, field `file` |
| `GET /api/tracks/:id/data` | `GET /api/jobs/{id}/stream` | `routes/jobs.rs` — streams result JSON |

### Additional endpoints (engine-specific)

| Route | Purpose |
|-------|---------|
| `GET /api/jobs/{id}` | Poll job status and progress |
| `DELETE /api/jobs/{id}` | Cancel a running job (best-effort, always 200) |
| `GET /api/presets` | List preset categories |
| `GET /api/presets/{category}` | Full detail for all presets in a category |
| `POST /api/genomes/upload` | Legacy JSON-body upload (server-side path) |
| `GET /api/genomes/` | List uploaded genomes |
| `GET /api/genomes/{id}` | Genome metadata |

## Method catalog generation

The method catalog (`GET /api/methods`) is **programmatically generated**
from the preset registry. No biological knowledge is hardcoded in the
catalog — it queries `PresetRegistry::select_options(category)` to build
the `options` arrays that GenomeHub renders as dropdowns.

```
presets/recognition/*.yaml  ──┐
                              ├──► PresetRegistry ──► method_catalog()
presets/tiling/*.yaml       ──┘         │
                                        ├──► GET /api/presets/{category}
                                        └──► GET /api/methods (select options)
```

### How a select parameter is built

A method parameter declared as `type: "select"` is populated at runtime:

1. The method declares which **preset category** the parameter draws from
   (e.g., `"recognition"` or `"tiling"`).
2. `PresetRegistry::select_options(category)` scans all embedded YAMLs in
   that category and builds `[{"value": id, "label": name, "description": ...}]`.
3. The first preset in the category becomes the `default`.
4. GenomeHub receives this in the method schema and renders a `<select>` dropdown.

Adding a new preset = one YAML file + one line in `ALL_PRESETS`. The method
catalog, API routes, and GenomeHub UI update automatically.

## Preset system

### Philosophy

Presets are organized by **axis of variation**, not by biological concept:

| Category | What it parameterizes | Examples |
|----------|----------------------|----------|
| `recognition` | How to find target sites — motif, flanking length, directionality | SpCas9, Cas12a, Cas12m, (future: Himar1, Tn5) |
| `tiling` | How to partition the genome into named regions based on gene proximity | saccer3 (yeast), (future: ecoli_k12, human_hg38) |

These axes are **orthogonal**. A pipeline run composes one preset from each
axis: `recognition × tiling`. No combinatorial explosion — you don't need
`spcas9_saccer3`, `spcas9_ecoli`, etc.

### YAML format

All presets follow SeqChain's YAML format — the same files work in both
systems. The `type` field in the YAML distinguishes categories:

**Recognition preset** (`presets/recognition/spcas9.yaml`):
```yaml
type: crispr
name: SpCas9
pam: NGG
spacer_len: 20
pam_direction: downstream
description: >
  Streptococcus pyogenes Cas9. Recognizes NGG PAM downstream of a 20bp
  spacer. The most widely used CRISPR nuclease.
```

**Tiling preset** (`presets/tiling/saccer3_features.yaml`):
```yaml
type: feature
organism: "S. cerevisiae"
features:
  gene_body:
    relation: overlap
    priority: 1
    anchor: five_prime
  promoter:
    relation: upstream
    max_distance: 500
    priority: 2
    anchor: five_prime
  terminator:
    relation: downstream
    max_distance: 200
    priority: 3
    anchor: three_prime
default_feature: intergenic
```

### Embedding strategy

Presets are embedded at compile time via `include_str!`. The binary is
self-contained — no runtime file access needed. The YAML files in `presets/`
remain the canonical source; the Rust code is a consumer, not an author.

### Registry internals

`PresetRegistry` (in `needletail-core/src/models/preset.rs`) provides:

| Method | Returns |
|--------|---------|
| `categories()` | All known category names |
| `list(category)` | Preset IDs in a category |
| `yaml(category, id)` | Raw YAML source for a preset |
| `select_options(category)` | `[{value, label, description}]` for method catalog |
| `detail(category)` | Full parsed YAML with injected `id` field |

The typed structs (`CRISPRPreset`, `FeatureConfig`) delegate to the registry
for lookup:

```rust
CRISPRPreset::by_name("spcas9")
// → PresetRegistry::yaml("recognition", "spcas9")
// → serde_yaml::from_str → CRISPRPreset

FeatureConfig::by_name("saccer3")
// → PresetRegistry::yaml("tiling", "saccer3")
// → serde_yaml::from_str → FeatureConfig
```

## Async job protocol

Needletail methods are async. The dispatch flow is:

```
POST /api/methods/design_library
  body: {"genome": "engine-file-id", "preset": "spcas9", "feature_config": "saccer3"}
  returns: {"job_id": "abc123"}

GET /api/jobs/abc123
  returns: {
    "status": "running",           // queued | running | complete | failed | cancelled
    "progress": {
      "pct_complete": 0.45,        // null if unknown
      "rate_per_sec": 20150.0,     // guides/sec, null if unknown
      "eta_seconds": 3             // null if unknown, 0 when complete
    },
    "error": null
  }

GET /api/jobs/abc123/stream        // once status = "complete"
  returns: JSON array of guide regions
```

GenomeHub polls `GET /api/jobs/{id}` until `status` reaches `complete` or
`failed`, then fetches the result via the stream endpoint. GenomeHub owns
file naming — the engine does not set filenames.

## File upload

`POST /api/files/upload` accepts multipart form data with a single field
named `file`. The server:

1. Writes the uploaded bytes to a temp file.
2. Detects format by extension (`.gb`/`.gbk` → GenBank, `.fa`/`.fasta` → FASTA).
3. For GenBank: parses genome, writes a temp FASTA, builds FM-Index from it.
4. For FASTA: builds FM-Index directly.
5. Returns `{"id": "uuid"}`.

The genome, FM-Index, and seed tiers are held in memory. The temp files
can be dropped after indexing completes.

Body size limit: 512 MB (configurable via `DefaultBodyLimit` in `main.rs`).

## Parameter override protocol

For recognition presets, the method also accepts override parameters that
modify the resolved preset before pipeline execution:

| Parameter | Type | Effect |
|-----------|------|--------|
| `pam` | string | Replace preset's PAM pattern (IUPAC) |
| `spacer_len` | string | Replace preset's spacer length |
| `mismatches` | string | Replace preset's mismatch tolerance (0-3) |

The resolution order is: load preset by name → apply overrides → run pipeline.
This allows users to select "SpCas9" but tweak the spacer length without
creating a new YAML file.
