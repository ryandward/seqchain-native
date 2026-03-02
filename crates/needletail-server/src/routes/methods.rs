//! Method catalog and generic dispatch.
//!
//! ```text
//! GET  /api/methods           → full catalog (GenomeHub builds UI from this)
//! GET  /api/methods/:id       → single descriptor
//! POST /api/methods/:id       → dispatch; async returns 202 {job_id}
//! GET  /api/presets/crispr    → list CRISPR nuclease presets
//! GET  /api/presets/features  → list feature config presets
//! ```
//!
//! Adding a new method:
//!   1. Add a descriptor to `method_catalog()`.
//!   2. Add a match arm to `dispatch_method`.

use std::sync::Arc;

use axum::extract::{Path, State};
use axum::http::StatusCode;
use axum::Json;
use serde_json::{json, Value};

use needletail_core::models::preset::{CRISPRPreset, FeatureConfig, PresetRegistry};

use crate::AppState;

// ---------------------------------------------------------------------------
// Catalog — programmatically generated from the preset registry.
// GenomeHub renders its UI entirely from this. No hub-side hardcoding.
// ---------------------------------------------------------------------------

fn method_catalog() -> Vec<Value> {
    vec![json!({
        "id": "design_library",
        "name": "Library Design",
        "description": "Design and score a CRISPR guide library targeting gene promoters.",
        "async": true,
        "parameters": [
            {
                "name": "genome",
                "type": "file",
                "required": true,
                "description": "Target genome",
                "accept": ["gb", "gbk", "gbff", "genbank", "fa", "fasta", "fna"],
            },
            {
                "name": "preset",
                "type": "select",
                "required": false,
                "description": "Recognition preset",
                "default": PresetRegistry::list("recognition").first().copied().unwrap_or("spcas9"),
                "options": PresetRegistry::select_options("recognition"),
            },
            {
                "name": "feature_config",
                "type": "select",
                "required": false,
                "description": "Tiling preset",
                "default": PresetRegistry::list("tiling").first().copied().unwrap_or("saccer3"),
                "options": PresetRegistry::select_options("tiling"),
            },
            {
                "name": "pam",
                "type": "string",
                "required": false,
                "description": "Override PAM pattern (IUPAC, e.g. NGG, TTTN)",
            },
            {
                "name": "spacer_len",
                "type": "string",
                "required": false,
                "description": "Override spacer length in bp",
            },
            {
                "name": "mismatches",
                "type": "string",
                "required": false,
                "description": "Max mismatches for off-target scoring (0-3)",
                "default": "0",
            },
        ],
        "returns": {
            "type": "file",
            "description": "Scored guide library",
        },
        "steps": [
            { "key": "indexing",   "label": "Building FM-Index" },
            { "key": "scanning",   "label": "Scanning genome" },
            { "key": "scoring",    "label": "Scoring & annotation" },
            { "key": "filtering",  "label": "Filtering library" },
        ],
    })]
}

fn method_index() -> std::collections::HashMap<String, Value> {
    method_catalog()
        .into_iter()
        .filter_map(|m| {
            let id = m.get("id")?.as_str()?.to_string();
            Some((id, m))
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Routes
// ---------------------------------------------------------------------------

/// Return the full method catalog.
pub async fn list_methods() -> Json<Value> {
    Json(json!(method_catalog()))
}

/// Return a single method descriptor.
pub async fn get_method(
    Path(method_id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let index = method_index();
    match index.get(&method_id) {
        Some(m) => Ok(Json(m.clone())),
        None => Err((
            StatusCode::NOT_FOUND,
            Json(json!({ "error": format!("Method '{}' not found", method_id) })),
        )),
    }
}

/// Dispatch a method.
pub async fn dispatch_method(
    State(state): State<Arc<AppState>>,
    Path(method_id): Path<String>,
    body: Option<Json<Value>>,
) -> Result<(StatusCode, Json<Value>), (StatusCode, Json<Value>)> {
    let index = method_index();
    if !index.contains_key(&method_id) {
        return Err((
            StatusCode::NOT_FOUND,
            Json(json!({ "error": format!("Method '{}' not found", method_id) })),
        ));
    }

    let body = body.map(|Json(v)| v).unwrap_or(json!({}));

    match method_id.as_str() {
        "design_library" => dispatch_design_library(&state, &body).await,
        _ => Err((
            StatusCode::NOT_IMPLEMENTED,
            Json(json!({ "error": format!("No dispatcher registered for '{}'", method_id) })),
        )),
    }
}

/// List all preset categories.
pub async fn list_preset_categories() -> Json<Value> {
    Json(json!(PresetRegistry::categories()))
}

/// List presets in a category with full detail.
pub async fn list_presets(
    Path(category): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let presets = PresetRegistry::detail(&category);
    if presets.is_empty() && !PresetRegistry::categories().contains(&category.as_str()) {
        return Err((
            StatusCode::NOT_FOUND,
            Json(json!({ "error": format!("Unknown preset category '{}'", category) })),
        ));
    }
    Ok(Json(json!(presets)))
}

// ---------------------------------------------------------------------------
// Dispatch handlers
// ---------------------------------------------------------------------------

async fn dispatch_design_library(
    state: &AppState,
    body: &Value,
) -> Result<(StatusCode, Json<Value>), (StatusCode, Json<Value>)> {
    let genome_id = body
        .get("genome")
        .and_then(|v| v.as_str())
        .ok_or_else(|| {
            (
                StatusCode::UNPROCESSABLE_ENTITY,
                Json(json!({ "error": "Missing required parameter: genome" })),
            )
        })?;

    let genome = state.genomes.get(genome_id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": format!("Genome '{}' not found", genome_id) })),
        )
    })?;

    // Resolve preset (default: spcas9)
    let preset_name = body
        .get("preset")
        .and_then(|v| v.as_str())
        .unwrap_or("spcas9");
    let mut preset = CRISPRPreset::by_name(preset_name).ok_or_else(|| {
        let available = CRISPRPreset::list().join(", ");
        (
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": format!("Unknown preset '{}'. Available: {}", preset_name, available) })),
        )
    })?;

    // Apply overrides
    if let Some(pam) = body.get("pam").and_then(|v| v.as_str()) {
        preset.pam = pam.to_string();
    }
    if let Some(sl) = body.get("spacer_len").and_then(|v| v.as_u64()) {
        preset.spacer_len = sl as usize;
    }
    if let Some(mm) = body.get("mismatches").and_then(|v| v.as_u64()) {
        preset.mismatches = mm as u8;
    }

    // Resolve feature config (default: saccer3)
    let config_name = body
        .get("feature_config")
        .and_then(|v| v.as_str())
        .unwrap_or("saccer3");
    let feature_config = FeatureConfig::by_name(config_name).ok_or_else(|| {
        let available = FeatureConfig::list().join(", ");
        (
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": format!("Unknown feature config '{}'. Available: {}", config_name, available) })),
        )
    })?;

    let job_id = state.jobs.submit(genome, preset, feature_config);

    Ok((StatusCode::ACCEPTED, Json(json!({ "job_id": job_id }))))
}
