//! Method catalog and generic dispatch.
//!
//! ```text
//! GET  /api/methods           → full catalog (GenomeHub builds UI from this)
//! GET  /api/methods/:id       → single descriptor
//! POST /api/methods/:id       → dispatch; async returns 202 {job_id}
//! ```
//!
//! Adding a new method:
//!   1. Add a descriptor to `METHOD_CATALOG`.
//!   2. Add a match arm to `dispatch_handler`.

use std::sync::Arc;

use axum::extract::{Path, State};
use axum::http::StatusCode;
use axum::Json;
use serde_json::{json, Value};

use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};

use crate::AppState;

// ---------------------------------------------------------------------------
// Catalog — GenomeHub renders its UI entirely from this. No hub-side hardcoding.
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
        ],
        "returns": {
            "type": "file",
            "description": "Scored guide library",
        },
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
///
/// GenomeHub builds its UI entirely from this response.
pub async fn list_methods() -> Json<Value> {
    Json(json!(method_catalog()))
}

/// Return a single method descriptor.
///
/// GenomeHub fetches this before every dispatch to get the current
/// parameter list.
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
///
/// Async methods (async: true) return 202 with `{"job_id": "..."}`.
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

// ---------------------------------------------------------------------------
// Dispatch handlers — one per method
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

    let preset = CRISPRPreset::spcas9();
    let feature_config = FeatureConfig::saccer3();

    let job_id = state.jobs.submit(genome, preset, feature_config);

    Ok((StatusCode::ACCEPTED, Json(json!({ "job_id": job_id }))))
}
