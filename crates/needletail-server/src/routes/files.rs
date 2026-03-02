//! File upload endpoint.
//!
//! ```text
//! POST /api/files/upload  → multipart form, field "file" → {"id": "uuid"}
//! ```
//!
//! GenomeHub uploads files here first, then references the returned ID
//! when dispatching methods.

use std::sync::Arc;

use axum::extract::State;
use axum::http::StatusCode;
use axum::Json;
use axum_extra::extract::Multipart;
use serde_json::{json, Value};

use crate::AppState;

/// Accept a multipart file upload, store it, and return a genome ID.
pub async fn upload(
    State(state): State<Arc<AppState>>,
    mut multipart: Multipart,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    // Extract the "file" field
    let mut file_bytes: Option<Vec<u8>> = None;
    let mut file_name: Option<String> = None;

    while let Some(field) = multipart.next_field().await.map_err(|e| {
        (
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": format!("Multipart error: {}", e) })),
        )
    })? {
        let name = field.name().unwrap_or("").to_string();
        if name == "file" {
            file_name = field.file_name().map(|s| s.to_string());
            file_bytes = Some(field.bytes().await.map_err(|e| {
                (
                    StatusCode::BAD_REQUEST,
                    Json(json!({ "error": format!("Failed to read file: {}", e) })),
                )
            })?.to_vec());
        }
    }

    let bytes = file_bytes.ok_or_else(|| {
        (
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": "Missing 'file' field in multipart form" })),
        )
    })?;

    let fname = file_name.unwrap_or_else(|| "upload.gb".into());

    // Write to a temp file so the genome store can process it
    let tmp_dir = tempfile::tempdir().map_err(|e| {
        (
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": format!("Failed to create temp dir: {}", e) })),
        )
    })?;

    let tmp_path = tmp_dir.path().join(&fname);
    std::fs::write(&tmp_path, &bytes).map_err(|e| {
        (
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": format!("Failed to write temp file: {}", e) })),
        )
    })?;

    let tmp_str = tmp_path.to_string_lossy().to_string();

    // Upload into genome store (parses genome + builds FM-index directly)
    let id = state
        .genomes
        .upload(&tmp_str, None)
        .map_err(|e| {
            (
                StatusCode::BAD_REQUEST,
                Json(json!({ "error": e })),
            )
        })?;

    Ok(Json(json!({ "id": id })))
}
