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

    // Determine extension for format detection
    let ext = std::path::Path::new(&fname)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

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

    // For GenBank files, we also need a FASTA for FM-index building.
    // Write a temp FASTA derived from the genome.
    let is_genbank = matches!(ext.as_str(), "gb" | "gbk" | "gbff" | "genbank");

    let fasta_path = if is_genbank {
        // Load genome to extract sequences, write temp FASTA
        let genome = needletail_core::io::genbank::load_genbank(&tmp_path).map_err(|e| {
            (
                StatusCode::BAD_REQUEST,
                Json(json!({ "error": format!("Failed to parse GenBank: {}", e) })),
            )
        })?;

        let fa_path = tmp_dir.path().join("genome.fa");
        let mut fa_writer = std::io::BufWriter::new(std::fs::File::create(&fa_path).map_err(|e| {
            (
                StatusCode::INTERNAL_SERVER_ERROR,
                Json(json!({ "error": format!("Failed to create FASTA: {}", e) })),
            )
        })?);

        use std::io::Write;
        for (name, seq) in &genome.sequences {
            writeln!(fa_writer, ">{}", name).unwrap();
            for chunk in seq.chunks(80) {
                fa_writer.write_all(chunk).unwrap();
                writeln!(fa_writer).unwrap();
            }
        }
        drop(fa_writer);

        Some(fa_path.to_string_lossy().to_string())
    } else {
        None
    };

    // Upload into genome store (builds FM-index)
    let fasta_ref = fasta_path.as_deref().or(Some(tmp_str.as_str()));
    let id = state
        .genomes
        .upload(&tmp_str, fasta_ref, None)
        .map_err(|e| {
            (
                StatusCode::BAD_REQUEST,
                Json(json!({ "error": e })),
            )
        })?;

    // Keep temp dir alive by leaking it — the genome store holds Arc refs
    // to the FM-index which was built from the temp files. The index data
    // is in memory so the files aren't needed after upload completes.
    // We can safely let the tempdir drop.

    Ok(Json(json!({ "id": id })))
}
