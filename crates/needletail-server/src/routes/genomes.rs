use std::sync::Arc;

use axum::extract::{Path, State};
use axum::http::StatusCode;
use axum::Json;
use serde::Deserialize;
use serde_json::{json, Value};

use crate::AppState;

#[derive(Deserialize)]
pub struct UploadRequest {
    pub file_path: String,
    pub index_path: Option<String>,
}

pub async fn upload(
    State(state): State<Arc<AppState>>,
    Json(req): Json<UploadRequest>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let index = req.index_path.as_deref();

    match state.genomes.upload(&req.file_path, index) {
        Ok(id) => Ok(Json(json!({ "id": id }))),
        Err(e) => Err((
            StatusCode::BAD_REQUEST,
            Json(json!({ "error": e })),
        )),
    }
}

pub async fn list(State(state): State<Arc<AppState>>) -> Json<Value> {
    let genomes: Vec<Value> = state
        .genomes
        .list()
        .into_iter()
        .map(|(id, name)| json!({ "id": id, "name": name }))
        .collect();
    Json(json!(genomes))
}

pub async fn get(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    match state.genomes.get(&id) {
        Some(stored) => {
            let chroms: Vec<Value> = stored
                .genome
                .chrom_lengths()
                .into_iter()
                .map(|(name, len)| json!({ "name": name, "length": len }))
                .collect();
            Ok(Json(json!({
                "id": stored.id,
                "name": stored.genome.name,
                "chromosomes": chroms,
                "n_features": stored.genome.features.len(),
                "n_genes": stored.genome.genes().len(),
            })))
        }
        None => Err((
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "genome not found" })),
        )),
    }
}
