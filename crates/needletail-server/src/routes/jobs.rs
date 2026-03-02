//! Job polling and streaming endpoints.
//!
//! ```text
//! GET    /api/jobs/:id        → status + progress
//! GET    /api/jobs/:id/stream → result stream (JSON)
//! DELETE /api/jobs/:id        → cancel (best-effort); always 200
//! ```

use std::sync::Arc;
use std::sync::atomic::Ordering;

use axum::body::Body;
use axum::extract::{Path, State};
use axum::http::{HeaderValue, StatusCode, header};
use axum::response::Response;
use axum::Json;
use serde_json::{json, Value};
use tokio_util::io::ReaderStream;

use crate::jobs::JobStatus;
use crate::AppState;

/// Poll job status and progress.
pub async fn status(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Json<Value>, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "Job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    let elapsed = job.elapsed_secs();

    let guides_complete = job.guides_complete.load(Ordering::Acquire);

    // Read step/stage/items from progress
    let step = job.progress.step.lock().unwrap().clone();
    let stage = job.progress.step_stage.lock().unwrap().clone();
    let items_complete = job.progress.items_complete.load(Ordering::Acquire);
    let items_total = job.progress.items_total.load(Ordering::Acquire);

    // Items object: only present when items_total > 0
    let items = if items_total > 0 {
        json!({ "complete": items_complete, "total": items_total })
    } else {
        Value::Null
    };

    // pct_complete: per-step, derived from items when available
    let pct: Option<f64> = if items_total > 0 {
        Some((items_complete as f64 / items_total as f64 * 1000.0).round() / 1000.0)
    } else {
        None
    };

    let rate: Option<f64> = match (elapsed, guides_complete) {
        (Some(e), gc) if e > 0.0 && gc > 0 => Some((gc as f64 / e * 10.0).round() / 10.0),
        _ => None,
    };

    let mut resp = json!({
        "status": status.as_str(),
        "step": step,
        "stage": stage,
        "items": items,
        "progress": {
            "pct_complete": pct,
            "rate_per_sec": rate,
        },
        "error": Value::Null,
    });

    if status == JobStatus::Failed {
        if let Some(ref err) = *job.error.lock().unwrap() {
            resp["error"] = json!(err);
        }
    }

    Ok(Json(resp))
}

/// Stream the completed job result.
///
/// Called once by GenomeHub after status reaches "complete".
/// Streams the result file directly from disk — zero RAM.
/// Deletes the file and evicts the job after streaming.
pub async fn stream(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Result<Response, (StatusCode, Json<Value>)> {
    let job = state.jobs.get(&id).ok_or_else(|| {
        (
            StatusCode::NOT_FOUND,
            Json(json!({ "error": "Job not found" })),
        )
    })?;

    let status = *job.status.lock().unwrap();
    if status != JobStatus::Complete {
        return Err((
            StatusCode::ACCEPTED,
            Json(json!({ "error": format!("Job not complete (status: {})", status.as_str()) })),
        ));
    }

    // Take the result path — file on disk written by FileSink.
    let result_path = job.result_path.lock().unwrap().take();
    let taken_result = job.result.lock().unwrap().take();

    match (result_path, taken_result) {
        (Some(path), Some(Ok(_lib_result))) => {
            // Open the file for async streaming
            let file = tokio::fs::File::open(&path).await.map_err(|e| {
                (
                    StatusCode::INTERNAL_SERVER_ERROR,
                    Json(json!({ "error": format!("Failed to open result file: {}", e) })),
                )
            })?;

            let stream = ReaderStream::new(file);
            let body = Body::from_stream(stream);

            // Evict the job from the store.
            drop(job);
            state.jobs.remove(&id);

            // Schedule file cleanup after response is sent.
            // The file stays on disk until the stream completes.
            let cleanup_path = path.clone();
            tokio::spawn(async move {
                // Give the response time to finish streaming
                tokio::time::sleep(std::time::Duration::from_secs(5)).await;
                tokio::fs::remove_file(&cleanup_path).await.ok();
            });

            let mut response = Response::new(body);
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json"),
            );
            Ok(response)
        }
        (_, Some(Err(e))) => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": e })),
        )),
        _ => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": "result already consumed or not available" })),
        )),
    }
}

/// Cancel a running job.
///
/// Best-effort: if the job has already completed or never existed,
/// returns 200 anyway.
pub async fn cancel(
    State(state): State<Arc<AppState>>,
    Path(id): Path<String>,
) -> Json<Value> {
    state.jobs.cancel(&id);
    Json(json!({ "ok": true }))
}
