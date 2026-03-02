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

use needletail_core::io::json::region_to_json;

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
/// Returns streaming JSON array with Content-Disposition header.
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

    // Take the result out of the job — frees the LibraryResult memory
    // after serialisation. Each result holds ~100 MB of guide Regions.
    let taken = job.result.lock().unwrap().take();
    match taken {
        Some(Ok(lib_result)) => {
            // Build streaming JSON array
            let mut chunks: Vec<String> = Vec::with_capacity(lib_result.guides.len() + 2);
            chunks.push("[".to_string());
            for (i, guide) in lib_result.guides.iter().enumerate() {
                let j = region_to_json(guide);
                if i > 0 {
                    chunks.push(",".to_string());
                }
                chunks.push(serde_json::to_string(&j).unwrap_or_default());
            }
            chunks.push("]".to_string());

            // Drop lib_result before building response — frees guides memory
            drop(lib_result);

            // Drop our Arc reference, then evict the job from the store.
            // This frees all job metadata + the now-empty result slot.
            drop(job);
            state.jobs.remove(&id);

            let body = chunks.join("");
            let mut response = Response::new(Body::from(body));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json"),
            );
            Ok(response)
        }
        Some(Err(e)) => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            Json(json!({ "error": e })),
        )),
        None => Err((
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
