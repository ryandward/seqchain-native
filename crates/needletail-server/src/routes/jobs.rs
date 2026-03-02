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
    let chroms_complete = job.chroms_complete.load(Ordering::Acquire);
    let chroms_total = job.chroms_total.load(Ordering::Acquire);

    // Compute progress metrics matching SeqChain's format
    let pct: Option<f64> = if chroms_total > 0 && chroms_complete > 0 {
        Some((chroms_complete as f64 / chroms_total as f64 * 1000.0).round() / 1000.0)
    } else {
        None
    };

    let rate: Option<f64> = match (elapsed, guides_complete) {
        (Some(e), gc) if e > 0.0 && gc > 0 => Some((gc as f64 / e * 10.0).round() / 10.0),
        _ => None,
    };

    let mut eta: Option<i64> = if let Some(e) = elapsed {
        if chroms_complete > 0 {
            let secs_per_chrom = e / chroms_complete as f64;
            Some((secs_per_chrom * (chroms_total - chroms_complete) as f64).round() as i64)
        } else {
            None
        }
    } else {
        None
    };
    if status == JobStatus::Complete {
        eta = Some(0);
    }

    let mut resp = json!({
        "status": status.as_str(),
        "progress": {
            "pct_complete": pct,
            "rate_per_sec": rate,
            "eta_seconds": eta,
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

    let result = job.result.lock().unwrap();
    match result.as_ref() {
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
            Json(json!({ "error": "no result available" })),
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
