//! Needletail HTTP Server — standalone Rust binary serving the CRISPR design API.
//!
//! Replaces SeqChain's Python HTTP layer with Axum.

mod jobs;
mod routes;
mod stores;

use std::sync::Arc;

use axum::extract::DefaultBodyLimit;
use axum::routing::{delete, get, post};
use axum::Router;
use tower_http::cors::CorsLayer;

use crate::jobs::JobManager;
use crate::stores::genome_store::GenomeStore;

/// Shared application state.
pub struct AppState {
    pub genomes: GenomeStore,
    pub jobs: JobManager,
}

#[tokio::main]
async fn main() {
    let state = Arc::new(AppState {
        genomes: GenomeStore::new(),
        jobs: JobManager::new(),
    });

    let app = Router::new()
        // Health
        .route("/api/health", get(routes::health::health))
        // File upload (GenomeHub multipart interface)
        .route("/api/files/upload", post(routes::files::upload))
        // Genomes
        .route("/api/genomes/upload", post(routes::genomes::upload))
        .route("/api/genomes/", get(routes::genomes::list))
        .route("/api/genomes/{id}", get(routes::genomes::get))
        // Methods — catalog + dynamic dispatch
        .route("/api/methods", get(routes::methods::list_methods))
        .route(
            "/api/methods/{method_id}",
            get(routes::methods::get_method).post(routes::methods::dispatch_method),
        )
        // Presets — generic category-based listing
        .route("/api/presets", get(routes::methods::list_preset_categories))
        .route("/api/presets/{category}", get(routes::methods::list_presets))
        // Jobs — status, stream, cancel
        .route("/api/jobs/{id}", get(routes::jobs::status))
        .route("/api/jobs/{id}/stream", get(routes::jobs::stream))
        .route("/api/jobs/{id}", delete(routes::jobs::cancel))
        .layer(DefaultBodyLimit::max(512 * 1024 * 1024)) // 512 MB
        .layer(CorsLayer::permissive())
        .with_state(state);

    let bind = std::env::var("BIND_ADDR").unwrap_or_else(|_| "0.0.0.0:8002".into());
    println!("needletail-server listening on {}", bind);

    let listener = tokio::net::TcpListener::bind(&bind).await.unwrap();
    axum::serve(listener, app).await.unwrap();
}
