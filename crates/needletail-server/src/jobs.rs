//! Background job management for pipeline execution.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use dashmap::DashMap;
use uuid::Uuid;

/// A Rayon thread pool that deliberately leaves cores free for the Tokio reactor.
///
/// Without this, `par_iter` inside `design_library` saturates every logical
/// core and starves the async executor — the 10Hz HTTP poll from GenomeHub
/// never gets scheduled.
static COMPUTE_POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();

fn compute_pool() -> &'static rayon::ThreadPool {
    COMPUTE_POOL.get_or_init(|| {
        // Leave 2 cores for Tokio (1 minimum if the machine has ≤ 3 cores).
        let total = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4);
        let compute = total.saturating_sub(2).max(1);
        eprintln!(
            "[pool] Rayon compute pool: {} threads ({} total cores, {} reserved for Tokio)",
            compute,
            total,
            total - compute,
        );
        rayon::ThreadPoolBuilder::new()
            .num_threads(compute)
            .thread_name(|i| format!("needletail-compute-{}", i))
            .build()
            .expect("failed to build Rayon compute pool")
    })
}

use needletail_core::io::parquet_file_sink::ParquetFileSink;
use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};
use needletail_core::pipeline::design::{self, LibraryResult, ProgressSink};

use crate::stores::genome_store::StoredGenome;

/// Job status: queued → running → complete | failed | cancelled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JobStatus {
    Queued,
    Running,
    Complete,
    Failed,
    Cancelled,
}

impl JobStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            JobStatus::Queued => "queued",
            JobStatus::Running => "running",
            JobStatus::Complete => "complete",
            JobStatus::Failed => "failed",
            JobStatus::Cancelled => "cancelled",
        }
    }
}

/// Shared progress state for a running job.
///
/// Supports the step/stage/items model for GenomeHub's stepper UI:
/// - `step`: current UI step key (matches method schema `steps[].key`)
/// - `step_stage`: free-text sublabel within the current step
/// - `items_complete`/`items_total`: discrete n-of-x counter
pub struct JobProgress {
    pub stage: std::sync::Mutex<String>,
    pub current: AtomicUsize,
    pub total: AtomicUsize,
    pub cancelled: AtomicBool,
    /// Current UI step key (e.g., "scanning", "scoring").
    pub step: std::sync::Mutex<Option<String>>,
    /// Free-text sublabel within the current step.
    pub step_stage: std::sync::Mutex<Option<String>>,
    /// Item-level progress: complete count.
    pub items_complete: AtomicUsize,
    /// Item-level progress: total count. 0 = no items.
    pub items_total: AtomicUsize,
}

impl ProgressSink for JobProgress {
    fn report(&self, stage: &str, current: usize, total: usize) {
        *self.stage.lock().unwrap() = stage.to_string();
        self.current.store(current, Ordering::Release);
        self.total.store(total, Ordering::Release);
    }

    fn set_step(&self, step: &str) {
        *self.step.lock().unwrap() = Some(step.to_string());
        // Reset stage and items on step change
        *self.step_stage.lock().unwrap() = None;
        self.items_complete.store(0, Ordering::Release);
        self.items_total.store(0, Ordering::Release);
    }

    fn set_stage(&self, stage: &str) {
        *self.step_stage.lock().unwrap() = Some(stage.to_string());
    }

    fn set_items(&self, complete: usize, total: usize) {
        self.items_complete.store(complete, Ordering::Release);
        self.items_total.store(total, Ordering::Release);
    }

    fn clear_items(&self) {
        self.items_complete.store(0, Ordering::Release);
        self.items_total.store(0, Ordering::Release);
    }

    fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Acquire)
    }
}

/// A completed or in-progress job.
pub struct Job {
    pub id: String,
    pub genome_id: String,
    pub status: std::sync::Mutex<JobStatus>,
    pub progress: Arc<JobProgress>,
    pub result: std::sync::Mutex<Option<Result<LibraryResult, String>>>,
    /// Path to the Parquet result file on disk.  The pipeline drains guides
    /// directly to this file via ParquetFileSink — zero RAM accumulation.
    pub result_path: std::sync::Mutex<Option<PathBuf>>,
    pub error: std::sync::Mutex<Option<String>>,
    pub guides_complete: AtomicUsize,
    pub chroms_complete: AtomicUsize,
    pub chroms_total: AtomicUsize,
    /// Monotonic start time (nanos since some epoch, for elapsed calculation).
    pub started_at: AtomicU64,
    pub finished_at: AtomicU64,
}

impl Job {
    /// Elapsed seconds since job started, or None if not started.
    pub fn elapsed_secs(&self) -> Option<f64> {
        let started = self.started_at.load(Ordering::Acquire);
        if started == 0 {
            return None;
        }
        let finished = self.finished_at.load(Ordering::Acquire);
        let end = if finished > 0 { finished } else {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_millis() as u64
        };
        if end >= started {
            Some((end - started) as f64 / 1000.0)
        } else {
            None
        }
    }
}

/// Thread-safe job manager.
pub struct JobManager {
    jobs: DashMap<String, Arc<Job>>,
}

impl JobManager {
    pub fn new() -> Self {
        JobManager {
            jobs: DashMap::new(),
        }
    }

    /// Submit a design_library job to run in the background.
    pub fn submit(
        &self,
        genome: Arc<StoredGenome>,
        preset: CRISPRPreset,
        feature_config: FeatureConfig,
    ) -> String {
        let job_id = Uuid::new_v4().simple().to_string()[..12].to_string();
        let progress = Arc::new(JobProgress {
            stage: std::sync::Mutex::new("queued".into()),
            current: AtomicUsize::new(0),
            total: AtomicUsize::new(0),
            cancelled: AtomicBool::new(false),
            step: std::sync::Mutex::new(None),
            step_stage: std::sync::Mutex::new(None),
            items_complete: AtomicUsize::new(0),
            items_total: AtomicUsize::new(0),
        });

        let chroms_total = genome.genome.chromosomes.len();

        let job = Arc::new(Job {
            id: job_id.clone(),
            genome_id: genome.id.clone(),
            status: std::sync::Mutex::new(JobStatus::Queued),
            progress: progress.clone(),
            result: std::sync::Mutex::new(None),
            result_path: std::sync::Mutex::new(None),
            error: std::sync::Mutex::new(None),
            guides_complete: AtomicUsize::new(0),
            chroms_complete: AtomicUsize::new(0),
            chroms_total: AtomicUsize::new(chroms_total),
            started_at: AtomicU64::new(0),
            finished_at: AtomicU64::new(0),
        });

        self.jobs.insert(job_id.clone(), job.clone());

        // Spawn blocking task
        let jid = job_id.clone();
        tokio::task::spawn_blocking(move || {
            let now_ms = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_millis() as u64;
            job.started_at.store(now_ms, Ordering::Release);
            *job.status.lock().unwrap() = JobStatus::Running;

            // Create results directory and FileSink.
            let results_dir = PathBuf::from("/tmp/needletail-results");
            std::fs::create_dir_all(&results_dir).ok();
            let file_path = results_dir.join(format!("{}.parquet", jid));

            let sink_result = ParquetFileSink::create(&file_path)
                .map_err(|e| format!("failed to create result file: {}", e));

            let result = match sink_result {
                Ok(mut sink) => {
                    // Run the CPU-bound pipeline inside the quarantined Rayon
                    // pool.  This ensures the global Rayon pool (and therefore
                    // the Tokio reactor's threads) are never crowded out.
                    let pipeline_result = compute_pool().install(|| {
                        design::design_library(
                            &genome.genome,
                            &genome.index,
                            genome.tier_small.as_ref(),
                            genome.tier_large.as_ref(),
                            &preset,
                            &feature_config,
                            &*progress,
                            &mut sink,
                        )
                    });

                    // Finish the Parquet file (flush remaining rows + close writer)
                    match pipeline_result {
                        Ok(lib) => {
                            match sink.finish() {
                                Ok(_) => Ok(lib),
                                Err(e) => Err(format!("failed to finalize result file: {}", e)),
                            }
                        }
                        Err(e) => {
                            // Clean up partial file on failure
                            drop(sink);
                            std::fs::remove_file(&file_path).ok();
                            Err(e)
                        }
                    }
                }
                Err(e) => Err(e),
            };

            let end_ms = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_millis() as u64;
            job.finished_at.store(end_ms, Ordering::Release);

            match &result {
                Ok(ref lib) => {
                    job.guides_complete
                        .store(lib.guides_written, Ordering::Release);
                    job.chroms_complete
                        .store(job.chroms_total.load(Ordering::Acquire), Ordering::Release);
                    *job.status.lock().unwrap() = JobStatus::Complete;
                    *job.result_path.lock().unwrap() = Some(file_path);
                }
                Err(e) => {
                    if e == "cancelled" {
                        *job.status.lock().unwrap() = JobStatus::Cancelled;
                    } else {
                        *job.status.lock().unwrap() = JobStatus::Failed;
                        *job.error.lock().unwrap() = Some(e.clone());
                    }
                }
            }

            *job.result.lock().unwrap() = Some(result);
        });

        job_id
    }

    /// Get a job by ID.
    pub fn get(&self, id: &str) -> Option<Arc<Job>> {
        self.jobs.get(id).map(|v| v.clone())
    }

    /// Cancel a job. Best-effort: always returns true if the job exists.
    pub fn cancel(&self, id: &str) -> bool {
        if let Some(job) = self.jobs.get(id) {
            let status = *job.status.lock().unwrap();
            if status == JobStatus::Queued || status == JobStatus::Running {
                job.progress.cancelled.store(true, Ordering::Release);
                *job.status.lock().unwrap() = JobStatus::Cancelled;
            }
            true
        } else {
            false
        }
    }

    /// Remove a job by ID (frees memory for completed/consumed jobs).
    pub fn remove(&self, id: &str) -> bool {
        self.jobs.remove(id).is_some()
    }

    /// List all job IDs with their statuses.
    pub fn list(&self) -> Vec<(String, String, JobStatus)> {
        self.jobs
            .iter()
            .map(|entry| {
                let job = entry.value();
                let status = *job.status.lock().unwrap();
                (job.id.clone(), job.genome_id.clone(), status)
            })
            .collect()
    }
}
