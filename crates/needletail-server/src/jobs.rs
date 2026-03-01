//! Background job management for pipeline execution.

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;


use dashmap::DashMap;
use uuid::Uuid;

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
pub struct JobProgress {
    pub stage: std::sync::Mutex<String>,
    pub current: AtomicUsize,
    pub total: AtomicUsize,
    pub cancelled: AtomicBool,
}

impl ProgressSink for JobProgress {
    fn report(&self, stage: &str, current: usize, total: usize) {
        *self.stage.lock().unwrap() = stage.to_string();
        self.current.store(current, Ordering::Release);
        self.total.store(total, Ordering::Release);
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
        });

        let chroms_total = genome.genome.sequences.len();

        let job = Arc::new(Job {
            id: job_id.clone(),
            genome_id: genome.id.clone(),
            status: std::sync::Mutex::new(JobStatus::Queued),
            progress: progress.clone(),
            result: std::sync::Mutex::new(None),
            error: std::sync::Mutex::new(None),
            guides_complete: AtomicUsize::new(0),
            chroms_complete: AtomicUsize::new(0),
            chroms_total: AtomicUsize::new(chroms_total),
            started_at: AtomicU64::new(0),
            finished_at: AtomicU64::new(0),
        });

        self.jobs.insert(job_id.clone(), job.clone());

        // Spawn blocking task
        tokio::task::spawn_blocking(move || {
            let now_ms = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_millis() as u64;
            job.started_at.store(now_ms, Ordering::Release);
            *job.status.lock().unwrap() = JobStatus::Running;

            let result = design::design_library(
                &genome.genome,
                &genome.index,
                genome.tier_small.as_ref(),
                genome.tier_large.as_ref(),
                &preset,
                &feature_config,
                &*progress,
            );

            let end_ms = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_millis() as u64;
            job.finished_at.store(end_ms, Ordering::Release);

            match &result {
                Ok(ref lib) => {
                    job.guides_complete
                        .store(lib.guides.len(), Ordering::Release);
                    job.chroms_complete
                        .store(job.chroms_total.load(Ordering::Acquire), Ordering::Release);
                    *job.status.lock().unwrap() = JobStatus::Complete;
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
