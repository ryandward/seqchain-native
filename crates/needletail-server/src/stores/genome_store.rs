//! Thread-safe in-memory genome store.

use std::path::Path;
use std::sync::Arc;

use dashmap::DashMap;
use uuid::Uuid;

use needletail_core::engine::fm_index::ChromInfo;
use needletail_core::io::genbank::{load_fasta, load_genbank};
use needletail_core::models::genome::Genome;
use needletail_core::{FmIndexSearcher, IndexHandle, SeedTier};

/// A stored genome with its FM-Index and seed tiers.
pub struct StoredGenome {
    pub id: String,
    pub genome: Genome,
    pub index: IndexHandle,
    pub tier_small: Option<SeedTier>,
    pub tier_large: Option<SeedTier>,
}

/// Thread-safe genome store.
pub struct GenomeStore {
    genomes: DashMap<String, Arc<StoredGenome>>,
}

impl GenomeStore {
    pub fn new() -> Self {
        GenomeStore {
            genomes: DashMap::new(),
        }
    }

    /// Upload a genome from a file path. Builds FM-Index directly from the
    /// genome's master buffer (zero-copy for GenBank, single-pass for FASTA).
    /// Returns the genome ID.
    pub fn upload(
        &self,
        file_path: &str,
        index_path: Option<&str>,
    ) -> Result<String, String> {
        use std::time::Instant;
        let t0 = Instant::now();

        let path = Path::new(file_path);
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");

        // Load genome
        let is_genbank = matches!(ext, "gb" | "gbk" | "gbff" | "genbank")
            || file_path.ends_with(".gb.gz");
        let is_fasta = matches!(ext, "fa" | "fasta" | "fna" | "fas");
        let mut genome = if is_genbank {
            load_genbank(path)?
        } else if is_fasta {
            load_fasta(path, None)?
        } else {
            return Err(format!("Unrecognized file extension: .{ext}"));
        };
        let t_parse = t0.elapsed();
        eprintln!("[upload] genome parse: {:.3}s", t_parse.as_secs_f64());

        // Build or load FM-Index
        let (index, tier_small, tier_large) = if let Some(idx_path) = index_path {
            let mapped = needletail_core::io::persist::load_index(Path::new(idx_path))
                .map_err(|e| e.to_string())?;
            let mapped = Arc::new(mapped);
            let handle = IndexHandle::Loaded(mapped);

            // Build seed tiers
            let text = handle.text();
            let ts = needletail_core::build_seed_tier_for_handle(&handle, text, idx_path, needletail_core::SEED_K_SMALL);
            let tl = needletail_core::build_seed_tier_for_handle(&handle, text, idx_path, needletail_core::SEED_K_LARGE);

            (handle, ts, tl)
        } else {
            // Build FM-Index directly from genome's master buffer
            let t_idx = Instant::now();
            let chroms: Vec<ChromInfo> = genome.chromosomes.iter()
                .map(|cm| ChromInfo { name: cm.name.clone(), start: cm.start, len: cm.len })
                .collect();
            let text = std::mem::take(&mut genome.text);
            let searcher = FmIndexSearcher::from_text(text, chroms)
                .map_err(|e| e.to_string())?;
            let searcher = Arc::new(searcher);
            eprintln!("[upload] FM-Index build: {:.3}s", t_idx.elapsed().as_secs_f64());

            let t_seed = Instant::now();
            let ts = needletail_core::build_seed_tier_for_handle(
                &IndexHandle::Built(searcher.clone()), searcher.text(), file_path,
                needletail_core::SEED_K_SMALL,
            );
            eprintln!("[upload] seed tier K={}: {:.3}s", needletail_core::SEED_K_SMALL, t_seed.elapsed().as_secs_f64());

            // K=14 tier skipped: only 8% faster scoring but +293 MB RAM.
            // mm=3 falls back to K=10 seeding automatically.
            (IndexHandle::Built(searcher), ts, None)
        };

        eprintln!("[upload] total: {:.3}s", t0.elapsed().as_secs_f64());
        let id = Uuid::new_v4().to_string();
        let stored = Arc::new(StoredGenome {
            id: id.clone(),
            genome,
            index,
            tier_small,
            tier_large,
        });

        self.genomes.insert(id.clone(), stored);
        Ok(id)
    }

    /// Get a genome by ID.
    pub fn get(&self, id: &str) -> Option<Arc<StoredGenome>> {
        self.genomes.get(id).map(|v| v.clone())
    }

    /// List all genome IDs and names.
    pub fn list(&self) -> Vec<(String, String)> {
        self.genomes
            .iter()
            .map(|entry| (entry.id.clone(), entry.genome.name.clone()))
            .collect()
    }

    /// Remove a genome by ID.
    pub fn remove(&self, id: &str) -> bool {
        self.genomes.remove(id).is_some()
    }
}
