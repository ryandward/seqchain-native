//! Thread-safe in-memory genome store.

use std::path::Path;
use std::sync::Arc;

use dashmap::DashMap;
use uuid::Uuid;

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

    /// Upload a genome from a file path. Builds FM-Index from the genome's FASTA.
    /// Returns the genome ID.
    pub fn upload(
        &self,
        file_path: &str,
        fasta_path: Option<&str>,
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
        let genome = if ext == "gb" || file_path.ends_with(".gb.gz") {
            load_genbank(path)?
        } else {
            load_fasta(path, None)?
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
        } else if let Some(fa_path) = fasta_path {
            let t_idx = Instant::now();
            let searcher = FmIndexSearcher::from_fasta(fa_path)
                .map_err(|e| e.to_string())?;
            let searcher = Arc::new(searcher);
            eprintln!("[upload] FM-Index build: {:.3}s", t_idx.elapsed().as_secs_f64());

            let t_seed = Instant::now();
            let (ts, tl) = needletail_core::build_seed_tiers(&*searcher, searcher.text(), fa_path)
                .map_err(|e| e.to_string())?;
            eprintln!("[upload] seed tiers: {:.3}s", t_seed.elapsed().as_secs_f64());

            (IndexHandle::Built(searcher), Some(ts), Some(tl))
        } else {
            return Err("Either fasta_path or index_path must be provided".into());
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
