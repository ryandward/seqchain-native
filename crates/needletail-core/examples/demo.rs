use std::path::Path;
use std::sync::Arc;
use needletail_core::engine::fm_index::ChromInfo;
use needletail_core::io::genbank::load_genbank;
use needletail_core::io::json::FileSink;
use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};
use needletail_core::pipeline::design::{design_library, NullProgress};
use needletail_core::{FmIndexSearcher, IndexHandle, build_seed_tiers};

fn main() {
    let gb_path = std::env::args().nth(1)
        .unwrap_or_else(|| "test/fixtures/zymomonas.gb".to_string());

    let t0 = std::time::Instant::now();
    let mut genome = load_genbank(Path::new(&gb_path)).unwrap();
    eprintln!("Loaded genome: {} chroms, {} genes ({:.2?})",
        genome.chromosomes.len(), genome.genes().len(), t0.elapsed());

    // Build FM-index directly from the genome's master buffer (zero-copy move)
    let t0 = std::time::Instant::now();
    let chroms: Vec<ChromInfo> = genome.chromosomes.iter()
        .map(|cm| ChromInfo { name: cm.name.clone(), start: cm.start, len: cm.len })
        .collect();
    let text = std::mem::take(&mut genome.text);
    let searcher = FmIndexSearcher::from_text(text, chroms).unwrap();
    let searcher = Arc::new(searcher);
    eprintln!("Built FM-index ({:.2?})", t0.elapsed());

    let t0 = std::time::Instant::now();
    let (ts, tl) = build_seed_tiers(&*searcher, searcher.text(), &gb_path).unwrap();
    eprintln!("Warmed seeds ({:.2?})", t0.elapsed());

    let handle = IndexHandle::Built(searcher);
    let preset = CRISPRPreset::by_name("spcas9").expect("spcas9 preset not found");
    let config = FeatureConfig::by_name("saccer3").expect("saccer3 feature config not found");

    let out_path = Path::new("/tmp/needletail-demo-output.json");
    let t0 = std::time::Instant::now();
    let mut sink = FileSink::create(out_path).expect("failed to create output file");
    let result = design_library(
        &genome, &handle, Some(&ts), Some(&tl),
        &preset, &config, &NullProgress, &mut sink,
    ).unwrap();
    let elapsed = t0.elapsed();
    sink.finish().expect("failed to finalize output file");

    eprintln!("Pipeline: {:.2?} — {} guides written from {} scored → {}",
        elapsed, result.guides_written, result.total_guides_scored,
        out_path.display());
}
