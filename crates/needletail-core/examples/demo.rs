use std::path::Path;
use std::sync::Arc;
use needletail_core::io::genbank::load_genbank;
use needletail_core::io::json::region_to_json;
use needletail_core::models::preset::{CRISPRPreset, FeatureConfig};
use needletail_core::pipeline::design::{design_library, NullProgress};
use needletail_core::{FmIndexSearcher, IndexHandle, build_seed_tiers};

fn main() {
    let gb_path = std::env::args().nth(1)
        .unwrap_or_else(|| "test/fixtures/zymomonas.gb".to_string());

    let t0 = std::time::Instant::now();
    let genome = load_genbank(Path::new(&gb_path)).unwrap();
    eprintln!("Loaded genome: {} chroms, {} genes ({:.2?})",
        genome.sequences.len(), genome.genes().len(), t0.elapsed());

    // Build FM-index from the genome sequences in memory
    // Write a temp FASTA for the FM-index builder
    let stem = Path::new(&gb_path).file_stem().unwrap().to_string_lossy();
    let tmp_fa_owned = format!("/tmp/{}_needletail.fa", stem);
    let tmp_fa = tmp_fa_owned.as_str();
    {
        use std::io::Write;
        let mut f = std::fs::File::create(tmp_fa).unwrap();
        for (name, seq) in &genome.sequences {
            writeln!(f, ">{}", name).unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
        }
    }

    let t0 = std::time::Instant::now();
    let searcher = FmIndexSearcher::from_fasta(tmp_fa).unwrap();
    let searcher = Arc::new(searcher);
    eprintln!("Built FM-index ({:.2?})", t0.elapsed());

    let t0 = std::time::Instant::now();
    let (ts, tl) = build_seed_tiers(&*searcher, searcher.text(), tmp_fa).unwrap();
    eprintln!("Warmed seeds ({:.2?})", t0.elapsed());

    let handle = IndexHandle::Built(searcher);
    let preset = CRISPRPreset::by_name("spcas9").expect("spcas9 preset not found");
    let config = FeatureConfig::by_name("saccer3").expect("saccer3 feature config not found");

    let t0 = std::time::Instant::now();
    let result = design_library(
        &genome, &handle, Some(&ts), Some(&tl),
        &preset, &config, &NullProgress,
    ).unwrap();
    let elapsed = t0.elapsed();

    eprintln!("Pipeline: {:.2?} — {} promoter guides from {} scored",
        elapsed, result.guides.len(), result.total_guides_scored);

    for guide in result.guides.iter().take(10) {
        let j = region_to_json(guide);
        println!("{}", serde_json::to_string_pretty(&j).unwrap());
    }
}
