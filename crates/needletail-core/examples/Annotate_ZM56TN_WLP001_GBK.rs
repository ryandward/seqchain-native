//! Patch ZM56TN_1_WLP001 GenBank file with gene names from the companion GFF.
//!
//! The GBK has AUGUSTUS-predicted gene features with no identifiers.
//! The GFF transcript lines carry `GN=` gene symbols for the same coordinates.
//!
//! This script rewrites the GBK, injecting `/gene=SYMBOL` into every `gene`
//! feature block whose coordinates match a GFF transcript entry.
//!
//! Matching: (LOCUS name, 1-based start, 1-based end) from the location string.
//! Location formats handled: `start..end` and `complement(start..end)`.
//! join() locations are left unmodified — no GFF entry will match them anyway.
//!
//! Output: ZM56TN_1_WLP001_cells_1_reference_named.gbk

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

const GBK_IN:  &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_reference.gbk";
const GFF:     &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_reference.gff";
const GBK_OUT: &str = "/home/ryan.ward/Downloads/ZM56TN_1_WLP001_cells_1_reference_named.gbk";

fn main() {
    // ── 1. Build GFF map: (contig, start1, end1) → gene symbol ───────────────

    let mut gff: HashMap<(String, i64, i64), String> = HashMap::new();

    for line in BufReader::new(File::open(GFF).expect("GFF not found")).lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.is_empty() { continue; }
        let f: Vec<&str> = line.splitn(9, '\t').collect();
        if f.len() < 9 || f[2] != "transcript" { continue; }

        let contig = f[0].to_string();
        let start1: i64 = f[3].parse().unwrap();
        let end1:   i64 = f[4].parse().unwrap();

        let name = extract_attr(f[8], "GN=")
            .or_else(|| extract_attr(f[8], "Name="))
            .unwrap_or_default();

        if !name.is_empty() {
            gff.insert((contig, start1, end1), name);
        }
    }

    eprintln!("[gff] {} named entries", gff.len());

    // ── 2. Rewrite GBK, injecting /gene= into matching gene features ──────────

    let reader = BufReader::new(File::open(GBK_IN).expect("GBK not found"));
    let mut writer = BufWriter::new(File::create(GBK_OUT).expect("cannot create output"));

    let mut current_locus = String::new();

    // State machine: are we inside a `gene` feature block waiting to inject?
    let mut pending_gene: Option<String> = None; // Some(symbol) while in a gene block
    let mut injected = 0usize;
    let mut unmatched = 0usize;

    for line in reader.lines() {
        let line = line.unwrap();

        // Track LOCUS name (contig name in this file)
        if line.starts_with("LOCUS ") || line.starts_with("LOCUS\t") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                current_locus = parts[1].to_string();
            }
            pending_gene = None;
            writeln!(writer, "{}", line).unwrap();
            continue;
        }

        // Detect start of a `gene` feature (5-space indent, then "gene")
        if is_feature_line(&line, "gene") {
            // Extract coords from location string on this line
            let loc = line.trim().trim_start_matches("gene").trim();
            if let Some((s1, e1)) = parse_simple_location(loc) {
                let key = (current_locus.clone(), s1, e1);
                if let Some(symbol) = gff.get(&key) {
                    pending_gene = Some(symbol.clone());
                } else {
                    pending_gene = None;
                    unmatched += 1;
                }
            } else {
                pending_gene = None; // join() or unparseable — skip
            }
            writeln!(writer, "{}", line).unwrap();
            continue;
        }

        // If we're inside a gene block with a pending name, watch for the
        // first qualifier line (/source=) and inject /gene= before it.
        if let Some(ref symbol) = pending_gene.clone() {
            if line.trim_start().starts_with('/') {
                // Inject /gene= immediately before the first qualifier
                writeln!(writer, "                     /gene=\"{}\"", symbol).unwrap();
                injected += 1;
                pending_gene = None;
            }
        }

        // Any line that starts a new feature (non-qualifier, non-continuation)
        // clears the pending state without injection (shouldn't happen, but safe).
        if is_any_feature_line(&line) && !is_feature_line(&line, "gene") {
            pending_gene = None;
        }

        writeln!(writer, "{}", line).unwrap();
    }

    eprintln!("[patch] injected /gene= into {} features, {} unmatched", injected, unmatched);
    eprintln!("[out]   {}", GBK_OUT);
}

/// True if `line` is a feature-table line for the given feature type.
/// Feature lines: exactly 5 leading spaces, then the feature keyword.
fn is_feature_line(line: &str, kind: &str) -> bool {
    if line.len() < 6 { return false; }
    let (indent, rest) = line.split_at(5);
    indent == "     " && rest.starts_with(kind)
        && rest[kind.len()..].starts_with(|c: char| c.is_whitespace())
}

/// True if this line starts any feature (5-space indent + word character).
fn is_any_feature_line(line: &str) -> bool {
    if line.len() < 6 { return false; }
    let (indent, rest) = line.split_at(5);
    indent == "     " && rest.starts_with(|c: char| c.is_alphabetic())
}

/// Parse `start..end` or `complement(start..end)` → (start_1based, end_1based).
/// Returns None for join(), fuzzy coords, or anything else.
fn parse_simple_location(loc: &str) -> Option<(i64, i64)> {
    let inner = if let Some(s) = loc.strip_prefix("complement(").and_then(|s| s.strip_suffix(')')) {
        s
    } else {
        loc
    };

    // Must be exactly "digits..digits"
    let mut parts = inner.splitn(2, "..");
    let s: i64 = parts.next()?.trim_matches(|c: char| !c.is_ascii_digit()).parse().ok()?;
    let e: i64 = parts.next()?.trim_matches(|c: char| !c.is_ascii_digit()).parse().ok()?;
    Some((s, e))
}

/// Extract `KEY=value` from a GFF9 attribute string.
fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    let start = attrs.find(key)?;
    let v_start = start + key.len();
    let v_end = attrs[v_start..].find(';').map(|i| v_start + i).unwrap_or(attrs.len());
    let val = attrs[v_start..v_end].trim().to_string();
    if val.is_empty() { None } else { Some(val) }
}
