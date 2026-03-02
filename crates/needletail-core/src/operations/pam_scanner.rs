//! PAM scanning and guide composition — hot loops over mmap'd genome text.
//!
//! Composes geometry (coordinates), chemistry (IUPAC masks, guide IDs),
//! and engine (ChromGeometry) to scan genomes for CRISPR target sites,
//! validate off-target PAM adjacency, and enrich hits with guide metadata.
//!
//! No Python, no FFI, no disk I/O. Pure computation over byte slices.

use crate::chemistry::{generate_guide_id, CompiledPam, PamDirection, BASE_MASK};
use crate::engine::simd_search::{ChromGeometry, HitAccumulator};
use crate::geometry::{fetch_sequence, interval_envelope, is_low_side, normalize};

// ═══════════════════════════════════════════════════════════════════════════════
//  PAM validation for off-target filtering
// ═══════════════════════════════════════════════════════════════════════════════

/// Check if the PAM adjacent to a spacer hit at genome-global position is valid.
///
/// Uses the same bitmask engine as `scan_region`: fwd_masks for forward-strand
/// hits, rc_masks for reverse-strand hits, applied directly to the genome text.
/// Handles circular chromosome wrapping via `rem_euclid`.
#[inline]
pub fn validate_pam_at(
    text: &[u8],
    chroms: &ChromGeometry,
    sa_pos: usize,
    spacer_len: usize,
    fwd: bool,
    compiled: &CompiledPam,
    direction: PamDirection,
    topologies: Option<&[bool]>,
) -> bool {
    let pam_len = compiled.len;

    // Global → local: binary search on sorted, non-overlapping chrom ranges.
    let ci = match chroms
        .ranges
        .binary_search_by(|&(start, _)| start.cmp(&sa_pos))
    {
        Ok(i) => i,
        Err(0) => return false,
        Err(i) => i - 1,
    };

    let (chrom_start, chrom_len) = chroms.ranges[ci];
    if sa_pos >= chrom_start + chrom_len {
        return false;
    }

    let local_pos = (sa_pos - chrom_start) as i64;
    let cl = chrom_len as i64;
    let circular = topologies.map_or(false, |t| t.get(ci).copied().unwrap_or(false));

    // is_low_side(fwd, downstream) ⇒ spacer on low side ⇒ PAM on high side.
    let pam_start = if is_low_side(fwd, direction.is_downstream()) {
        local_pos + spacer_len as i64
    } else {
        local_pos - pam_len as i64
    };

    // Boundary check for linear chromosomes.
    let pam_end = pam_start + pam_len as i64;
    if !circular && (pam_start < 0 || pam_end > cl) {
        return false;
    }

    // Forward hits → fwd_masks, reverse hits → rc_masks (same as scan_region).
    let masks = if fwd {
        &compiled.fwd_masks
    } else {
        &compiled.rc_masks
    };

    for j in 0..pam_len {
        let genome_idx = chrom_start + ((pam_start + j as i64).rem_euclid(cl)) as usize;
        if BASE_MASK[text[genome_idx] as usize] & masks[j] == 0 {
            return false;
        }
    }

    true
}

/// Filter a HitAccumulator in-place, keeping only hits with valid adjacent PAM.
pub fn filter_hits_by_pam(
    acc: &mut HitAccumulator,
    text: &[u8],
    chroms: &ChromGeometry,
    compiled: &CompiledPam,
    direction: PamDirection,
    spacer_len: usize,
    topologies: Option<&[bool]>,
) {
    let n = acc.query_id.len();
    let mut kept = 0;

    for i in 0..n {
        let sa_pos = acc.position[i] as usize;
        let fwd = acc.strand[i];

        if validate_pam_at(
            text, chroms, sa_pos, spacer_len, fwd, compiled, direction, topologies,
        ) {
            // Unconditional writes — redundant self-copy is cheaper than
            // an unpredictable branch on `kept != i`.
            acc.query_id[kept] = acc.query_id[i];
            acc.position[kept] = acc.position[i];
            acc.strand[kept] = acc.strand[i];
            acc.score[kept] = acc.score[i];
            kept += 1;
        }
    }

    acc.query_id.truncate(kept);
    acc.position.truncate(kept);
    acc.strand.truncate(kept);
    acc.score.truncate(kept);
}

// ═══════════════════════════════════════════════════════════════════════════════
//  SoA output
// ═══════════════════════════════════════════════════════════════════════════════

pub struct PamHits {
    pub chrom_id: Vec<u32>,
    pub position: Vec<u32>,
    pub strand: Vec<i8>,
    pub spacers: Vec<u8>,
    pub pam_seqs: Vec<u8>,
    pub spacer_len: usize,
    pub pam_len: usize,
}

impl PamHits {
    fn with_capacity(cap: usize, spacer_len: usize, pam_len: usize) -> Self {
        PamHits {
            chrom_id: Vec::with_capacity(cap),
            position: Vec::with_capacity(cap),
            strand: Vec::with_capacity(cap),
            spacers: Vec::with_capacity(cap * spacer_len),
            pam_seqs: Vec::with_capacity(cap * pam_len),
            spacer_len,
            pam_len,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Scan kernel
// ═══════════════════════════════════════════════════════════════════════════════

#[inline]
fn scan_region(
    text: &[u8],
    chrom_start: usize,
    chrom_len: usize,
    masks: &[u8],
    strand_val: i8,
    chrom_id: u32,
    spacer_len: usize,
    pam_direction: PamDirection,
    circular: bool,
    hits: &mut PamHits,
) {
    let pat_len = masks.len();
    if chrom_len < pat_len {
        return;
    }
    let scan_end = chrom_start + chrom_len - pat_len + 1;
    let fwd = strand_val > 0;
    let low = is_low_side(fwd, pam_direction.is_downstream());
    let rc = !fwd;

    let mut pos = chrom_start;
    while pos < scan_end {
        let mut ok = true;
        for j in 0..pat_len {
            if BASE_MASK[text[pos + j] as usize] & masks[j] == 0 {
                ok = false;
                break;
            }
        }

        if ok {
            let pam_local = (pos - chrom_start) as i64;
            let pam_end = pam_local + pat_len as i64;

            // sp_start = (start - spacer_len) if is_low_side else end
            let sp_start = if low { pam_local - spacer_len as i64 } else { pam_end };

            if !fetch_sequence(text, chrom_start, chrom_len, sp_start, spacer_len, circular, rc, &mut hits.spacers) {
                pos += 1;
                continue;
            }

            hits.chrom_id.push(chrom_id);
            hits.position.push(pam_local as u32);
            hits.strand.push(strand_val);
            fetch_sequence(text, chrom_start, chrom_len, pam_local, pat_len, circular, rc, &mut hits.pam_seqs);
        }

        pos += 1;
    }
}

fn scan_junction(
    text: &[u8],
    chrom_start: usize,
    chrom_len: usize,
    masks: &[u8],
    strand_val: i8,
    chrom_id: u32,
    spacer_len: usize,
    pam_direction: PamDirection,
    hits: &mut PamHits,
) {
    let pat_len = masks.len();
    if pat_len <= 1 || chrom_len < pat_len {
        return;
    }
    let ov = pat_len - 1;

    let mut junction = Vec::with_capacity(2 * ov);
    junction.extend_from_slice(&text[chrom_start + chrom_len - ov..chrom_start + chrom_len]);
    junction.extend_from_slice(&text[chrom_start..chrom_start + ov]);

    for pos in 0..=junction.len() - pat_len {
        let mut ok = true;
        for j in 0..pat_len {
            if BASE_MASK[junction[pos + j] as usize] & masks[j] == 0 {
                ok = false;
                break;
            }
        }

        if ok && pos < ov && pos + pat_len > ov {
            let chrom_pos = (chrom_len - ov + pos) as i64;
            let fwd = strand_val > 0;
            let pam_end = chrom_pos + pat_len as i64;
            let sp_start = if is_low_side(fwd, pam_direction.is_downstream()) {
                chrom_pos - spacer_len as i64
            } else {
                pam_end
            };

            hits.chrom_id.push(chrom_id);
            hits.position.push((chrom_pos % chrom_len as i64) as u32);
            hits.strand.push(strand_val);

            fetch_sequence(text, chrom_start, chrom_len, sp_start, spacer_len, true, !fwd, &mut hits.spacers);
            fetch_sequence(text, chrom_start, chrom_len, chrom_pos, pat_len, true, !fwd, &mut hits.pam_seqs);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Public entry point
// ═══════════════════════════════════════════════════════════════════════════════

pub fn find_pam_sites(
    text: &[u8],
    chroms: &ChromGeometry,
    pam: &str,
    spacer_len: usize,
    pam_direction: PamDirection,
    topologies: Option<&[bool]>,
) -> Result<PamHits, String> {
    let compiled = CompiledPam::compile(pam)?;
    let total_bases: usize = chroms.ranges.iter().map(|&(_, len)| len).sum();
    let est = total_bases / 8;
    let mut hits = PamHits::with_capacity(est, spacer_len, compiled.len);

    for (ci, &(chrom_start, chrom_len)) in chroms.ranges.iter().enumerate() {
        let cid = ci as u32;
        let circular = topologies.map_or(false, |t| t.get(ci).copied().unwrap_or(false));

        scan_region(text, chrom_start, chrom_len, &compiled.fwd_masks, 1, cid,
                    spacer_len, pam_direction, circular, &mut hits);
        scan_region(text, chrom_start, chrom_len, &compiled.rc_masks, -1, cid,
                    spacer_len, pam_direction, circular, &mut hits);

        if circular {
            scan_junction(text, chrom_start, chrom_len, &compiled.fwd_masks, 1, cid,
                          spacer_len, pam_direction, &mut hits);
            scan_junction(text, chrom_start, chrom_len, &compiled.rc_masks, -1, cid,
                          spacer_len, pam_direction, &mut hits);
        }
    }

    Ok(hits)
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Enriched guide output
// ═══════════════════════════════════════════════════════════════════════════════

pub struct GuideHits {
    pub chrom_ids: Vec<u32>,
    pub pam_positions: Vec<u32>,
    pub guide_starts: Vec<u32>,
    pub guide_ends: Vec<u32>,
    pub strands: Vec<i8>,
    pub spacers: Vec<u8>,
    pub pam_seqs: Vec<u8>,
    pub guide_seqs: Vec<u8>,
    pub guide_ids: Vec<u8>,
    pub spacer_len: usize,
    pub pam_len: usize,
    pub count: usize,
}

pub fn enrich_hits(
    hits: PamHits,
    chrom_names: &[String],
    pam_direction: PamDirection,
    topologies: Option<&[bool]>,
    chrom_ranges: &[(usize, usize)],
) -> GuideHits {
    let n = hits.chrom_id.len();
    let spacer_len = hits.spacer_len;
    let pam_len = hits.pam_len;
    let guide_len = spacer_len + pam_len;
    let pam_first = pam_direction == PamDirection::Upstream;

    let mut out = GuideHits {
        chrom_ids: Vec::with_capacity(n),
        pam_positions: Vec::with_capacity(n),
        guide_starts: Vec::with_capacity(n),
        guide_ends: Vec::with_capacity(n),
        strands: Vec::with_capacity(n),
        spacers: Vec::with_capacity(n * spacer_len),
        pam_seqs: Vec::with_capacity(n * pam_len),
        guide_seqs: Vec::with_capacity(n * guide_len),
        guide_ids: Vec::with_capacity(n * 8),
        spacer_len,
        pam_len,
        count: n,
    };

    for i in 0..n {
        let cid = hits.chrom_id[i];
        let pam_pos = hits.position[i] as i64;
        let strand = hits.strand[i];
        let fwd = strand > 0;
        let circular = topologies.map_or(false, |t| t.get(cid as usize).copied().unwrap_or(false));
        let cl = chrom_ranges[cid as usize].1 as i64;

        // Mirrors extract_spacer exactly:
        //   sp_start = (start - spacer_len) if is_low_side else end
        let pam_end = pam_pos + pam_len as i64;
        let sp_start = if is_low_side(fwd, pam_direction.is_downstream()) { pam_pos - spacer_len as i64 } else { pam_end };
        let sp_end = sp_start + spacer_len as i64;

        // interval_envelope(pam_start, pam_end, sp_start, sp_end)
        let (raw_start, raw_end) = interval_envelope(pam_pos, pam_end, sp_start, sp_end);

        // normalize both endpoints
        let guide_start = normalize(raw_start, cl, circular) as u32;
        let guide_end = normalize(raw_end, cl, circular) as u32;

        let sp = &hits.spacers[i * spacer_len..(i + 1) * spacer_len];
        let pm = &hits.pam_seqs[i * pam_len..(i + 1) * pam_len];

        // guide_seq = (matched + spacer) if pam_first else (spacer + matched)
        if pam_first {
            out.guide_seqs.extend_from_slice(pm);
            out.guide_seqs.extend_from_slice(sp);
        } else {
            out.guide_seqs.extend_from_slice(sp);
            out.guide_seqs.extend_from_slice(pm);
        }

        // generate_guide_id(chrom, guide_start, strand, pam_seq)
        let strand_ch = if fwd { "+" } else { "-" };
        let pam_str = std::str::from_utf8(pm).unwrap_or("");
        let id = generate_guide_id(&chrom_names[cid as usize], guide_start, strand_ch, pam_str);
        out.guide_ids.extend_from_slice(&id);

        out.chrom_ids.push(cid);
        out.pam_positions.push(pam_pos as u32);
        out.guide_starts.push(guide_start);
        out.guide_ends.push(guide_end);
        out.strands.push(strand);
        out.spacers.extend_from_slice(sp);
        out.pam_seqs.extend_from_slice(pm);
    }

    out
}
