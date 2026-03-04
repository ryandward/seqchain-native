//! Tensor pivot + splice-aware affine gap extension.
//!
//! Two primitives, one physical law:
//!
//! **`pivot_reads`** — The Tensor Pivot.
//! Transforms K row-major reads (`K × L` bytes) into L column-major SIMD
//! lanes (`L × [i32; 8]`). The spatial axis of the sequence becomes the
//! temporal axis of SIMD execution. The O(K) gather latency is paid exactly
//! once; the extension loop consumes contiguous array loads.
//!
//! **`extend_batch`** — The Wind.
//! Extends up to 8 BWT seed anchors simultaneously. Each SIMD lane carries
//! one (read, reference_position) pair through the affine splice automaton.
//! All conditional state transitions are collapsed via bit-resonance masking:
//! `blendv`, `cmpeq`, `max`. No branches enter the inner loop.
//!
//! # State tensor
//!
//! At reference step `t` (advancing all 8 reference positions by 1):
//!
//! ```text
//! S_t = (M_t, E_t, J_t) ∈ ℤ^{8×3}
//!
//! M_t  — best score at main diagonal: read[k][t] vs ref[anchor_k + t]
//! E_t  — best score in ref-gap state (ref advances; read held)
//! J_t  — best score in intron state  (large ref-gap, splice-gated)
//! ```
//!
//! # Recurrence (branchless)
//!
//! ```text
//! σ_t   = blendv(-mismatch, +match, cmpeq(read_col, ref_bases))
//! M_t+1 = max(M_t + σ, E_t + σ, (J_t - γ_ic + σ) ⊗ M_AG, 0)
//! E_t+1 = max(M_t+1 - γ_o, E_t - γ_e)
//! J_t+1 = blendv(J_t - γ_ie, M_t+1 - γ_io, M_GT) ⊗ ¬M_AG
//! ```
//!
//! where `⊗` = bitwise AND (Hadamard on Z₂), `M_GT`/`M_AG` are the
//! characteristic masks of GT-donor / AG-acceptor splice signals.
//!
//! # Import hierarchy
//! No biology. No Python. Depends only on `std`.

// ── GT = b'G' | (b'T' << 8) in little-endian memory layout
const GT_LE16: i32 = (b'T' as i32) << 8 | b'G' as i32; // 0x5447
// ── AG = b'A' | (b'G' << 8) in little-endian memory layout
const AG_LE16: i32 = (b'G' as i32) << 8 | b'A' as i32; // 0x4741

// ═══════════════════════════════════════════════════════════════════════════
//  SpliceParams — scalar penalty constants, platform-agnostic
// ═══════════════════════════════════════════════════════════════════════════

/// Affine gap + splice-junction scoring parameters.
///
/// Constructed once per pipeline; broadcast into SIMD registers by
/// `SpliceParamVec::from(params)`.
#[derive(Debug, Clone)]
pub struct SpliceParams {
    pub match_score:   i32,  // reward for a matching base          (default  2)
    pub mismatch_pen:  i32,  // penalty for a substitution          (default  3)
    pub gap_open:      i32,  // cost to open a small gap            (default  6)
    pub gap_ext:       i32,  // cost to extend a small gap per base (default  1)
    pub intron_open:   i32,  // cost to open an intron jump         (default 20)
    pub intron_ext:    i32,  // cost to extend intron per base      (default  0)
    pub intron_close:  i32,  // cost to close an intron at AG       (default  5)
}

impl Default for SpliceParams {
    fn default() -> Self {
        Self {
            match_score:  2,
            mismatch_pen: 3,
            gap_open:     6,
            gap_ext:      1,
            intron_open:  20,
            intron_ext:   0,
            intron_close: 5,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tensor Pivot — public, platform-agnostic
// ═══════════════════════════════════════════════════════════════════════════

/// The Tensor Pivot: K row-major reads → L column-major SIMD lanes.
///
/// `pivot_reads(reads, len)` returns a `Vec` of length `len`, where each
/// element is `[i32; 8]` and `result[t][k] = reads[k][t] as i32`.
///
/// Reads with index ≥ `reads.len()` are padded with `0x00` (guaranteed
/// mismatch against any real nucleotide; their lanes score to zero).
///
/// On x86-64 with AVX2, the extension kernel will load these directly as
/// `_mm256_loadu_si256`, costing one L1 load per column.
///
/// # Panics
/// If `reads.len() > 8` or any `reads[k].len() < len`.
pub fn pivot_reads(reads: &[&[u8]], len: usize) -> Vec<[i32; 8]> {
    assert!(reads.len() <= 8, "pivot_reads: at most 8 reads per SIMD batch");
    debug_assert!(
        reads.iter().all(|r| r.len() >= len),
        "pivot_reads: all reads must have length >= len"
    );

    let k = reads.len();
    (0..len)
        .map(|t| {
            let mut col = [0i32; 8];
            for lane in 0..k {
                col[lane] = reads[lane][t] as i32;
            }
            col
        })
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════
//  extend_batch — public entry point, dispatches to AVX2 or scalar
// ═══════════════════════════════════════════════════════════════════════════

/// Align up to 8 reads to their BWT seed anchors simultaneously.
///
/// Returns one `(best_score, best_ref_endpoint)` per lane. Lanes with
/// no real read (padding) will have score ≤ 0 and should be discarded.
///
/// # Arguments
/// * `genome`  — reference text (must have ≥ `max(anchors) + max_ext + 4` bytes)
/// * `reads`   — up to 8 read sequences, all `>= max_ext` bytes long
/// * `anchors` — BWT-identified reference positions, one per read
/// * `params`  — scoring parameters
/// * `max_ext` — maximum extension length (capped at min read length)
pub fn extend_batch(
    genome:  &[u8],
    reads:   &[&[u8]],
    anchors: &[u32],
    params:  &SpliceParams,
    max_ext: usize,
) -> Vec<(i32, u32)> {
    assert_eq!(reads.len(), anchors.len());
    assert!(reads.len() <= 8);

    let n = reads.len();
    if n == 0 {
        return Vec::new();
    }

    let ext = max_ext.min(reads.iter().map(|r| r.len()).min().unwrap_or(0));
    if ext == 0 {
        return vec![(0, 0); n];
    }

    // Pivot: pays the O(K × L) gather cost once before the inner loop.
    let cols = pivot_reads(reads, ext);

    #[cfg(target_arch = "x86_64")]
    if is_x86_feature_detected!("avx2") {
        // Safety: feature check above guarantees AVX2 availability.
        let (scores, ends) = unsafe { extend_avx2(genome, &cols, anchors, params, ext) };
        return (0..n).map(|k| (scores[k], ends[k])).collect();
    }

    extend_scalar(genome, &cols, anchors, params, ext)
}

// ═══════════════════════════════════════════════════════════════════════════
//  AVX2 kernel
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(target_arch = "x86_64")]
mod avx2_impl {
    use std::arch::x86_64::*;
    use super::{GT_LE16, AG_LE16, SpliceParams};

    // ── Broadcast constants ────────────────────────────────────────────────
    pub(super) struct ParamVec {
        pub k_match:       __m256i,
        pub k_mismatch:    __m256i,
        pub k_gap_open:    __m256i,
        pub k_gap_ext:     __m256i,
        pub k_intron_open: __m256i,
        pub k_intron_ext:  __m256i,
        pub k_intron_cls:  __m256i,
        pub k_neg_inf:     __m256i,
        pub k_zero:        __m256i,
        pub k_byte_mask:   __m256i,
    }

    impl ParamVec {
        #[target_feature(enable = "avx2")]
        pub(super) unsafe fn new(p: &SpliceParams) -> Self {
            // NEG_INF: far below any real score, won't overflow with int32.
            const NEG_INF: i32 = -1_000_000;
            Self {
                k_match:       _mm256_set1_epi32(p.match_score),
                k_mismatch:    _mm256_set1_epi32(-p.mismatch_pen),
                k_gap_open:    _mm256_set1_epi32(p.gap_open),
                k_gap_ext:     _mm256_set1_epi32(p.gap_ext),
                k_intron_open: _mm256_set1_epi32(p.intron_open),
                k_intron_ext:  _mm256_set1_epi32(p.intron_ext),
                k_intron_cls:  _mm256_set1_epi32(p.intron_close),
                k_neg_inf:     _mm256_set1_epi32(NEG_INF),
                k_zero:        _mm256_setzero_si256(),
                k_byte_mask:   _mm256_set1_epi32(0xFF),
            }
        }
    }

    // ── SIMD state ─────────────────────────────────────────────────────────
    pub(super) struct State {
        pub m:   __m256i,  // 8 × i32 match scores
        pub gap: __m256i,  // 8 × i32 ref-gap scores
        pub jmp: __m256i,  // 8 × i32 intron-jump scores
    }

    /// Detect a dinucleotide pattern at 8 arbitrary genome positions.
    ///
    /// `_mm256_i32gather_epi32(base, offsets, 1)` reads 4 bytes per lane
    /// (byte-addressed, hardware handles unaligned access). We AND the low
    /// 16 bits to get the first two bytes, then compare to `pattern_le16`.
    ///
    /// Returns `0xFFFFFFFF` per lane on match, `0x00000000` otherwise.
    #[target_feature(enable = "avx2")]
    #[inline]
    pub(super) unsafe fn dinucleotide_mask(
        genome:       *const u8,
        positions:    __m256i,
        pattern_le16: i32,
    ) -> __m256i {
        let words   = _mm256_i32gather_epi32(genome as *const i32, positions, 1);
        let low16   = _mm256_and_si256(words, _mm256_set1_epi32(0x0000_FFFF));
        _mm256_cmpeq_epi32(low16, _mm256_set1_epi32(pattern_le16))
    }

    /// One column of the affine + splice-aware DP.
    ///
    /// Advances all 8 reference positions by 1. All transitions are branchless:
    ///
    /// ```text
    /// σ       = blendv(-mismatch, +match, cmpeq(read_col, ref_bases))
    /// M_new   = max(M + σ, E + σ, (J - γ_ic + σ) ⊗ M_AG, 0)
    /// E_new   = max(M_new - γ_o, E - γ_e)          ← uses M_new: gap opened NOW
    /// J_new   = blendv(J - γ_ie, M_new - γ_io, M_GT) ⊗ ¬M_AG
    /// ```
    ///
    /// The E and J transitions use `M_new` rather than `M_old`; this is the
    /// standard "column-after-update" approximation used in BLAST and SIMD
    /// Smith-Waterman implementations. It accurately models the common case
    /// (gap opened at the current position) and slightly overestimates for
    /// off-diagonal paths — acceptable for seed verification.
    #[target_feature(enable = "avx2")]
    #[inline]
    pub(super) unsafe fn step(
        s:            &State,
        read_col:     __m256i,      // pivot column: 8 × i32 read bases
        ref_positions:__m256i,      // 8 × i32 reference byte offsets
        genome:       *const u8,
        p:            &ParamVec,
    ) -> State {
        // ── Splice signal detection ──────────────────────────────────────
        let is_donor    = dinucleotide_mask(genome, ref_positions, GT_LE16);
        let accept_pos  = _mm256_sub_epi32(ref_positions, _mm256_set1_epi32(1));
        let is_acceptor = dinucleotide_mask(genome, accept_pos, AG_LE16);

        // ── σ_t: substitution score (single blendv, zero branches) ──────
        let ref_words = _mm256_i32gather_epi32(genome as *const i32, ref_positions, 1);
        let ref_bases = _mm256_and_si256(ref_words, p.k_byte_mask);
        let is_match  = _mm256_cmpeq_epi32(read_col, ref_bases);
        let sigma     = _mm256_blendv_epi8(p.k_mismatch, p.k_match, is_match);

        // ── M_new: four sources, collapsed by max ────────────────────────
        // Source 1 — diagonal (M + σ): both read and ref consumed.
        let from_m = _mm256_add_epi32(s.m, sigma);
        // Source 2 — close ref-gap (E + σ): gap closed, current base matched.
        let from_e = _mm256_add_epi32(s.gap, sigma);
        // Source 3 — close intron (J - γ_ic + σ): only where AG resonates.
        let j_closed  = _mm256_add_epi32(
            _mm256_sub_epi32(s.jmp, p.k_intron_cls),
            sigma,
        );
        let from_j = _mm256_blendv_epi8(p.k_neg_inf, j_closed, is_acceptor);

        // ⊕ = OR-of-maxima; local floor at 0 (Smith-Waterman style).
        let m_new = _mm256_max_epi32(
            _mm256_max_epi32(from_m, from_e),
            _mm256_max_epi32(from_j, p.k_zero),
        );

        // ── E_new: open or extend ref-gap (uses M_new) ──────────────────
        let e_open   = _mm256_sub_epi32(m_new,  p.k_gap_open);
        let e_extend = _mm256_sub_epi32(s.gap, p.k_gap_ext);
        let gap_new  = _mm256_max_epi32(e_open, e_extend);

        // ── J_new: open on GT, extend, silence on AG ─────────────────────
        let j_open_cand   = _mm256_sub_epi32(m_new,  p.k_intron_open);
        let j_extend_cand = _mm256_sub_epi32(s.jmp,  p.k_intron_ext);
        // blendv: M_GT gates opening; non-donor lanes extend.
        let jmp_raw = _mm256_blendv_epi8(j_extend_cand, j_open_cand, is_donor);
        // Silence J on acceptor: score was consumed by from_j → m_new.
        let jmp_new = _mm256_blendv_epi8(jmp_raw, p.k_neg_inf, is_acceptor);

        State { m: m_new, gap: gap_new, jmp: jmp_new }
    }
} // mod avx2_impl

// ── AVX2 top-level driver ──────────────────────────────────────────────────

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn extend_avx2(
    genome:  &[u8],
    cols:    &[[i32; 8]],
    anchors: &[u32],
    params:  &SpliceParams,
    max_ext: usize,
) -> ([i32; 8], [u32; 8]) {
    use std::arch::x86_64::*;
    use avx2_impl::{ParamVec, State, step};

    let p = ParamVec::new(params);

    // Pack anchor positions into one register.
    // Panics on non-8 batches are caught by the caller (extend_batch).
    let mut anchor_arr = [0i32; 8];
    for (k, &a) in anchors.iter().enumerate() {
        anchor_arr[k] = a as i32;
    }
    let mut ref_pos = _mm256_loadu_si256(anchor_arr.as_ptr() as *const __m256i);
    let one         = _mm256_set1_epi32(1);

    // Initial state: all M/E/J at NEG_INF except M = 0 (local alignment seed).
    let mut state = State {
        m:   p.k_zero,
        gap: p.k_neg_inf,
        jmp: p.k_neg_inf,
    };

    let mut best_m   = p.k_zero;
    let mut best_pos = ref_pos;

    let genome_ptr = genome.as_ptr();

    for t in 0..max_ext {
        // The pivot column: one contiguous 32-byte load. Zero gather stalls.
        let read_col = _mm256_loadu_si256(cols[t].as_ptr() as *const __m256i);

        state = step(&state, read_col, ref_pos, genome_ptr, &p);

        // Track per-lane best score and its reference endpoint.
        // _mm256_cmpgt_epi32: 0xFFFFFFFF where state.m > best_m, else 0.
        let improved = _mm256_cmpgt_epi32(state.m, best_m);
        best_m   = _mm256_blendv_epi8(best_m,   state.m,  improved);
        best_pos = _mm256_blendv_epi8(best_pos,  ref_pos,  improved);

        // Advance all 8 reference positions simultaneously: O(1) instruction.
        ref_pos = _mm256_add_epi32(ref_pos, one);
    }

    let mut scores: [i32; 8] = [0; 8];
    let mut ends:   [u32; 8] = [0; 8];
    _mm256_storeu_si256(scores.as_mut_ptr() as *mut __m256i, best_m);
    _mm256_storeu_si256(ends.as_mut_ptr()   as *mut __m256i, best_pos);

    (scores, ends)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Scalar fallback (all platforms, no gaps — mismatch-only diagonal score)
// ═══════════════════════════════════════════════════════════════════════════

/// Scalar diagonal extension: processes lanes sequentially, no gap states.
///
/// Each lane scores `read[k][t]` vs `genome[anchor_k + t]` with `+match`
/// or `-mismatch`, floored at 0 (local alignment). No intron detection.
///
/// Used on non-AVX2 hardware and as a correctness reference.
fn extend_scalar(
    genome:  &[u8],
    cols:    &[[i32; 8]],
    anchors: &[u32],
    params:  &SpliceParams,
    max_ext: usize,
) -> Vec<(i32, u32)> {
    let n = anchors.len();
    let mut best_score = vec![0i32; n];
    let mut best_end   = anchors.to_vec();
    let mut score      = vec![0i32; n];

    for t in 0..max_ext {
        let col = &cols[t];
        for k in 0..n {
            let ref_idx = anchors[k] as usize + t;
            if ref_idx >= genome.len() {
                break;
            }
            let ref_base = genome[ref_idx] as i32;
            let read_base = col[k];
            let s = if read_base == ref_base {
                params.match_score
            } else {
                -params.mismatch_pen
            };
            score[k] = (score[k] + s).max(0);
            if score[k] > best_score[k] {
                best_score[k] = score[k];
                best_end[k]   = ref_idx as u32;
            }
        }
    }

    (0..n).map(|k| (best_score[k], best_end[k])).collect()
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pivot_layout() {
        // pivot_reads[t][k] = reads[k][t]
        let r0: &[u8] = b"ACGT";
        let r1: &[u8] = b"TTTT";
        let cols = pivot_reads(&[r0, r1], 4);
        assert_eq!(cols[0][0], b'A' as i32);
        assert_eq!(cols[0][1], b'T' as i32);
        assert_eq!(cols[1][0], b'C' as i32);
        assert_eq!(cols[2][0], b'G' as i32);
        assert_eq!(cols[3][0], b'T' as i32);
        // Padded lanes are zero.
        assert_eq!(cols[0][2], 0);
    }

    #[test]
    fn pivot_padding() {
        let r: &[u8] = b"A";
        let cols = pivot_reads(&[r], 1);
        // Lanes 1-7 are 0.
        for lane in 1..8 {
            assert_eq!(cols[0][lane], 0);
        }
    }

    #[test]
    fn extend_batch_exact_match() {
        // Read identical to genome → max score = match_score × read_len.
        let genome: Vec<u8> = b"ACGTACGTACGT".to_vec();
        let read: &[u8] = b"ACGT";
        let params = SpliceParams::default();
        let results = extend_batch(&genome, &[read], &[0], &params, 4);
        assert_eq!(results.len(), 1);
        let (score, _end) = results[0];
        // Score = 4 × match_score = 8.
        assert_eq!(score, 4 * params.match_score);
    }

    #[test]
    fn extend_batch_all_mismatch() {
        let genome: Vec<u8> = b"TTTTTTTT".to_vec();
        let read: &[u8] = b"AAAAAAAA";
        let params = SpliceParams::default();
        let results = extend_batch(&genome, &[read], &[0], &params, 8);
        // Local alignment: floor at 0.
        assert_eq!(results[0].0, 0);
    }

    #[test]
    fn gt_ag_constants() {
        // Verify the GT/AG little-endian constants are correct.
        let bytes = GT_LE16.to_le_bytes();
        assert_eq!(bytes[0], b'G');
        assert_eq!(bytes[1], b'T');
        let bytes = AG_LE16.to_le_bytes();
        assert_eq!(bytes[0], b'A');
        assert_eq!(bytes[1], b'G');
    }
}
