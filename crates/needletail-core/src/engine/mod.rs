//! The Mechanics — mechanical data structures for exact-match search.
//!
//! This layer owns the BWT, suffix arrays, rank structures, k-mer seed
//! tables, and the width-first SIMD traversal engine. It knows how to
//! navigate a compressed text index but has no concept of biology,
//! chemistry, or Python.
//!
//! Import hierarchy: `geometry` (never chemistry, operations, io, or lib).

pub mod affine;
pub mod fm_index;
pub mod kmer_index;
pub mod simd_search;
