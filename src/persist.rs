//! Zero-copy mmap persistence for seqchain FM-Index data.
//!
//! Serializes the FM-Index (interleaved rank array, suffix array, less table,
//! chromosome metadata, and genome text) to a `.seqchain` file using rkyv.
//! On load, the file is mmap'd and the archived data is accessed directly —
//! no deserialization for the large rank array or suffix array.
//!
//! File format:
//! ```text
//! Offset 0:   "SQCH" (4 bytes magic)
//! Offset 4:   1u32 LE (format version)
//! Offset 8:   payload_len u64 LE
//! Offset 16:  rkyv payload (16-byte aligned)
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;
use rkyv::rancor::Error as RkyvError;
use rkyv::{Archive, Deserialize, Serialize};

use crate::error::SearchError;
use crate::fm_index::FmIndexSearcher;
use crate::simd_search::{ChromGeometry, FmOcc, BASES};

// ═══════════════════════════════════════════════════════════════════════════════
//  File format constants
// ═══════════════════════════════════════════════════════════════════════════════

const MAGIC: &[u8; 4] = b"SQCH";
const FORMAT_VERSION: u32 = 1;
const HEADER_SIZE: usize = 16;

// ═══════════════════════════════════════════════════════════════════════════════
//  Stored types (rkyv-serializable)
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Archive, Serialize, Deserialize)]
#[rkyv(compare(PartialEq), derive(Debug))]
struct StoredChromInfo {
    name: String,
    start: u64,
    len: u64,
}

#[derive(Archive, Serialize, Deserialize)]
#[rkyv(compare(PartialEq), derive(Debug))]
struct SeqchainIndex {
    rank_data: Vec<u32>,
    less: Vec<u64>,
    sa: Vec<u64>,
    chroms: Vec<StoredChromInfo>,
    text: Vec<u8>,
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Save
// ═══════════════════════════════════════════════════════════════════════════════

/// Serialize an in-memory FM-Index to a `.seqchain` file.
pub fn save_index(searcher: &FmIndexSearcher, path: &Path) -> Result<(), SearchError> {
    let rank_data: Vec<u32> = searcher.to_interleaved_rank_data();
    let less: Vec<u64> = searcher.less_raw().iter().map(|&v| v as u64).collect();
    let sa: Vec<u64> = searcher.sa_raw().iter().map(|&v| v as u64).collect();
    let chroms: Vec<StoredChromInfo> = searcher
        .chroms()
        .iter()
        .map(|c| StoredChromInfo {
            name: c.name.clone(),
            start: c.start as u64,
            len: c.len as u64,
        })
        .collect();
    let text: Vec<u8> = searcher.text().to_vec();

    let index = SeqchainIndex {
        rank_data,
        less,
        sa,
        chroms,
        text,
    };

    let payload = rkyv::to_bytes::<RkyvError>(&index)
        .map_err(|e| SearchError::Other(anyhow::anyhow!("rkyv serialize failed: {}", e)))?;

    let file = File::create(path).map_err(SearchError::Io)?;
    let mut w = BufWriter::new(file);

    // Write header
    w.write_all(MAGIC).map_err(SearchError::Io)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())
        .map_err(SearchError::Io)?;
    w.write_all(&(payload.len() as u64).to_le_bytes())
        .map_err(SearchError::Io)?;

    // Write rkyv payload
    w.write_all(&payload).map_err(SearchError::Io)?;
    w.flush().map_err(SearchError::Io)?;

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════════════
//  MappedIndex — zero-copy mmap'd FM-Index
// ═══════════════════════════════════════════════════════════════════════════════

/// An FM-Index loaded from a `.seqchain` file via mmap.
///
/// The large arrays (rank_data, sa, text) are accessed directly from the
/// memory-mapped file — no deserialization or copying. The small `less` table
/// and chromosome metadata are converted once on load.
pub struct MappedIndex {
    /// The memory map — must outlive the archived pointer.
    _mmap: Mmap,
    /// Pointer into the mmap'd data. Valid for the lifetime of `_mmap`.
    archived: *const ArchivedSeqchainIndex,
    /// Pre-converted less table: `less[byte_value]` = count of chars < byte in BWT.
    less: Box<[usize; 256]>,
    /// Pre-converted chromosome metadata.
    chroms: Vec<(String, usize, usize)>,
}

// SAFETY: The Mmap is read-only and the archived pointer is derived from it.
// No mutation occurs after construction.
unsafe impl Send for MappedIndex {}
unsafe impl Sync for MappedIndex {}

impl MappedIndex {
    /// Access the archived index.
    #[inline]
    fn archived(&self) -> &ArchivedSeqchainIndex {
        // SAFETY: pointer is valid for the lifetime of _mmap, which we own.
        unsafe { &*self.archived }
    }

    /// Chromosome names in FASTA order.
    pub fn chrom_names(&self) -> Vec<String> {
        self.chroms.iter().map(|(name, _, _)| name.clone()).collect()
    }

    /// Chromosome geometry for the search engine.
    pub fn chrom_geometry(&self) -> ChromGeometry {
        ChromGeometry {
            ranges: self.chroms.iter().map(|&(_, start, len)| (start, len)).collect(),
        }
    }

    /// Concatenated genome text.
    pub fn text(&self) -> &[u8] {
        let archived = self.archived();
        &archived.text[..]
    }
}

// ─── FmOcc impl ──────────────────────────────────────────────────────────────

impl FmOcc for MappedIndex {
    #[inline]
    fn less(&self, c: u8) -> usize {
        self.less[c as usize]
    }

    #[inline]
    fn occ(&self, pos: usize, c: u8) -> usize {
        let bi = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
        let archived = self.archived();
        archived.rank_data[pos * 4 + bi].to_native() as usize
    }

    #[inline]
    fn sa(&self, idx: usize) -> usize {
        let archived = self.archived();
        archived.sa[idx].to_native() as usize
    }

    #[inline]
    fn sa_len(&self) -> usize {
        let archived = self.archived();
        archived.sa.len()
    }

    #[inline]
    fn lf_map(&self, l: usize, r: usize) -> [(u32, u32); 4] {
        let archived = self.archived();
        let data = &archived.rank_data[..];

        let occ_l = if l > 0 {
            let b = (l - 1) * 4;
            [
                data[b].to_native(),
                data[b + 1].to_native(),
                data[b + 2].to_native(),
                data[b + 3].to_native(),
            ]
        } else {
            [0u32; 4]
        };

        let b = r * 4;
        let occ_r = [
            data[b].to_native(),
            data[b + 1].to_native(),
            data[b + 2].to_native(),
            data[b + 3].to_native(),
        ];

        let mut result = [(0u32, 0u32); 4];
        for bi in 0..4 {
            let less_b = self.less[BASES[bi] as usize] as u32;
            result[bi] = (less_b + occ_l[bi], less_b + occ_r[bi]);
        }
        result
    }

    #[inline]
    fn prefetch_lf(&self, l: usize, r: usize) {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
            let archived = self.archived();
            let data = &archived.rank_data[..];
            if l > 0 {
                let ptr_l = (data.as_ptr() as *const u8).add((l - 1) * 16) as *const i8;
                _mm_prefetch(ptr_l, _MM_HINT_T0);
            }
            let ptr_r = (data.as_ptr() as *const u8).add(r * 16) as *const i8;
            _mm_prefetch(ptr_r, _MM_HINT_T0);
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (l, r);
        }
    }

    #[inline]
    fn rank_data(&self) -> Option<(&[u32], &[usize; 256])> {
        let archived = self.archived();
        // On little-endian, ArchivedVec<u32> is bit-identical to &[u32].
        // We transmute the slice of archived u32 to native u32.
        let data: &[rkyv::rend::u32_le] = &archived.rank_data[..];
        // SAFETY: u32_le has the same layout as u32 on LE platforms.
        #[cfg(target_endian = "little")]
        {
            let native_data: &[u32] =
                unsafe { std::slice::from_raw_parts(data.as_ptr() as *const u32, data.len()) };
            Some((native_data, &self.less))
        }
        #[cfg(not(target_endian = "little"))]
        {
            let _ = data;
            None // Fall back to scalar lf_map on big-endian
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
//  Load
// ═══════════════════════════════════════════════════════════════════════════════

/// Load a `.seqchain` file via mmap. Returns a zero-copy `MappedIndex`.
pub fn load_index(path: &Path) -> Result<MappedIndex, SearchError> {
    let file = File::open(path).map_err(SearchError::Io)?;
    let mmap = unsafe { Mmap::map(&file).map_err(SearchError::Io)? };

    // Validate header
    if mmap.len() < HEADER_SIZE {
        return Err(SearchError::Other(anyhow::anyhow!(
            "file too small for header: {} bytes",
            mmap.len()
        )));
    }

    if &mmap[0..4] != MAGIC {
        return Err(SearchError::Other(anyhow::anyhow!(
            "invalid magic bytes (not a .seqchain file)"
        )));
    }

    let version = u32::from_le_bytes(mmap[4..8].try_into().unwrap());
    if version != FORMAT_VERSION {
        return Err(SearchError::Other(anyhow::anyhow!(
            "unsupported format version: {} (expected {})",
            version,
            FORMAT_VERSION
        )));
    }

    let payload_len = u64::from_le_bytes(mmap[8..16].try_into().unwrap()) as usize;
    if mmap.len() < HEADER_SIZE + payload_len {
        return Err(SearchError::Other(anyhow::anyhow!(
            "file truncated: expected {} payload bytes, have {}",
            payload_len,
            mmap.len() - HEADER_SIZE
        )));
    }

    let payload = &mmap[HEADER_SIZE..HEADER_SIZE + payload_len];

    // Validate and access the archived data
    let archived: &ArchivedSeqchainIndex =
        rkyv::access::<ArchivedSeqchainIndex, RkyvError>(payload)
            .map_err(|e| SearchError::Other(anyhow::anyhow!("rkyv validation failed: {}", e)))?;

    // Convert less table (256 × u64 → 256 × usize): 2KB, once
    let mut less = Box::new([0usize; 256]);
    for (i, val) in archived.less.iter().enumerate() {
        if i < 256 {
            less[i] = val.to_native() as usize;
        }
    }

    // Convert chromosome metadata (tiny)
    let chroms: Vec<(String, usize, usize)> = archived
        .chroms
        .iter()
        .map(|c| {
            (
                c.name.as_str().to_string(),
                c.start.to_native() as usize,
                c.len.to_native() as usize,
            )
        })
        .collect();

    let ptr = archived as *const ArchivedSeqchainIndex;

    Ok(MappedIndex {
        _mmap: mmap,
        archived: ptr,
        less,
        chroms,
    })
}
