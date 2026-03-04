//! Streaming FASTQ parser: O(C × L) memory, never materialises the full file.
//!
//! `ChunkedFastq<R>` is an iterator over `io::Result<Vec<FastqRecord>>`.
//! Each `Vec` is bounded by `chunk_size` records. The previous chunk is freed
//! by the caller before `next()` is called — Rust ownership enforces this.
//!
//! Import hierarchy: no engine, no biology. Pure I/O boundary.

use std::io::{self, BufRead};

// ═══════════════════════════════════════════════════════════════════════════
//  FastqRecord
// ═══════════════════════════════════════════════════════════════════════════

/// One FASTQ record with owned, ASCII-uppercase sequence bytes.
///
/// Owning the data (rather than borrowing from a line buffer) means chunks
/// can be passed to parallel threads without lifetime entanglement.
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Read identifier — everything after `@` up to the first whitespace.
    pub id:   Vec<u8>,
    /// Nucleotide sequence, ASCII-uppercased. Safe to pass to BWT search as `&[u8]`.
    pub seq:  Vec<u8>,
    /// Phred+33 quality string. Same length as `seq`.
    pub qual: Vec<u8>,
}

impl FastqRecord {
    /// Read length in bases.
    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  ChunkedFastq
// ═══════════════════════════════════════════════════════════════════════════

/// Streaming FASTQ parser.
///
/// Yields chunks of at most `chunk_size` records. Exactly one line buffer
/// lives inside the parser; all other memory is owned by the yielded chunks.
///
/// # Memory invariant
/// At any moment the parser holds O(max_line_len) bytes internally.
/// The caller holds O(chunk_size × read_len) bytes in the live chunk.
/// When the chunk is dropped (end of loop body), those bytes are freed.
pub struct ChunkedFastq<R: BufRead> {
    inner:      R,
    line:       Vec<u8>,  // re-used line buffer: O(max_line_len)
    chunk_size: usize,
    exhausted:  bool,
}

impl<R: BufRead> ChunkedFastq<R> {
    /// Wrap any `BufRead` source. `chunk_size` must be > 0.
    pub fn new(inner: R, chunk_size: usize) -> Self {
        assert!(chunk_size > 0, "chunk_size must be > 0");
        Self {
            inner,
            line: Vec::with_capacity(512),
            chunk_size,
            exhausted: false,
        }
    }
}

impl<R: BufRead> Iterator for ChunkedFastq<R> {
    type Item = io::Result<Vec<FastqRecord>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }

        // Allocate chunk storage once; previous chunk was freed by caller.
        let mut chunk = Vec::with_capacity(self.chunk_size);

        for _ in 0..self.chunk_size {
            // ── @id line ────────────────────────────────────────────────────
            self.line.clear();
            match self.inner.read_until(b'\n', &mut self.line) {
                Ok(0) => {
                    self.exhausted = true;
                    break;
                }
                Err(e) => return Some(Err(e)),
                Ok(_) => {}
            }
            let raw = trim_newline(&self.line);
            // Skip blank lines or non-record lines.
            if raw.is_empty() {
                continue;
            }
            if raw[0] != b'@' {
                // Malformed FASTQ — skip to resync.
                continue;
            }
            // QNAME: id up to the first whitespace, strip leading '@'.
            let id = qname_from_id_line(&raw[1..]);

            // ── sequence line ────────────────────────────────────────────
            self.line.clear();
            if let Err(e) = self.inner.read_until(b'\n', &mut self.line) {
                return Some(Err(e));
            }
            let mut seq = trim_newline(&self.line).to_vec();
            seq.make_ascii_uppercase();

            // ── '+' separator (discard) ──────────────────────────────────
            self.line.clear();
            if let Err(e) = self.inner.read_until(b'\n', &mut self.line) {
                return Some(Err(e));
            }

            // ── quality line ─────────────────────────────────────────────
            self.line.clear();
            if let Err(e) = self.inner.read_until(b'\n', &mut self.line) {
                return Some(Err(e));
            }
            let qual = trim_newline(&self.line).to_vec();

            if seq.is_empty() {
                continue;
            }
            chunk.push(FastqRecord { id, seq, qual });
        }

        if chunk.is_empty() {
            None
        } else {
            Some(Ok(chunk))
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Private helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Strip trailing `\r` and `\n` without allocating.
#[inline]
fn trim_newline(buf: &[u8]) -> &[u8] {
    let end = buf
        .iter()
        .rposition(|&b| b != b'\n' && b != b'\r')
        .map(|i| i + 1)
        .unwrap_or(0);
    &buf[..end]
}

/// Extract QNAME: everything up to the first ASCII whitespace.
/// SAM §1.4: QNAME must not contain whitespace.
#[inline]
fn qname_from_id_line(line: &[u8]) -> Vec<u8> {
    let end = line
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(line.len());
    line[..end].to_vec()
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const FASTQ: &[u8] = b"\
@read1 some comment\n\
ACGTACGT\n\
+\n\
IIIIIIII\n\
@read2\n\
TTTTGGGG\n\
+\n\
HHHHHHHH\n\
@read3\n\
CCCCAAAA\n\
+\n\
FFFFFFFF\n";

    #[test]
    fn parses_all_records() {
        let reader = Cursor::new(FASTQ);
        let chunks: Vec<_> = ChunkedFastq::new(reader, 10)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].len(), 3);
    }

    #[test]
    fn qname_strips_comment() {
        let reader = Cursor::new(FASTQ);
        let mut parser = ChunkedFastq::new(reader, 10);
        let chunk = parser.next().unwrap().unwrap();
        assert_eq!(chunk[0].id, b"read1");
        assert_eq!(chunk[1].id, b"read2");
    }

    #[test]
    fn sequence_uppercased() {
        let fq = b"@r\nacgt\n+\nIIII\n";
        let reader = Cursor::new(fq);
        let mut parser = ChunkedFastq::new(reader, 10);
        let chunk = parser.next().unwrap().unwrap();
        assert_eq!(chunk[0].seq, b"ACGT");
    }

    #[test]
    fn respects_chunk_size() {
        let reader = Cursor::new(FASTQ);
        let chunks: Vec<_> = ChunkedFastq::new(reader, 2)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].len(), 2);
        assert_eq!(chunks[1].len(), 1);
    }
}
