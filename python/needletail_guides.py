from __future__ import annotations
from collections import Counter
from dataclasses import replace
from itertools import islice
from typing import TYPE_CHECKING, Iterable, Iterator

if TYPE_CHECKING:
    import needletail
    from seqchain.recipes.crispr import CRISPRPreset
    from seqchain.io.genome import Genome

from seqchain.region import Region


def scan_guides(
    idx: needletail.FmIndex,
    preset: CRISPRPreset,
    genome: Genome,
    *,
    chrom: str | None = None,
) -> Iterator[Region]:
    """Drop-in replacement for regex_map() + interpret_guides().

    All coordinate math, sequence assembly, and hashing done in Rust.
    Python only does: chrom name lookup + Region construction.
    """
    names = idx.chrom_names()
    topo_flags = [genome.topologies.get(n, "linear") == "circular" for n in names]

    buf = idx.scan_guides(preset.pam, preset.spacer_len, preset.pam_direction, topo_flags)

    chrom_filter = None
    if chrom is not None:
        chrom_filter = {i for i, n in enumerate(names) if n == chrom}

    pam_pattern = preset.pam
    pam_len = buf.pam_len
    for i in range(len(buf)):
        cid, pam_pos, gs, ge, strand, spacer_b, pam_b, gseq_b, gid = buf[i]
        if chrom_filter is not None and cid not in chrom_filter:
            continue
        yield Region(
            chrom=names[cid],
            start=gs,
            end=ge,
            strand=strand,
            name=gid,
            tags={
                "pattern": pam_pattern,
                "matched": pam_b.decode(),
                "guide_id": gid,
                "pam_start": pam_pos,
                "pam_end": pam_pos + pam_len,
                "spacer": spacer_b.decode(),
                "guide_seq": gseq_b.decode(),
            },
        )


_CHUNK_SIZE = 100_000


def score_off_targets_fast(
    regions: Iterable[Region],
    index: needletail.FmIndex,
    pam: str,
    *,
    spacer_len: int = 20,
    pam_direction: str = "downstream",
    mismatches: int = 2,
    topologies: dict[str, str] | None = None,
    max_width: int = 8,
    chunk_size: int = _CHUNK_SIZE,
    **kwargs,
) -> Iterator[Region]:
    """Streaming Rust-native off-target scoring with macro-chunking.

    Pulls ``chunk_size`` regions at a time from the upstream generator,
    deduplicates spacers within the chunk, issues a single search_batch
    call with in-Rust PAM validation, and yields scored Regions one by one.

    Memory: O(chunk_size).  Axis 8 compliant — never materializes the
    full stream.  The upstream generator is consumed lazily via islice.

    Drop-in replacement for score_off_targets / score_off_targets_native.
    Extra keyword arguments (e.g. ``sequences``) are accepted and ignored
    for call-site compatibility.
    """
    topo = topologies or {}
    chrom_names = index.chrom_names()
    topo_flags = [topo.get(n, "linear") == "circular" for n in chrom_names]

    it = iter(regions)
    while True:
        chunk = list(islice(it, chunk_size))
        if not chunk:
            return

        # Dedup spacers within this chunk.
        spacer_to_idx: dict[str, int] = {}
        unique_spacers: list[str] = []
        spacer_map: list[int] = []
        for r in chunk:
            spacer = r.tags.get("spacer", "")
            if spacer not in spacer_to_idx:
                spacer_to_idx[spacer] = len(unique_spacers)
                unique_spacers.append(spacer)
            spacer_map.append(spacer_to_idx[spacer])

        # Single Rust call per chunk — GIL released, PAM validated in Rust.
        qi, _pos, _strand, _scores = index.search_batch(
            unique_spacers,
            mismatches=mismatches,
            max_width=max_width,
            pam=pam,
            pam_direction=pam_direction,
            topologies=topo_flags,
        )

        # Count PAM-valid hits per unique spacer.
        counts = Counter(qi)

        for i, r in enumerate(chunk):
            total = counts.get(spacer_map[i], 0)
            off_targets = max(0, total - 1)
            yield replace(
                r,
                score=float(off_targets),
                tags={**r.tags, "total_sites": total, "off_targets": off_targets},
            )

        # Explicit cleanup — release chunk memory before next pull.
        del chunk, spacer_to_idx, unique_spacers, spacer_map, qi, counts
