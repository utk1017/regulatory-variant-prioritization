"""Variant identity helpers for Layer 1.

canon_key normalises both alleles so that the same indel written in a
non-minimal representation by two different sources (GTEx ``variantId``,
gnomAD ``variant_id``, a local MPRA table) still produces the same lookup key.

plan_windows clusters variant positions into dense genomic sub-windows so
region-based API queries (gnomAD GraphQL, GWAS Catalog REST) only span
actually occupied regions rather than the full flank interval.
"""
from __future__ import annotations

from collections import namedtuple
from typing import List, Sequence, Tuple

import pandas as pd


Variant = namedtuple("Variant", "chrom pos ref alt key")


def canon_key(
    chrom: str | int,
    pos: int | str,
    ref: str,
    alt: str,
) -> Tuple[str, int, str, str]:
    """Return a normalised (chrom, pos, ref, alt) 4-tuple.

    Trims shared suffix then shared prefix so non-minimal indels match their
    minimal representations.  Does *not* left-align (that requires the
    reference FASTA), so indel matching remains best-effort.
    """
    ref, alt, pos = str(ref).upper(), str(alt).upper(), int(pos)
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return (str(chrom).replace("chr", ""), pos, ref, alt)


def plan_windows(
    positions: Sequence[int], max_span: int
) -> List[Tuple[int, int]]:
    """Cluster sorted genomic positions into windows of at most *max_span* bp.

    Returns a list of ``(start, stop)`` pairs where every position is covered
    by exactly one window.
    """
    windows: list[list[int]] = []
    for p in sorted(set(positions)):
        if windows and p - windows[-1][0] < max_span:
            windows[-1][1] = p
        else:
            windows.append([p, p])
    return [(s, e) for s, e in windows]


def to_variants(df: pd.DataFrame) -> List[Variant]:
    """Convert a Layer 0 candidate DataFrame into a list of :class:`Variant`."""
    out = []
    for row in df[["chrom", "pos", "ref", "alt"]].itertuples(index=False):
        chrom = str(row.chrom).replace("chr", "")
        out.append(
            Variant(
                chrom,
                int(row.pos),
                str(row.ref),
                str(row.alt),
                canon_key(chrom, row.pos, row.ref, row.alt),
            )
        )
    return out
