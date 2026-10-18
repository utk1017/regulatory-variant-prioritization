"""Input candidate intake for Layer 1.

Discovers Layer 0 candidate truth-sets from either:
1. Explicit CSV path in config
2. Repo layout: ``data/outputs/<gene_symbol>/layer0_regulatory_truthset.csv``
3. Legacy notebook layout: ``layer0_outputs/<gene_symbol>_layer0_regulatory_truthset.csv``

Reads gene metadata from Layer 0's manifest when present to ensure both layers
remain pinned to the exact same gene definition and version.
"""
from __future__ import annotations

import json
import os
from typing import Any, Callable, Dict, Optional, Tuple

import pandas as pd

from regvar.common.provenance import _file_hash
from regvar.layer0.gene_meta import fetch_gene_metadata


def find_layer0_inputs(
    gene_symbol: str, explicit_csv: Optional[str] = None
) -> Tuple[str, Optional[str], str]:
    """Find Layer 0 output CSV and optional manifest path.

    Returns ``(csv_path, manifest_path, layout_label)``.
    Raises FileNotFoundError if no candidate CSV exists.
    """
    candidates = []
    if explicit_csv:
        candidates.append(
            (
                explicit_csv,
                os.path.join(os.path.dirname(explicit_csv), "manifest.json"),
                "config path",
            )
        )
    candidates.append(
        (
            f"data/outputs/{gene_symbol}/layer0_regulatory_truthset.csv",
            f"data/outputs/{gene_symbol}/manifest.json",
            "repo layout",
        )
    )
    candidates.append(
        (
            f"layer0_outputs/{gene_symbol}_layer0_regulatory_truthset.csv",
            None,
            "layer-0 notebook layout",
        )
    )
    for csv_path, manifest_path, label in candidates:
        if os.path.exists(csv_path):
            return csv_path, manifest_path, label
    raise FileNotFoundError(
        "No Layer 0 candidate CSV found. Looked in: "
        + ", ".join(c[0] for c in candidates)
    )


def read_candidate_csv(csv_path: str) -> pd.DataFrame:
    """Read Layer 0 candidate variants and enforce string chrom and integer pos."""
    df = pd.read_csv(csv_path, dtype={"chrom": str})
    need = ["chrom", "pos", "ref", "alt"]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{csv_path} is missing required column(s): {missing}")
    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = df["pos"].astype(int)
    return (
        df.drop_duplicates(subset=need)
        .sort_values(["chrom", "pos"])
        .reset_index(drop=True)
    )


def load_layer0_output(
    gene_symbol: str,
    explicit_csv: Optional[str] = None,
    fetch: Optional[Callable] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any], Dict[str, Any]]:
    """Load Layer 0 variants and pinned gene metadata.

    Returns ``(variants_df, gene_meta, info_dict)``.
    """
    csv_path, manifest_path, layout = find_layer0_inputs(gene_symbol, explicit_csv)
    variants_df = read_candidate_csv(csv_path)
    info = {
        "csv_path": csv_path,
        "csv_sha256": _file_hash(csv_path),
        "layout": layout,
        "manifest_path": None,
        "n_variants": len(variants_df),
    }

    if manifest_path and os.path.exists(manifest_path):
        with open(manifest_path, "r", encoding="utf-8") as fh:
            l0_manifest = json.load(fh)
        gene_meta = l0_manifest["gene_metadata"]
        info["manifest_path"] = manifest_path
        info["clinvar_snapshot_date"] = l0_manifest.get("clinvar_snapshot_date")
    else:
        print(
            "WARNING: no Layer 0 manifest found -- fetching gene metadata fresh from Ensembl. "
            "This is NOT pinned to the gene build/version Layer 0 used."
        )
        gene_meta = fetch_gene_metadata(gene_symbol)
        info["gene_metadata_source"] = "fresh Ensembl lookup (not pinned to Layer 0)"

    return variants_df, gene_meta, info
