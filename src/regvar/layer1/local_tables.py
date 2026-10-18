"""Local evidence table loaders for MPRA and CRISPR screens.

No live query API is used for these modalities (they are downloadable catalogs
or experimental files). The local table defines what was tested:

- Variant in table, p <= threshold -> observed_positive
- Variant in table, p > threshold -> observed_null (tested, not significant)
- Variant not in table -> absent (not_tested_in_provided_table)
- No path configured -> absent (no_data_source_configured)
- Configured path that does not exist -> raises FileNotFoundError
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd

from regvar.layer1.schema import make_record
from regvar.layer1.variants import Variant, canon_key


MPRA_COLUMNS = [
    "chrom", "pos", "ref", "alt",
    "cell_type", "log2_allelic_skew", "pvalue", "dataset_name",
]

CRISPR_COLUMNS = [
    "chrom", "pos", "ref", "alt",
    "target_gene", "cell_type", "effect_size", "pvalue",
    "screen_type", "dataset_name",
]


def load_local_evidence_table(
    path: Optional[str],
    expected_columns: List[str],
    label: str,
    numeric_columns: Tuple[str, ...] = (),
) -> Optional[Dict[Tuple[str, int, str, str], List[dict]]]:
    """Load a local CSV/TSV table and index it by canon_key.

    Returns ``{canon_key: [row_dict, ...]}`` or ``None`` when path is None.
    Raises FileNotFoundError if path is configured but missing.
    Raises ValueError if expected columns are missing.
    """
    if not path:
        return None
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{label} evidence path is configured but the file does not exist: {path}"
        )
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    missing = [c for c in expected_columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} evidence table at {path} is missing expected column(s): {missing}. "
            f"Expected: {expected_columns}"
        )

    df["chrom"] = df["chrom"].astype(str).str.replace("chr", "", regex=False)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    n_bad = int(df["pos"].isna().sum())
    if n_bad:
        print(f"WARNING: {label} table has {n_bad} row(s) with a non-numeric pos; dropped.")
        df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    df["ref"], df["alt"] = df["ref"].str.upper(), df["alt"].str.upper()

    for col in ["pvalue", *numeric_columns]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    index: Dict[Tuple[str, int, str, str], List[dict]] = defaultdict(list)
    for row in df.to_dict("records"):
        index[canon_key(row["chrom"], row["pos"], row["ref"], row["alt"])].append(row)
    return dict(index)


def build_local_table_records(
    modality: str,
    variants: List[Variant],
    index: Optional[Dict[Tuple[str, int, str, str], List[dict]]],
    p_threshold: float,
    value_fn: Callable[[dict], dict],
):
    """Generate EvidenceRecords for local tables."""
    records = []
    for v in variants:
        if index is None:
            records.append(
                make_record(
                    v, modality, "experimental", "absent",
                    reason="no_data_source_configured",
                )
            )
            continue
        rows = index.get(v.key)
        if not rows:
            records.append(
                make_record(
                    v, modality, "experimental", "absent",
                    reason="not_tested_in_provided_table",
                )
            )
            continue
        for row in rows:
            p, ctx = row.get("pvalue"), row.get("cell_type")
            if p is None or p != p:
                records.append(
                    make_record(
                        v, modality, "experimental", "absent",
                        context=ctx,
                        reason="p_value_missing",
                    )
                )
            elif p <= p_threshold:
                records.append(
                    make_record(
                        v, modality, "experimental", "observed_positive",
                        context=ctx,
                        **value_fn(row),
                    )
                )
            else:
                records.append(
                    make_record(
                        v, modality, "experimental", "observed_null",
                        context=ctx,
                        reason="tested_not_significant",
                        pvalue=p,
                    )
                )
    return records


def mpra_value(row: dict) -> dict:
    return {
        "log2_allelic_skew": row.get("log2_allelic_skew"),
        "pvalue": row.get("pvalue"),
        "dataset": row.get("dataset_name"),
    }


def crispr_value(row: dict) -> dict:
    return {
        "target_gene": row.get("target_gene"),
        "effect_size": row.get("effect_size"),
        "pvalue": row.get("pvalue"),
        "screen_type": row.get("screen_type"),
        "dataset": row.get("dataset_name"),
    }
