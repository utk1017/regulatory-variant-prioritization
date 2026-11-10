"""Evidence stratification and calibration tier assignment for Layer 2.

Assigns variants into evidence-availability tiers (Rich, Moderate, Sparse)
based on the count of tested experimental regulatory modalities, and performs
balance and label coverage checks needed for Layer 6 calibration.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from regvar.layer2.reasons import aggregate_modality_status, filter_records_to_targets
from regvar.layer2.schema import EVIDENCE_TIERS, LAYER2_AVAILABILITY_COLUMNS

DEFAULT_TIER_THRESHOLDS: Dict[str, int] = {
    "tier_1_min_modalities": 2,  # >= 2 experimental regulatory modalities tested -> Tier 1
    "tier_2_modalities": 1,      # exactly 1 modality tested -> Tier 2
    "tier_3_modalities": 0,      # 0 modalities tested -> Tier 3
}


def compute_evidence_tiers(
    variants: Sequence[Any],
    evidence_by_key: Dict[Tuple[str, int, str, str], List[Dict[str, Any]]],
    clinvar_labels: Optional[Dict[Tuple[str, int, str, str], Dict[str, Any]]] = None,
    regulatory_modalities: Sequence[str] = ("eqtl", "mpra", "crispr"),
    target_tissues: Optional[Sequence[str]] = None,
    target_cell_types: Optional[Sequence[str]] = None,
    tier_thresholds: Optional[Dict[str, int]] = None,
) -> pd.DataFrame:
    """Assemble per-variant availability vector and assign evidence tiers.

    Tier counts use allele-resolved experimental modalities (eQTL, MPRA, CRISPR).
    GWAS Catalog evidence is retained as locus/position-level context, but is not counted
    as allele-specific evidence because the Layer 1 match is positional, not allele-specific.
    """
    thresholds = tier_thresholds or DEFAULT_TIER_THRESHOLDS
    rows: List[Dict[str, Any]] = []

    for v in variants:
        records_for_var = filter_records_to_targets(
            evidence_by_key.get(v.key, []),
            target_tissues=target_tissues,
            target_cell_types=target_cell_types,
        )
        by_modality: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for r in records_for_var:
            by_modality[r["modality"]].append(r)

        status_dict: Dict[str, str] = {}
        all_reasons: List[str] = []
        n_pos = 0
        n_null = 0
        n_abs = 0
        n_tested_experimental = 0

        for mod in ("gwas", "eqtl", "mpra", "crispr", "gnomad_af"):
            records = by_modality.get(mod, [])
            status, cat, reasons, _ = aggregate_modality_status(records)
            status_dict[f"{mod}_status"] = status

            if mod in regulatory_modalities:
                if status == "observed_positive":
                    n_pos += 1
                    n_tested_experimental += 1
                elif status == "observed_null":
                    n_null += 1
                    n_tested_experimental += 1
                else:
                    n_abs += 1

            if reasons and status == "absent":
                all_reasons.extend(reasons)

        # Tier assignment:
        if n_tested_experimental >= thresholds.get("tier_1_min_modalities", 2):
            tier = "Tier 1 (Evidence-Rich)"
        elif n_tested_experimental == thresholds.get("tier_2_modalities", 1):
            tier = "Tier 2 (Moderate Evidence)"
        else:
            tier = "Tier 3 (Evidence-Sparse)"

        row: Dict[str, Any] = {
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "variant_key": f"{v.chrom}:{v.pos}:{v.ref}:{v.alt}",
            "gwas_status": status_dict.get("gwas_status", "absent"),
            "eqtl_status": status_dict.get("eqtl_status", "absent"),
            "mpra_status": status_dict.get("mpra_status", "absent"),
            "crispr_status": status_dict.get("crispr_status", "absent"),
            "gnomad_af_status": status_dict.get("gnomad_af_status", "absent"),
            "n_positive": n_pos,
            "n_null": n_null,
            "n_absent": n_abs,
            "n_tested_experimental": n_tested_experimental,
            "evidence_tier": tier,
            "primary_gap_reasons": ";".join(sorted(set(all_reasons))) if all_reasons else "none",
        }

        if clinvar_labels and v.key in clinvar_labels:
            cl = clinvar_labels[v.key]
            row["CLNSIG"] = cl.get("CLNSIG")
            row["label"] = cl.get("label")
            row["review_stars"] = cl.get("review_stars")

        rows.append(row)

    df = pd.DataFrame(rows)
    extra_cols = [c for c in ["CLNSIG", "label", "review_stars"] if clinvar_labels]
    for c in extra_cols:
        if c not in df.columns:
            df[c] = None

    return df[LAYER2_AVAILABILITY_COLUMNS + extra_cols]


def tier_balance_warnings(df: pd.DataFrame, min_per_tier: int = 5) -> List[str]:
    """Flag empty or sparse tiers prior to downstream Layer 6 calibration."""
    out: List[str] = []
    counts = df["evidence_tier"].value_counts()
    for tier in EVIDENCE_TIERS:
        n = int(counts.get(tier, 0))
        if n == 0:
            out.append(f"{tier}: 0 variants -- any Layer 6 comparison involving this tier is impossible on this run.")
        elif n < min_per_tier:
            out.append(f"{tier}: only {n} variant(s) -- too few for stable per-tier calibration (want >= {min_per_tier}).")
    return out


def label_coverage_warnings(df: pd.DataFrame) -> List[str]:
    """Ensure ClinVar pathogenic and benign ground truth labels exist within each tier."""
    if "label" not in df.columns or df["label"].isna().all():
        return [
            "No ClinVar pathogenic/benign labels attached to any variant -- Layer 6 cannot calibrate. "
            "Check that Layer 0's CSV has a 'label' column and that its variant keys match Layer 1's."
        ]
    out: List[str] = []
    for tier, sub in df.groupby("evidence_tier"):
        lab = sub["label"].dropna()
        for cls in ("pathogenic", "benign"):
            if (lab == cls).sum() == 0:
                out.append(f"{tier}: no '{cls}' variants -- per-tier calibration needs both classes.")
    return out
