"""Canonical output schema for Layer 0 (and the shape later layers should converge on).

CONTRIBUTING.md rule #1 names this file as load-bearing: every layer should produce output
conforming to one shared shape rather than inventing its own column names. Layer 0 populates
"computational" (VEP/interval overlap) and "curated_clinical" (ClinVar) observations;
Layer 1 (GWAS/eQTL/MPRA/CRISPR) and Layer 3 (AlphaGenome) emit records conforming to this
unified EvidenceRecord schema with explicit dataset_id and dataset_version provenance.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import pandas as pd

# Columns every Layer 0 output DataFrame must carry. Extend this list (don't replace it) when a
# new layer adds a column that should become part of the shared contract.
LAYER0_REQUIRED_COLUMNS = [
    "chrom", "pos", "ref", "alt",
    "label", "CLNSIG", "CLNREVSTAT", "review_stars",
    "interval_evidence",
    "vep_canonical_terms", "vep_all_transcript_terms", "vep_regulatory_terms",
    "regulatory_evidence", "is_regulatory",
]

SOURCE_TYPES = ("experimental", "computational", "curated_clinical")
MISSINGNESS_FLAGS = ("observed_positive", "observed_null", "absent")
COMMON_MODALITIES = (
    "clinvar",
    "ccre_overlap",
    "vep_consequence",
    "gwas",
    "eqtl",
    "gnomad_af",
    "mpra",
    "crispr",
    "alphagenome",
)


@dataclass
class EvidenceRecord:
    """One (variant, evidence-source) observation, in the canonical framework shape.

    Unified schema carrying explicit dataset_id and dataset_version provenance.
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    modality: str
    source_type: str
    evidence_value: str
    context: Optional[str] = None
    missingness_flag: str = "observed_positive"
    dataset_id: Optional[str] = None
    dataset_version: Optional[str] = None

    def __post_init__(self) -> None:
        if self.source_type not in SOURCE_TYPES:
            raise ValueError(f"source_type must be one of {SOURCE_TYPES}, got {self.source_type!r}")
        if self.missingness_flag not in MISSINGNESS_FLAGS:
            raise ValueError(f"missingness_flag must be one of {MISSINGNESS_FLAGS}, got {self.missingness_flag!r}")


def validate_layer0_output(df: pd.DataFrame) -> None:
    """Raise loudly if a Layer 0 output DataFrame is missing a required column.

    Call this at the end of run_layer0() so schema drift is caught at the source, not discovered
    later by a KeyError three layers downstream.
    """
    missing = [c for c in LAYER0_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Layer 0 output missing required column(s): {missing}")


def layer0_to_evidence_records(df: pd.DataFrame, clinvar_version: str = "latest") -> List[EvidenceRecord]:
    """Convert Layer 0 wide DataFrame rows into canonical EvidenceRecord rows.

    Emits standardized records for:
    - ClinVar clinical classification (curated_clinical / clinvar)
    - Interval / cCRE overlaps (computational / ccre_overlap)
    - Ensembl VEP regulatory consequence terms (computational / vep_consequence)
    """
    records: List[EvidenceRecord] = []
    for row in df.to_dict("records"):
        chrom = str(row["chrom"]).replace("chr", "")
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])

        # 1. ClinVar clinical label
        if pd.notna(row.get("label")):
            cv_val = json.dumps({
                "label": row.get("label"),
                "clnsig": row.get("CLNSIG"),
                "clnrevstat": row.get("CLNREVSTAT"),
                "review_stars": row.get("review_stars"),
            }, sort_keys=True)
            records.append(
                EvidenceRecord(
                    chrom=chrom, pos=pos, ref=ref, alt=alt,
                    modality="clinvar",
                    source_type="curated_clinical",
                    evidence_value=cv_val,
                    context="clinvar",
                    missingness_flag="observed_positive",
                    dataset_id="clinvar",
                    dataset_version=str(clinvar_version),
                )
            )

        # 2. cCRE / Interval evidence
        intervals = row.get("interval_evidence") or []
        if isinstance(intervals, str):
            try:
                intervals = json.loads(intervals)
            except Exception:
                intervals = [intervals]
        if intervals:
            records.append(
                EvidenceRecord(
                    chrom=chrom, pos=pos, ref=ref, alt=alt,
                    modality="ccre_overlap",
                    source_type="computational",
                    evidence_value=json.dumps({"intervals": intervals}, sort_keys=True),
                    context="encode_ccre",
                    missingness_flag="observed_positive",
                    dataset_id="encode_ccre",
                    dataset_version="v4",
                )
            )

        # 3. VEP Regulatory Consequence
        vep_terms = row.get("vep_regulatory_terms") or row.get("vep_canonical_terms") or []
        if isinstance(vep_terms, str):
            try:
                vep_terms = json.loads(vep_terms)
            except Exception:
                vep_terms = [vep_terms]
        if vep_terms:
            records.append(
                EvidenceRecord(
                    chrom=chrom, pos=pos, ref=ref, alt=alt,
                    modality="vep_consequence",
                    source_type="computational",
                    evidence_value=json.dumps({"terms": vep_terms}, sort_keys=True),
                    context="ensembl_vep",
                    missingness_flag="observed_positive",
                    dataset_id="ensembl_vep",
                    dataset_version="release_110",
                )
            )
    return records
