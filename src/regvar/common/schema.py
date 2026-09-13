"""Canonical output schema for Layer 0 (and the shape later layers should converge on).

CONTRIBUTING.md rule #1 names this file as load-bearing: every layer should produce output
conforming to one shared shape rather than inventing its own column names. Layer 0 today only
populates the "computational_annotation" (VEP/interval overlap) and "curated_clinical" (ClinVar)
source types; Layer 1 (GWAS/eQTL/MPRA/CRISPR) and Layer 3 (AlphaGenome) are expected to emit
records identifiable by the same source_type/modality vocabulary -- see
docs/framework_overview.md, Layers 1 and 3.
"""
from dataclasses import dataclass
from typing import Optional

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


@dataclass
class EvidenceRecord:
    """One (variant, evidence-source) observation, in the shape later layers should also emit.

    This is intentionally a plain dataclass, not a full ORM/pydantic model -- the goal is a single
    place that names the fields, not a validation framework. `validate_layer0_output` below is
    the actual enforcement point.
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    modality: str                 # e.g. "clinvar", "ccre_overlap", "vep_consequence"
    source_type: str              # one of SOURCE_TYPES
    evidence_value: str           # the tag/value itself, e.g. "dELS", "5_prime_UTR_variant"
    context: Optional[str] = None
    missingness_flag: str = "observed_positive"

    def __post_init__(self):
        if self.source_type not in SOURCE_TYPES:
            raise ValueError(f"source_type must be one of {SOURCE_TYPES}, got {self.source_type!r}")
        if self.missingness_flag not in MISSINGNESS_FLAGS:
            raise ValueError(f"missingness_flag must be one of {MISSINGNESS_FLAGS}, got {self.missingness_flag!r}")


def validate_layer0_output(df) -> None:
    """Raise loudly if a Layer 0 output DataFrame is missing a required column.

    Call this at the end of run_layer0() so schema drift is caught at the source, not discovered
    later by a KeyError three layers downstream.
    """
    missing = [c for c in LAYER0_REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Layer 0 output missing required column(s): {missing}")
