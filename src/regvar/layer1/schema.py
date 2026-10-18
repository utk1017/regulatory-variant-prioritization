"""Layer 1 output schema and validation.

Layer 1 produces one EvidenceRecord per (variant, modality) observation.
All callers must use :func:`make_record` rather than constructing
:class:`EvidenceRecord` directly -- this ensures the ``evidence_value`` JSON
is always produced by the NaN-safe :func:`ev` serialiser.

MODALITY_META tells downstream layers what they may do with each modality:

- Layer 2 (Rules): reads ``completable`` to decide which absent modalities
  should be sent to AlphaGenome for computational completion.
- Layer 5 (Bayesian fusion): reads ``use_in_fusion`` to decide which
  modalities contribute to the fused regulatory score.
- Layer 6 (Calibration): uses ``gnomad_af`` (``use_in_fusion=False``) for
  stratification only, to avoid label-calibration leakage.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any, List, Optional

import pandas as pd


SOURCE_TYPES = ("experimental", "computational", "curated_clinical")
MISSINGNESS_FLAGS = ("observed_positive", "observed_null", "absent")
MODALITIES = ("gwas", "eqtl", "gnomad_af", "mpra", "crispr")

MODALITY_META = {
    "gwas": {
        "kind": "association",
        "completable": False,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
    },
    "eqtl": {
        "kind": "expression_effect",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
    },
    "mpra": {
        "kind": "reporter_activity",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
    },
    "crispr": {
        "kind": "perturbation_effect",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
    },
    "gnomad_af": {
        "kind": "population_frequency",
        "completable": False,
        "use_in_fusion": False,
        "direction": "benign_leaning_when_observed",
    },
}

LAYER1_COLUMNS = [
    "chrom", "pos", "ref", "alt",
    "modality", "source_type",
    "evidence_value", "context", "missingness_flag",
]


# ---------------------------------------------------------------------------
# JSON serialisation helper
# ---------------------------------------------------------------------------

def _clean(obj: Any) -> Any:
    """Recursively make *obj* JSON-safe (NaN/inf -> None, numpy scalars -> Python)."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_clean(v) for v in obj]
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        try:
            obj = obj.item()
        except Exception:
            pass
    if isinstance(obj, float) and (obj != obj or obj in (float("inf"), float("-inf"))):
        return None
    return obj


def ev(**kw) -> str:
    """Serialise keyword arguments as compact JSON, safely handling NaN/inf."""
    return json.dumps(_clean(kw), sort_keys=True, allow_nan=False)


# ---------------------------------------------------------------------------
# Record dataclass
# ---------------------------------------------------------------------------


@dataclass
class EvidenceRecord:
    """One (variant, modality) observation produced by Layer 1.

    ``evidence_value`` must be compact JSON (always produced via :func:`ev`).
    Non-positive rows carry a ``reason`` key inside ``evidence_value``.
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

    def __post_init__(self) -> None:
        if self.modality not in MODALITIES:
            raise ValueError(
                f"modality must be one of {MODALITIES}, got {self.modality!r}"
            )
        if self.source_type not in SOURCE_TYPES:
            raise ValueError(
                f"source_type must be one of {SOURCE_TYPES}, got {self.source_type!r}"
            )
        if self.missingness_flag not in MISSINGNESS_FLAGS:
            raise ValueError(
                f"missingness_flag must be one of {MISSINGNESS_FLAGS}, "
                f"got {self.missingness_flag!r}"
            )


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def make_record(variant, modality, source_type, flag, context=None, **value) -> EvidenceRecord:
    """Construct an :class:`EvidenceRecord` from a :class:`~regvar.layer1.variants.Variant`."""
    return EvidenceRecord(
        chrom=variant.chrom,
        pos=variant.pos,
        ref=variant.ref,
        alt=variant.alt,
        modality=modality,
        source_type=source_type,
        evidence_value=ev(**value),
        context=context,
        missingness_flag=flag,
    )


# ---------------------------------------------------------------------------
# DataFrame helpers
# ---------------------------------------------------------------------------


def records_to_dataframe(records: List[EvidenceRecord]) -> pd.DataFrame:
    """Convert a list of :class:`EvidenceRecord` to a tidy DataFrame."""
    if not records:
        return pd.DataFrame(columns=LAYER1_COLUMNS)
    return pd.DataFrame([asdict(r) for r in records])[LAYER1_COLUMNS]


def validate_layer1_output(df: pd.DataFrame) -> None:
    """Raise loudly if the evidence DataFrame is malformed.

    Checks required columns, vocabulary constraints, and that every
    ``evidence_value`` cell is valid JSON.
    """
    missing_cols = [c for c in LAYER1_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Layer 1 output missing required column(s): {missing_cols}")
    for col, allowed in (
        ("source_type", SOURCE_TYPES),
        ("missingness_flag", MISSINGNESS_FLAGS),
        ("modality", MODALITIES),
    ):
        bad = set(df[col].dropna()) - set(allowed)
        if bad:
            raise ValueError(f"Layer 1 output has invalid {col} value(s): {bad}")
    for cell in df["evidence_value"]:
        json.loads(cell)  # raises on malformed JSON


def validate_matrix_complete(
    evidence_df: pd.DataFrame,
    variants,
    modalities=MODALITIES,
) -> None:
    """Assert every candidate variant has >= 1 record for every modality.

    Layer 2 depends on a complete matrix -- a missing cell means a bug in the
    collection step, not a legitimate absence.
    """
    from regvar.layer1.variants import canon_key  # local import to avoid circular

    expected = {v.key for v in variants}
    problems = []
    for mod in modalities:
        sub = evidence_df[evidence_df["modality"] == mod]
        got = {
            canon_key(r.chrom, r.pos, r.ref, r.alt)
            for r in sub[["chrom", "pos", "ref", "alt"]].itertuples(index=False)
        }
        missing = expected - got
        if missing:
            problems.append(f"{mod}: {len(missing)} candidate variant(s) have no record")
    if problems:
        raise ValueError("Evidence matrix incomplete -- " + "; ".join(problems))
