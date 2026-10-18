"""Layer 1: Multi-source evidence collection and harmonization.

Extracts experimental and computational evidence across gnomAD, GWAS Catalog,
GTEx Portal, MPRA, and CRISPR screens.
"""
from regvar.layer1.pipeline import run_layer1
from regvar.layer1.schema import (
    EvidenceRecord,
    MODALITIES,
    MODALITY_META,
    SOURCE_TYPES,
    MISSINGNESS_FLAGS,
    validate_layer1_output,
    validate_matrix_complete,
)
from regvar.layer1.variants import Variant, canon_key, plan_windows, to_variants

__all__ = [
    "run_layer1",
    "EvidenceRecord",
    "Variant",
    "canon_key",
    "plan_windows",
    "to_variants",
    "MODALITIES",
    "MODALITY_META",
    "SOURCE_TYPES",
    "MISSINGNESS_FLAGS",
    "validate_layer1_output",
    "validate_matrix_complete",
]
