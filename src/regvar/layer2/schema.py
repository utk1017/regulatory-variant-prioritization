"""Layer 2 Schema and Data Contracts.

Defines output schemas, modality metadata, reason categorization, and
strict validation routines for Layer 2 gap detection, AlphaGenome triggers,
agreement benchmark checks, and stratification tiers.
"""
from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Sequence, Set

import pandas as pd

from regvar.common.schema import MISSINGNESS_FLAGS

LAYER2_AVAILABILITY_COLUMNS: List[str] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "gwas_status",
    "eqtl_status",
    "mpra_status",
    "crispr_status",
    "gnomad_af_status",
    "n_positive",
    "n_null",
    "n_absent",
    "n_tested_experimental",
    "evidence_tier",
    "primary_gap_reasons",
]

LAYER2_TRIGGER_COLUMNS: List[str] = [
    "trigger_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "modality",
    "context",
    "target_gene",
    "alphagenome_output_types",
    "alphagenome_proxy_quality",
    "alphagenome_binding_status",
    "receptive_field_bp",
    "trigger_reason",
    "source_evidence_value",
]

LAYER2_AGREEMENT_COLUMNS: List[str] = [
    "check_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "modality",
    "context",
    "observed_flag",
    "observed_value",
    "validation_purpose",
]

EVIDENCE_TIERS: Sequence[str] = (
    "Tier 1 (Evidence-Rich)",
    "Tier 2 (Moderate Evidence)",
    "Tier 3 (Evidence-Sparse)",
)

MODALITY_STATUSES: Sequence[str] = MISSINGNESS_FLAGS

AG_OUTPUT_TYPES: Sequence[str] = (
    "gene_expression",
    "transcription_initiation",
    "chromatin_accessibility",
    "histone_modification",
    "tf_binding",
    "contact_map",
    "splicing",
)

MODALITY_META: Dict[str, Dict[str, Any]] = {
    "gwas": {
        "kind": "association",
        "completable": False,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
        "description": "GWAS Catalog association / fine-mapping",
    },
    "eqtl": {
        "kind": "expression_effect",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
        "description": "GTEx expression QTL in target tissues",
        "alphagenome_output_types": ["gene_expression"],
        "proxy_quality": "direct",
        "receptive_field_bp": 1048576,
    },
    "mpra": {
        "kind": "reporter_activity",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
        "description": "Massively Parallel Reporter Assay allelic skew",
        "alphagenome_output_types": ["chromatin_accessibility", "transcription_initiation"],
        "proxy_quality": "indirect_proxy",
        "receptive_field_bp": 2048,
    },
    "crispr": {
        "kind": "perturbation_effect",
        "completable": True,
        "use_in_fusion": True,
        "direction": "supports_regulatory_relevance",
        "description": "CRISPRi / flowFISH perturbation screen",
        "alphagenome_output_types": ["gene_expression"],
        "proxy_quality": "indirect_proxy",
        "receptive_field_bp": 1048576,
    },
    "gnomad_af": {
        "kind": "population_frequency",
        "completable": False,
        "use_in_fusion": False,
        "direction": "benign_leaning_when_observed",
        "description": "gnomAD genome allele frequency (annotation only)",
    },
}

REASON_CATEGORIES: Dict[str, str] = {
    # 1. True biological or dataset coverage gaps -> Trigger AlphaGenome
    "not_testable_rare_or_unknown_frequency": "GAP_BIOLOGICAL",
    "no_significant_eqtl_association_returned_test_status_unverified": "GAP_BIOLOGICAL",
    "not_tested_in_provided_table": "GAP_BIOLOGICAL",
    "not_in_catalog": "GAP_BIOLOGICAL",
    "not_observed_in_genome_callset": "GAP_BIOLOGICAL",
    "no_data_source_configured": "GAP_BIOLOGICAL",
    "no_usable_association_records": "GAP_BIOLOGICAL",
    "p_value_missing": "GAP_BIOLOGICAL",

    # 2. Tested null -> Experimental negative measurement (NEVER trigger AlphaGenome)
    "common_variant_not_significant": "OBSERVED_NULL",
    "tested_not_significant": "OBSERVED_NULL",
    "subthreshold_association": "OBSERVED_NULL",

    # 3. Technical / query failures -> Actionable transient warnings (NEVER trigger AlphaGenome)
    "query_failed": "ERROR_TRANSIENT",
    "incomplete_region_query": "ERROR_TRANSIENT",
    "incomplete_tissue_query": "ERROR_TRANSIENT",
}

NEVER_COMPLETE_CATEGORIES: Set[str] = {"ERROR_TRANSIENT", "UNKNOWN_REASON"}


def validate_layer2_availability_output(df: pd.DataFrame) -> None:
    """Validate completeness and allowable values of the Modality Availability table."""
    missing = [c for c in LAYER2_AVAILABILITY_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Layer 2 availability DataFrame missing required column(s): {missing}")

    for col in ["gwas_status", "eqtl_status", "mpra_status", "crispr_status", "gnomad_af_status"]:
        invalid = set(df[col].dropna()) - set(MODALITY_STATUSES)
        if invalid:
            raise ValueError(f"Invalid status value(s) in {col}: {invalid}")

    invalid_tiers = set(df["evidence_tier"].dropna()) - set(EVIDENCE_TIERS)
    if invalid_tiers:
        raise ValueError(f"Invalid evidence tier value(s): {invalid_tiers}")


def validate_layer2_triggers_output(
    df: pd.DataFrame, completable_modalities: Optional[Sequence[str]] = None
) -> None:
    """Validate schema and allowable types/proxies of the AlphaGenome triggers table."""
    missing = [c for c in LAYER2_TRIGGER_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Layer 2 triggers DataFrame missing required column(s): {missing}")

    if completable_modalities is None:
        completable_modalities = [m for m, meta in MODALITY_META.items() if meta.get("completable")]

    if not df.empty:
        non_completable = set(df["modality"]) - set(completable_modalities)
        if non_completable:
            raise ValueError(f"AlphaGenome triggers contain non-completable modality: {non_completable}")

        empty_types = df[df["alphagenome_output_types"].astype(str).str.strip() == ""]
        if not empty_types.empty:
            bad_mods = sorted(set(empty_types["modality"]))
            raise ValueError(
                f"AlphaGenome triggers for modality {bad_mods} have no alphagenome_output_types "
                f"-- every completable modality in MODALITY_META must declare at least one."
            )

        bad_types = {t for s in df["alphagenome_output_types"] for t in str(s).split(";")} - set(AG_OUTPUT_TYPES)
        if bad_types:
            raise ValueError(f"AlphaGenome triggers contain unknown output type(s): {bad_types}")

        for s in df["source_evidence_value"]:
            val = json.loads(s) if isinstance(s, str) else s
            if not isinstance(val, dict):
                raise ValueError("source_evidence_value must serialize a JSON object (double-encoded?)")


def validate_layer2_agreement_output(df: pd.DataFrame) -> None:
    """Validate schema and observed flags of the internal agreement benchmark table."""
    missing = [c for c in LAYER2_AGREEMENT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Layer 2 agreement checks DataFrame missing required column(s): {missing}")

    if not df.empty:
        invalid = set(df["observed_flag"]) - {"observed_positive", "observed_null"}
        if invalid:
            raise ValueError(f"Agreement check contains non-observed flag: {invalid}")


def validate_layer1_input(df: pd.DataFrame) -> None:
    """Validate the Layer 1 evidence table BEFORE Layer 2 processes it.

    Ensures required columns exist, no missingness_flags are NaN, and all flags
    match valid enum values to prevent silent gap misclassification.
    """
    required = ["chrom", "pos", "ref", "alt", "modality", "missingness_flag"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Layer 1 evidence is missing required column(s): {missing} "
            f"(found: {sorted(df.columns)}). Was this produced by the Layer 1 pipeline?"
        )

    if df["missingness_flag"].isna().any():
        n_na = int(df["missingness_flag"].isna().sum())
        raise ValueError(
            f"Layer 1 evidence has {n_na} row(s) with an empty missingness_flag; "
            f"refusing to guess whether they are gaps."
        )

    bad_flags = set(df["missingness_flag"]) - set(MODALITY_STATUSES)
    if bad_flags:
        raise ValueError(
            f"Layer 1 evidence has invalid missingness_flag value(s): {sorted(bad_flags)} "
            f"(allowed: {list(MODALITY_STATUSES)}). Refusing to treat them as gaps."
        )

    unknown_mods = sorted(set(df["modality"].dropna()) - set(MODALITY_META))
    if unknown_mods:
        print(
            f"NOTE: Layer 1 contains modality/modalities not in MODALITY_META: {unknown_mods}. "
            f"Layer 2 ignores them (not counted toward tiers, never completed)."
        )
