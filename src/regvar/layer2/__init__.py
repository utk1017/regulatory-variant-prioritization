"""Layer 2: Gap & Missingness Detection.

Determines which absent modalities trigger computational completion models
(AlphaGenome) vs. experimental negative results vs. transient errors, constructs
an agreement benchmark queue, and assigns variants into evidence-availability tiers.
"""
from regvar.layer2.pipeline import run_layer2
from regvar.layer2.reasons import (
    aggregate_modality_status,
    classify_evidence_reason,
    filter_records_to_targets,
    parse_evidence_value,
)
from regvar.layer2.schema import (
    AG_OUTPUT_TYPES,
    EVIDENCE_TIERS,
    LAYER2_AGREEMENT_COLUMNS,
    LAYER2_AVAILABILITY_COLUMNS,
    LAYER2_TRIGGER_COLUMNS,
    MODALITY_META,
    MODALITY_STATUSES,
    NEVER_COMPLETE_CATEGORIES,
    REASON_CATEGORIES,
    validate_layer1_input,
    validate_layer2_agreement_output,
    validate_layer2_availability_output,
    validate_layer2_triggers_output,
)
from regvar.layer2.tiers import (
    compute_evidence_tiers,
    label_coverage_warnings,
    tier_balance_warnings,
)
from regvar.layer2.triggers import (
    compile_agreement_check_queue,
    compile_alphagenome_triggers,
)

__all__ = [
    "run_layer2",
    "LAYER2_AVAILABILITY_COLUMNS",
    "LAYER2_TRIGGER_COLUMNS",
    "LAYER2_AGREEMENT_COLUMNS",
    "EVIDENCE_TIERS",
    "MODALITY_STATUSES",
    "AG_OUTPUT_TYPES",
    "MODALITY_META",
    "REASON_CATEGORIES",
    "NEVER_COMPLETE_CATEGORIES",
    "validate_layer2_availability_output",
    "validate_layer2_triggers_output",
    "validate_layer2_agreement_output",
    "validate_layer1_input",
    "classify_evidence_reason",
    "aggregate_modality_status",
    "filter_records_to_targets",
    "parse_evidence_value",
    "compute_evidence_tiers",
    "tier_balance_warnings",
    "label_coverage_warnings",
    "compile_alphagenome_triggers",
    "compile_agreement_check_queue",
]
