import pandas as pd
import pytest


@pytest.fixture
def hbb_like_gene_meta():
    """Minimal-strand gene metadata shaped like the real HBB lookup, without hitting Ensembl."""
    return {
        "gene_symbol": "HBB",
        "gene_id": "ENSG00000244734",
        "chrom": "11",
        "start": 5225464,
        "end": 5229395,
        "strand": -1,
    }


@pytest.fixture
def synthetic_variants_df():
    """A handful of synthetic ClinVar-shaped variants: one inside a distal 'LCR-like' interval,
    one elsewhere, matching the columns extract_clinvar_labels would normally produce."""
    return pd.DataFrame([
        {"chrom": "11", "pos": 5280000, "ref": "A", "alt": "G",  # inside the synthetic LCR interval
         "CLNSIG": "Pathogenic", "CLNREVSTAT": "criteria_provided,_single_submitter",
         "review_stars": 1, "label": "pathogenic"},
        {"chrom": "11", "pos": 5225000, "ref": "C", "alt": "T",  # not in any interval
         "CLNSIG": "Benign", "CLNREVSTAT": "criteria_provided,_single_submitter",
         "review_stars": 1, "label": "benign"},
    ])


@pytest.fixture
def synthetic_reg_intervals_df():
    """A distal enhancer-like interval (standing in for the HBB LCR) plus a promoter window,
    matching the Chromosome/Start/End/ccre_label shape flag_interval_overlap expects."""
    return pd.DataFrame([
        {"Chromosome": "11", "Start": 5269925, "End": 5304186, "ccre_label": "dELS"},
        {"Chromosome": "11", "Start": 5229395, "End": 5231895, "ccre_label": "computed_promoter_window"},
    ])


def make_vep_result(chrom, pos, ref, alt, transcript_consequences=None, regulatory_feature_consequences=None):
    """Build a minimal synthetic Ensembl VEP REST response for one variant."""
    return {
        "seq_region_name": chrom,
        "start": pos,
        "allele_string": f"{ref}/{alt}",
        "transcript_consequences": transcript_consequences or [],
        "regulatory_feature_consequences": regulatory_feature_consequences or [],
    }
