"""HBB sits in the beta-globin gene cluster alongside HBD and the HBBP1 pseudogene. With a wide
flank, VEP's region query can legitimately return transcript consequences belonging to those
neighboring genes -- an unrestricted 'union across all transcripts' would silently attribute a
neighboring gene's consequence to what's supposed to be a gene-specific regulatory truth-set.
This test locks in the gene_id restriction that prevents that."""
from regvar.layer0.vep import parse_vep_results
from tests.conftest import make_vep_result

TARGET_GENE = "ENSG00000244734"       # HBB
NEIGHBOR_GENE = "ENSG00000188536"     # HBA-like neighbor, standing in for HBD/HBBP1
CANONICAL_TX = "ENST00000335295"


def test_neighboring_gene_terms_are_excluded():
    vep_results = [make_vep_result(
        chrom="11", pos=5230000, ref="A", alt="T",
        transcript_consequences=[
            {"transcript_id": CANONICAL_TX, "gene_id": TARGET_GENE,
             "consequence_terms": ["intron_variant"]},
            {"transcript_id": "ENST00000999999", "gene_id": NEIGHBOR_GENE,
             "consequence_terms": ["3_prime_UTR_variant"]},   # belongs to the neighbor, not HBB
        ],
    )]

    df = parse_vep_results(vep_results, CANONICAL_TX, TARGET_GENE)
    row = df.iloc[0]

    assert "3_prime_UTR_variant" not in row["vep_all_transcript_terms"], (
        "A neighboring gene's consequence term leaked into the target gene's regulatory evidence"
    )
    assert "intron_variant" in row["vep_all_transcript_terms"]


def test_regulatory_feature_terms_remain_gene_unrestricted():
    """Regulatory Build features aren't transcript/gene-scoped, so they should NOT be filtered
    by gene_id the way transcript_consequences are."""
    vep_results = [make_vep_result(
        chrom="11", pos=5230000, ref="A", alt="T",
        transcript_consequences=[],
        regulatory_feature_consequences=[
            {"consequence_terms": ["regulatory_region_variant"]},
        ],
    )]

    df = parse_vep_results(vep_results, CANONICAL_TX, TARGET_GENE)
    row = df.iloc[0]

    assert "regulatory_region_variant" in row["vep_regulatory_terms"]
