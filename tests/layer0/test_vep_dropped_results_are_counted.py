"""These two warnings existed in the finalized notebook and were dropped during the port to
src/regvar/layer0/vep.py -- this test exists so that regression can't happen silently again."""
from regvar.layer0.vep import parse_vep_results
from tests.conftest import make_vep_result

TARGET_GENE = "ENSG00000244734"
CANONICAL_TX = "ENST00000335295"


def test_missing_allele_string_is_counted_not_silent(capsys):
    good = make_vep_result(chrom="11", pos=5230000, ref="A", alt="T")
    bad = make_vep_result(chrom="11", pos=5230100, ref="A", alt="T")
    bad["allele_string"] = None  # simulate a malformed VEP response

    df = parse_vep_results([good, bad], CANONICAL_TX, TARGET_GENE)

    assert len(df) == 1, "Only the well-formed result should produce a row"
    captured = capsys.readouterr()
    assert "1" in captured.out and "allele_string" in captured.out, (
        "A dropped VEP result must be counted/printed, not silently discarded"
    )
