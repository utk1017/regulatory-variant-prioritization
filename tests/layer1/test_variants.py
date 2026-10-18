import pandas as pd
from regvar.layer1.variants import canon_key, plan_windows, to_variants
from regvar.layer1.gtex import parse_gtex_variant_id


def test_canon_key_normalization():
    assert canon_key("chr11", 100, "A", "G") == ("11", 100, "A", "G")
    assert canon_key("11", 100, "a", "g") == ("11", 100, "A", "G")
    # Non-minimal indel normalization (shared suffix and prefix)
    assert canon_key("11", 5227100, "G", "GA") == canon_key("11", 5227099, "AG", "AGA")
    assert canon_key("chrX", 1000, "ACGT", "AG") == ("X", 1001, "CGT", "G")


def test_plan_windows():
    positions = [100, 200, 300, 30000, 30500, 80000]
    windows = plan_windows(positions, max_span=25000)
    assert windows == [(100, 300), (30000, 30500), (80000, 80000)]
    assert plan_windows([], max_span=25000) == []
    assert plan_windows([500], max_span=25000) == [(500, 500)]


def test_to_variants():
    df = pd.DataFrame([
        {"chrom": "chr11", "pos": 100, "ref": "A", "alt": "G"},
        {"chrom": "11", "pos": 200, "ref": "C", "alt": "T"},
    ])
    variants = to_variants(df)
    assert len(variants) == 2
    assert variants[0].chrom == "11"
    assert variants[0].pos == 100
    assert variants[0].key == ("11", 100, "A", "G")


def test_parse_gtex_variant_id():
    assert parse_gtex_variant_id("chr11_5227002_T_C_b38") == ("11", 5227002, "T", "C")
    assert parse_gtex_variant_id("chr1_55516888_G_GA_b38") == ("1", 55516888, "G", "GA")
    assert parse_gtex_variant_id("malformed_id") is None
    assert parse_gtex_variant_id(None) is None
