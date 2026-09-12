import os
import pandas as pd
import pyranges as pr
from regvar.common.utils import _run, to_ucsc_chrom, to_plain_chrom

def fetch_bigbedtobed_tool():
    tool_path = "scripts/bigBedToBed"
    if not os.path.exists(tool_path):
        _run(f"wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed -O {tool_path}")
        _run(f"chmod +x {tool_path}")
    return f"./{tool_path}"

def fetch_ccres(chrom, start, end, out_bed="data/cache/ccre_region.bed"):
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    tool = fetch_bigbedtobed_tool()
    remote_bb = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"
    cmd = f"{tool} {remote_bb} -chrom={to_ucsc_chrom(chrom)} -start={start} -end={end} {out_bed}"
    r = _run(cmd, check=False)
    if r.returncode != 0 or not os.path.exists(out_bed) or os.path.getsize(out_bed) == 0:
        print("WARNING: cCRE fetch failed or returned empty for this region")
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "ccre_label"])
    
    cols = ["Chromosome", "Start", "End", "name", "score", "strand", "thickStart", "thickEnd",
            "reserved", "ccre_group", "ccre_group2", "zscore", "ucscLabel", "accession", "ccre_label"]
    df = pd.read_csv(out_bed, sep="\t", header=None, names=cols)
    df["Chromosome"] = to_plain_chrom(chrom)
    return df[["Chromosome", "Start", "End", "ccre_label"]]

def flag_interval_overlap(variants_df, reg_intervals_df):
    variants_df = variants_df.reset_index(drop=True).copy()
    variants_df["_row_id"] = variants_df.index

    if reg_intervals_df.empty:
        variants_df["interval_evidence"] = [[] for _ in range(len(variants_df))]
        return variants_df.drop(columns=["_row_id"])

    var_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": variants_df["chrom"],
        "Start": variants_df["pos"] - 1,
        "End": variants_df["pos"],
        "row_id": variants_df["_row_id"],
    }))
    reg_pr = pr.PyRanges(reg_intervals_df)

    joined = var_pr.join(reg_pr).df
    if joined.empty or "row_id" not in joined.columns:
        variants_df["interval_evidence"] = [[] for _ in range(len(variants_df))]
        return variants_df.drop(columns=["_row_id"])

    evidence_map = joined.groupby("row_id")["ccre_label"].apply(list).to_dict()
    variants_df["interval_evidence"] = variants_df["_row_id"].map(lambda i: evidence_map.get(i, []))
    return variants_df.drop(columns=["_row_id"])
