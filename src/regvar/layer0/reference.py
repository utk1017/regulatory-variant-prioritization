import os
from regvar.common.utils import _run

def download_reference(chrom, out_dir="data/cache", out_prefix="ref"):
    os.makedirs(out_dir, exist_ok=True)
    fa_gz = os.path.join(out_dir, f"{out_prefix}_chr{chrom}.fa.gz")
    fa = os.path.join(out_dir, f"{out_prefix}_chr{chrom}.fa")
    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr{chrom}.fa.gz"

    _run(f"wget -nc {url} -O {fa_gz}")
    integrity = _run(f"gzip -t {fa_gz}", check=False)
    if integrity.returncode != 0:
        print(f"Corrupt download detected, re-fetching {fa_gz}")
        _run(f"rm -f {fa_gz}")
        _run(f"wget {url} -O {fa_gz}")
        _run(f"gzip -t {fa_gz}")

    _run(f"gunzip -kf {fa_gz}")
    _run(f"sed -i 's/^>chr{chrom}/>{chrom}/' {fa}")
    _run(f"samtools faidx {fa}")
    if not os.path.exists(f"{fa}.fai"):
        raise RuntimeError(f"samtools faidx did not produce an index for {fa}")
    return fa
