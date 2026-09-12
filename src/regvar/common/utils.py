import subprocess

def _run(cmd, check=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        msg = f"Command failed (exit {r.returncode}): {cmd}\nSTDERR:\n{r.stderr}"
        if check:
            raise RuntimeError(msg)
        else:
            print(msg)
    return r

def to_ucsc_chrom(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith("chr") else f"chr{chrom}"

def to_plain_chrom(chrom):
    chrom = str(chrom)
    return chrom[3:] if chrom.startswith("chr") else chrom
