# Regulatory Variant Prioritization Framework

A multimodal, provenance-aware evidence framework for prioritizing regulatory non-coding variants. 

Fuses GWAS, eQTL, MPRA, and CRISPR evidence with AlphaGenome-based computational completion for missing modalities, calibrated against ClinVar reference sets in an ACMG/AMP-style Bayesian framework. 

## Status
- **Layer 0** (variant intake & regulatory-candidate flagging): Finalized and tested.
- **Layers 1–7**: In progress.

## Installation

```bash
# 1. Create the conda environment for bioconda dependencies (bcftools, samtools)
conda env create -f environment.yml
conda activate regvar-env

# 2. Install the python package
pip install -e .

# 3. Setup external binaries (UCSC tools)
bash scripts/setup_env.sh
```

## Quickstart

Run Layer 0 for the HBB locus:
```bash
regvar run --layer 0 --config configs/hbb.yaml
```

## Contributing
Please see `CONTRIBUTING.md` for our strict design rules regarding schema stability, data provenance, and module scoping.
