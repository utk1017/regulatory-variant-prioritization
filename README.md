# Regulatory Variant Prioritization Framework

A multimodal, provenance-aware evidence framework for prioritizing regulatory non-coding variants. 

This pipeline fuses GWAS, eQTL, MPRA, and CRISPR evidence with AlphaGenome-based computational completion for missing modalities. The results are calibrated against ClinVar reference sets using an ACMG/AMP-style Bayesian framework. 

## 📊 Project Status
- **Layer 0 (Variant Intake & Harmonization):** ✅ Finalized and tested. (Flags whether a variant sits in a regulatory region).
- **Layers 1–7 (Experimental & Computational Integration):** 🚧 In progress.

For a detailed breakdown of the architecture, see [docs/framework_overview.md](docs/framework_overview.md).

---

## ⚙️ Requirements

> [!IMPORTANT]
> **Operating System:** This pipeline relies heavily on bioinformatics binaries (`bcftools`, `samtools`, `tabix`) and UCSC command-line tools. It requires **Linux** or **macOS (Intel/Rosetta or via Conda)**. Windows users should use WSL2 (Windows Subsystem for Linux).

- Python >= 3.9
- Conda (Miniconda or Anaconda)

## 🚀 Installation

Setting up the environment is streamlined using Conda:

```bash
# 1. Create the conda environment (installs bioconda deps + pip installs the local package)
conda env create -f environment.yml
conda activate regvar-env

# 2. Setup external binaries (e.g., UCSC bigBedToBed)
bash scripts/setup_env.sh
```

## 💻 Quickstart

Run Layer 0 for the `HBB` locus using the provided configuration:

```bash
regvar run --layer 0 --config configs/hbb.yaml
```

The pipeline outputs to `data/outputs/HBB/layer0_regulatory_truthset.csv`. Alongside the data, a `manifest.json` is generated to record exact data provenance (ClinVar snapshot date, reference build, git commit, and config used).

## 🧪 Testing

The test suite runs against small synthetic fixtures and requires no network access. 

```bash
# The 'dev' dependencies are already installed via environment.yml
pytest tests/ -v
```

CI (`.github/workflows/tests.yml`) enforces these tests on every push. Every bug found during development must include a regression test.

## 🤝 Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) before submitting pull requests. We maintain strict design rules regarding schema stability, data provenance, and module scoping.

## 📄 License

Copyright (c) 2026 Utkarsh Gupta , Kamakshi Mudgal. All rights reserved. 
See the [LICENSE](LICENSE) file for complete details and restrictions.
