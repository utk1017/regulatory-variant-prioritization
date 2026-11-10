# Regulatory Variant Prioritization Framework

A multimodal, provenance-aware evidence framework for prioritizing regulatory non-coding variants. 

This pipeline fuses GWAS, eQTL, MPRA, and CRISPR evidence with AlphaGenome-based computational completion for missing modalities. The results are calibrated against ClinVar reference sets using an ACMG/AMP-style Bayesian framework. 

## 📊 Project Status
- **Layer 0 (Variant Intake & Harmonization):** ✅ Finalized and tested. (Flags whether a variant sits in a regulatory region).
- **Layer 1 (Multi-Source Evidence Collection):** ✅ Finalized and tested. (Standardizes evidence from gnomAD, GWAS, GTEx, MPRA, and CRISPR).
- **Layer 2 (Evidence Rules & Completability):** ✅ Finalized and tested. (Disambiguates gaps, compiles AlphaGenome trigger queue, isolates agreement benchmarks, and computes stratification tiers).
- **Layers 3–7 (Computational Completion & Bayesian Integration):** 🚧 In progress.

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

Run Layer 1 to collect multimodal evidence:

```bash
regvar run --layer 1 --config configs/hbb.yaml
```

Run Layer 2 to detect missingness gaps and compile AlphaGenome triggers:

```bash
regvar run --layer 2 --config configs/hbb.yaml
```

The pipeline outputs:
- `data/outputs/HBB/layer2_modality_availability.csv`: Per-variant availability vector and evidence tiers (Tiers 1–3).
- `data/outputs/HBB/layer2_alphagenome_triggers.csv`: Inference queue for missing completable modalities.
- `data/outputs/HBB/layer2_agreement_checks.csv`: Observed experimental measurements isolated for model benchmarking.
- `data/outputs/HBB/layer2_manifest.json`: Full provenance audit trail, trigger distributions, and pipeline warnings.

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
