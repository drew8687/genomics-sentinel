# 🧬 GenomicsSentinel

> **Production-grade multi-omics bioinformatics pipeline**  
> Reproducible · Containerized · ML-powered · One command to run

[![CI](https://github.com/driss-eloifi/genomics-sentinel/actions/workflows/ci.yml/badge.svg)](https://github.com/driss-eloifi/genomics-sentinel/actions)
[![Docker](https://img.shields.io/badge/docker-ready-2496ED?logo=docker)](https://github.com/driss-eloifi/genomics-sentinel/pkgs/container/pipeline)
[![Snakemake](https://img.shields.io/badge/workflow-Snakemake-009639?logo=snakemake)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

---

## Overview

**GenomicsSentinel** is a complete, reproducible multi-omics pipeline that integrates:

| Module | Tools | Output |
|--------|-------|--------|
| Quality Control | FastQC · Trimmomatic · MultiQC | QC reports, trimmed reads |
| Alignment (RNA) | STAR · featureCounts | BAM files, count matrix |
| Differential Expression | DESeq2 · apeglm · clusterProfiler | DE table, volcano, GSEA |
| scRNA-Seq | Seurat v5 · Harmony · SingleR | UMAP, cell types, markers |
| Variant Calling | GATK4 · BCFtools · SnpEff | Annotated VCF |
| Machine Learning | Random Forest · XGBoost · SHAP · Optuna | Biomarker model, ROC, SHAP |
| Reproducibility | Snakemake · Docker Compose · GitHub CI | Fully portable |

---

## Quick Start

### Prerequisites
- Docker Desktop ≥ 24.0
- Docker Compose v2
- 16 GB RAM, 8 CPU cores recommended

### 1. Clone & configure

```bash
git clone https://github.com/driss-eloifi/genomics-sentinel.git
cd genomics-sentinel
cp config/config.yaml.example config/config.yaml
# Edit config/config.yaml with your sample names and genome paths
```

### 2. Launch all services

```bash
docker compose up -d rstudio jupyter dashboard pgdb
```

| Service | URL | Credentials |
|---------|-----|-------------|
| RStudio Server | http://localhost:8787 | user: `rstudio` / pass: `genomics2024` |
| JupyterLab | http://localhost:8888 | token: `genomics2024` |
| Dashboard | http://localhost:3000 | — |
| DB Admin (Adminer) | http://localhost:8080 | — |

### 3. Run the full pipeline

```bash
docker compose run --rm pipeline \
  snakemake --cores all --use-conda --latency-wait 30
```

Or step-by-step:

```bash
# QC only
snakemake results/qc/multiqc/multiqc_report.html --cores 8

# RNA-Seq DE
snakemake results/rnaseq/deseq2/DE_results.csv --cores 8

# ML + SHAP
snakemake results/ml/model_metrics.json --cores 8
```

---

## Architecture

```
genomics-sentinel/
├── Snakefile                     # Master orchestrator
├── docker-compose.yml            # All services
├── config/
│   └── config.yaml               # Single configuration file
├── docker/
│   ├── Dockerfile.pipeline       # STAR · GATK · Snakemake
│   ├── Dockerfile.jupyter        # Python ML stack
│   └── envs/pipeline.yml         # Conda environment
├── snakemake/rules/
│   ├── qc.smk                    # FastQC → MultiQC
│   ├── alignment.smk             # STAR + Picard + GATK
│   ├── rnaseq.smk                # featureCounts + DESeq2
│   ├── scrnaseq.smk              # Seurat v5
│   ├── variants.smk              # HaplotypeCaller → SnpEff
│   └── ml.smk                    # XGBoost + SHAP + Optuna
├── scripts/
│   ├── rnaseq/deseq2_analysis.R  # Full DE + GSEA
│   ├── scrnaseq/seurat_analysis.R# scRNA-Seq clustering
│   └── ml/biomarker_ml.py        # ML + SHAP + MLflow
├── notebooks/                    # Exploratory notebooks
├── tests/                        # Unit tests (pytest + testthat)
└── .github/workflows/ci.yml      # CI: lint → test → build Docker
```

---

## Pipeline DAG

The Snakemake DAG shows all rule dependencies. Generate it with:

```bash
snakemake --dag | dot -Tsvg > docs/dag.svg
```

---

## Outputs

```
results/
├── qc/multiqc/multiqc_report.html      # Aggregate QC
├── alignment/*.sorted.bam               # Aligned reads
├── rnaseq/
│   ├── deseq2/DE_results.csv           # DE genes table
│   ├── deseq2/volcano_plot.png         # Volcano plot
│   ├── deseq2/MA_plot.png              # MA plot
│   └── gsea/enrichment_results.csv     # GO / KEGG GSEA
├── scrnaseq/seurat/
│   ├── umap_clusters.png               # UMAP by cluster
│   ├── umap_celltypes.png              # UMAP by cell type
│   └── marker_genes.csv               # Marker genes per cluster
├── variants/
│   ├── filtered.vcf.gz                 # Filtered variants
│   └── annotation/annotated.vcf        # SnpEff annotation
└── ml/
    ├── model_metrics.json              # AUC, MCC all models
    ├── shap_summary.png                # SHAP feature importance
    └── roc_curves.png                  # ROC comparison
```

---

## Customization

### Adding a new sample

```yaml
# config/config.yaml
samples:
  - SRR_myNewSample   # add here
```

### Changing reference genome

```yaml
genome:
  fasta: "resources/genome/GRCm39.fa"    # mouse, for example
  gtf:   "resources/genome/GRCm39.gtf"
  star_index: "resources/genome/STAR_index_mm39"
```

### Running only one module

```bash
snakemake results/scrnaseq/seurat/umap_clusters.png --cores 16
```

---

## Scientific Background

This pipeline was designed by **Driss El Oifi**, mathematician and bioinformatician, combining:

- 11+ years of teaching mathematical reasoning (hypothesis testing, statistical modeling)
- Research in **Q-learning applied to protein folding** (HP model, Markov decision processes)
- Bioinformatics pipelines in R (DESeq2/RNA-Seq, GEO datasets) and Python (Random Forest + SHAP, XGBoost microbiome)

The ML module implements the approach from:  
*El Oifi, D. (2024). Reinforcement learning strategies for 3D protein structure prediction using the HP model. Internal manuscript.*

---

## Citation

```bibtex
@software{eloifi2024genomicssentinel,
  author  = {El Oifi, Driss},
  title   = {GenomicsSentinel: A Reproducible Multi-Omics Pipeline},
  year    = {2024},
  url     = {https://github.com/driss-eloifi/genomics-sentinel}
}
```

---

## License

MIT — See [LICENSE](LICENSE)
