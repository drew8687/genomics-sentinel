"""
GenomicsSentinel — Multi-Omics Reproducible Pipeline
Author: Driss El Oifi
Description: End-to-end pipeline: QC → Alignment → RNA-Seq DE → scRNA-Seq → Variant Calling → ML
"""

configfile: "config/config.yaml"

include: "snakemake/rules/qc.smk"
include: "snakemake/rules/alignment.smk"
include: "snakemake/rules/rnaseq.smk"
include: "snakemake/rules/scrnaseq.smk"
include: "snakemake/rules/variants.smk"
include: "snakemake/rules/ml.smk"
include: "snakemake/rules/report.smk"

samples = config["samples"]
conditions = config["conditions"]

rule all:
    input:
        # QC
        expand("results/qc/fastqc/{sample}_R1_fastqc.html", sample=samples),
        expand("results/qc/multiqc/multiqc_report.html"),
        # Alignment
        expand("results/alignment/{sample}.sorted.bam.bai", sample=samples),
        # RNA-Seq DE
        "results/rnaseq/deseq2/DE_results.csv",
        "results/rnaseq/deseq2/MA_plot.png",
        "results/rnaseq/deseq2/volcano_plot.png",
        "results/rnaseq/gsea/enrichment_results.csv",
        # scRNA-Seq
        "results/scrnaseq/seurat/umap_clusters.png",
        "results/scrnaseq/seurat/marker_genes.csv",
        # Variant Calling
        "results/variants/filtered.vcf.gz",
        "results/variants/annotation/annotated.vcf",
        # ML
        "results/ml/model_metrics.json",
        "results/ml/shap_summary.png",
        # Final report
        "results/report/genomics_sentinel_report.html"
