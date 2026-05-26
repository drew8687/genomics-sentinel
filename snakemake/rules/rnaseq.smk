# ──────────────────────────────────────────────────────────
#  rnaseq.smk — featureCounts + DESeq2 + clusterProfiler
# ──────────────────────────────────────────────────────────
rule feature_counts:
    input:
        bams = expand("results/alignment/{sample}.sorted.bam", sample=config["samples"]),
        gtf  = config["genome"]["gtf"]
    output:
        counts = "results/rnaseq/counts/raw_counts.txt",
        summary = "results/rnaseq/counts/raw_counts.txt.summary"
    threads: 8
    shell:
        """
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -p --countReadPairs \
            -s 0 \
            {input.bams}
        """

rule deseq2:
    input:
        counts  = "results/rnaseq/counts/raw_counts.txt",
        coldata = "config/sample_metadata.csv"
    output:
        de_results = "results/rnaseq/deseq2/DE_results.csv",
        ma_plot    = "results/rnaseq/deseq2/MA_plot.png",
        volcano    = "results/rnaseq/deseq2/volcano_plot.png",
        gsea       = "results/rnaseq/gsea/enrichment_results.csv",
        heatmap    = "results/rnaseq/deseq2/heatmap_top50.png"
    script:
        "../../scripts/rnaseq/deseq2_analysis.R"


# ──────────────────────────────────────────────────────────
#  scrnaseq.smk — Seurat v5
# ──────────────────────────────────────────────────────────
rule seurat_analysis:
    input:
        cellranger = config["scrnaseq"]["cellranger_out"]
    output:
        umap    = "results/scrnaseq/seurat/umap_clusters.png",
        markers = "results/scrnaseq/seurat/marker_genes.csv",
        seurat_h5 = "results/scrnaseq/seurat/seurat_object.h5seurat"
    script:
        "../../scripts/scrnaseq/seurat_analysis.R"


# ──────────────────────────────────────────────────────────
#  ml.smk — Feature matrix → RF + XGB + SVM + SHAP
# ──────────────────────────────────────────────────────────
rule build_feature_matrix:
    input:
        de      = "results/rnaseq/deseq2/DE_results.csv",
        counts  = "results/rnaseq/counts/raw_counts.txt",
        coldata = "config/sample_metadata.csv"
    output:
        features = "results/ml/feature_matrix.csv",
        labels   = "results/ml/labels.csv"
    run:
        import pandas as pd
        de     = pd.read_csv(input.de)
        sig    = de[de["significance"] != "NS"]["gene_id"].tolist()
        counts = pd.read_csv(input.counts, sep="\t", comment="#", index_col=0)
        counts = counts.iloc[:, 5:]
        counts.columns = [c.split("/")[-1].replace(".sorted.bam", "") for c in counts.columns]
        feat   = counts.loc[counts.index.isin(sig)].T
        feat.to_csv(output.features)
        meta   = pd.read_csv(input.coldata, index_col=0)
        meta[["condition"]].to_csv(output.labels)

rule ml_biomarker:
    input:
        features = "results/ml/feature_matrix.csv",
        labels   = "results/ml/labels.csv"
    output:
        metrics   = "results/ml/model_metrics.json",
        shap_plot = "results/ml/shap_summary.png",
        roc_plot  = "results/ml/roc_curves.png"
    script:
        "../../scripts/ml/biomarker_ml.py"


# ──────────────────────────────────────────────────────────
#  report.smk — Quarto HTML report
# ──────────────────────────────────────────────────────────
rule generate_report:
    input:
        de_results = "results/rnaseq/deseq2/DE_results.csv",
        gsea       = "results/rnaseq/gsea/enrichment_results.csv",
        ml_metrics = "results/ml/model_metrics.json",
        ma_plot    = "results/rnaseq/deseq2/MA_plot.png",
        volcano    = "results/rnaseq/deseq2/volcano_plot.png",
        shap_plot  = "results/ml/shap_summary.png",
        umap       = "results/scrnaseq/seurat/umap_clusters.png",
        multiqc    = "results/qc/multiqc/multiqc_report.html"
    output:
        report = "results/report/genomics_sentinel_report.html"
    shell:
        """
        quarto render docs/report_template.qmd \
            --to html \
            --output {output.report} \
            -P de_results:{input.de_results} \
            -P ml_metrics:{input.ml_metrics}
        """
