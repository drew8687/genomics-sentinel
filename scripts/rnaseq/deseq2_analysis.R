#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════
#  GenomicsSentinel — Differential Expression Analysis
#  Tools: DESeq2, apeglm shrinkage, clusterProfiler GSEA
#  Output: DE table, MA plot, volcano, GSEA enrichment
# ═══════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(yaml)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Config ────────────────────────────────────────────────
cfg        <- yaml::read_yaml(snakemake@config[["configfile"]])
count_mat  <- snakemake@input[["counts"]]
coldata_f  <- snakemake@input[["coldata"]]
out_de     <- snakemake@output[["de_results"]]
out_ma     <- snakemake@output[["ma_plot"]]
out_vol    <- snakemake@output[["volcano"]]
out_gsea   <- snakemake@output[["gsea"]]
out_heat   <- snakemake@output[["heatmap"]]

padj_thr   <- cfg$deseq2$padj_threshold
lfc_thr    <- cfg$deseq2$lfc_threshold

cat("═══ DESeq2 Analysis — GenomicsSentinel ═══\n")

# ── Load data ─────────────────────────────────────────────
counts  <- read.csv(count_mat, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_f, row.names = 1)
coldata$condition <- factor(coldata$condition, levels = c("control", "treated"))

stopifnot(all(colnames(counts) == rownames(coldata)))

# ── DESeq2 object ─────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~ condition
)

# Pre-filter: keep genes with ≥10 reads in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
cat(sprintf("  Genes after filtering: %d\n", nrow(dds)))

# ── Run DESeq2 ────────────────────────────────────────────
dds <- DESeq(dds, parallel = TRUE)

# ── apeglm Log2FC shrinkage ───────────────────────────────
res_raw <- results(dds, contrast = c("condition", "treated", "control"),
                   alpha = padj_thr)
res     <- lfcShrink(dds, coef = "condition_treated_vs_control",
                     type = "apeglm", quiet = TRUE)

cat(sprintf("  Significant DE genes (padj<%.2f, |lfc|>%.1f): %d\n",
            padj_thr, lfc_thr,
            sum(res$padj < padj_thr & abs(res$log2FoldChange) > lfc_thr, na.rm = TRUE)))

# ── Export DE table ───────────────────────────────────────
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  mutate(
    significance = case_when(
      padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up",
      padj < padj_thr & log2FoldChange < -lfc_thr ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

write.csv(res_df, out_de, row.names = FALSE)

# ── MA Plot ───────────────────────────────────────────────
top_genes <- res_df %>% filter(significance != "NS") %>% head(20)

ma <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange,
                          color = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c("Up" = "#D85A30", "Down" = "#378ADD", "NS" = "#B4B2A9")) +
  geom_text_repel(data = top_genes, aes(label = gene_id),
                  size = 2.5, max.overlaps = 15) +
  labs(title = "MA Plot — DESeq2 (apeglm shrinkage)",
       x = "Mean normalized counts (log10)",
       y = "Log2 Fold Change",
       color = "DE status") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(out_ma, ma, width = 8, height = 5, dpi = 300)

# ── Volcano Plot ──────────────────────────────────────────
volcano <- ggplot(res_df, aes(x = log2FoldChange,
                               y = -log10(padj + 1e-300),
                               color = significance)) +
  geom_point(alpha = 0.6, size = 0.9) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  scale_color_manual(values = c("Up" = "#D85A30", "Down" = "#378ADD", "NS" = "#B4B2A9")) +
  geom_text_repel(data = top_genes, aes(label = gene_id),
                  size = 2.5, max.overlaps = 20) +
  labs(title = "Volcano Plot — Treated vs Control",
       x = "Log2 Fold Change", y = "-log10(padj)",
       color = "DE status") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(out_vol, volcano, width = 8, height = 6, dpi = 300)

# ── GSEA with clusterProfiler ─────────────────────────────
sig_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(entrez = mapIds(org.Hs.eg.db, gene_id,
                         "ENTREZID", "SYMBOL", multiVals = "first")) %>%
  filter(!is.na(entrez))

ranked <- setNames(sig_genes$log2FoldChange, sig_genes$entrez)
ranked <- sort(ranked, decreasing = TRUE)

gsea_bp <- gseGO(
  geneList     = ranked,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  seed         = 42
)

gsea_kegg <- gseKEGG(
  geneList     = ranked,
  organism     = "hsa",
  minGSSize    = 10,
  pvalueCutoff = 0.05,
  verbose      = FALSE,
  seed         = 42
)

gsea_combined <- rbind(
  as.data.frame(gsea_bp)   %>% mutate(ontology = "GO:BP"),
  as.data.frame(gsea_kegg) %>% mutate(ontology = "KEGG")
)
write.csv(gsea_combined, out_gsea, row.names = FALSE)

# ── Heatmap of top 50 DE genes ───────────────────────────
top50 <- res_df %>%
  filter(significance != "NS") %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene_id)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top50, ]
mat <- t(scale(t(mat)))

annotation_col <- data.frame(
  condition = coldata$condition,
  row.names = rownames(coldata)
)

ann_colors <- list(
  condition = c(control = "#378ADD", treated = "#D85A30")
)

png(out_heat, width = 10, height = 12, units = "in", res = 300)
pheatmap(mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize = 8,
  main = "Top 50 DE Genes — VST normalized"
)
dev.off()

cat("✓ DESeq2 + GSEA analysis complete.\n")
