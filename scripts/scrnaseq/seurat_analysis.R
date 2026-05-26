#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════
#  GenomicsSentinel — scRNA-Seq Analysis with Seurat v5
#  Pipeline: Load CellRanger → QC → Normalize → Cluster →
#            UMAP → Marker genes → Cell type annotation
# ═══════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(yaml)
  library(harmony)      # batch correction
  library(SingleR)      # automated cell type annotation
  library(celldex)      # reference datasets for SingleR
})

cfg       <- yaml::read_yaml(snakemake@config[["configfile"]])
cr_out    <- snakemake@input[["cellranger"]]
out_umap  <- snakemake@output[["umap"]]
out_mark  <- snakemake@output[["markers"]]
out_h5    <- snakemake@output[["seurat_h5"]]

sc_cfg    <- cfg$scrnaseq

cat("═══ scRNA-Seq — Seurat v5 Analysis ═══\n")

# ── Load CellRanger output ────────────────────────────────
counts <- Read10X(data.dir = cr_out)
sobj   <- CreateSeuratObject(
  counts    = counts,
  min.cells = sc_cfg$min_cells,
  min.features = sc_cfg$min_features,
  project   = "GenomicsSentinel"
)

# ── QC metrics ────────────────────────────────────────────
sobj[["pct_mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
sobj[["pct_rb"]] <- PercentageFeatureSet(sobj, pattern = "^RPS|^RPL")

qc_plot <- VlnPlot(sobj,
  features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
  ncol = 3, pt.size = 0.1
) & theme(axis.title.x = element_blank())

ggsave(gsub("umap", "qc_violin", out_umap), qc_plot,
       width = 14, height = 5, dpi = 200)

# ── Filter cells ──────────────────────────────────────────
sobj <- subset(sobj,
  nFeature_RNA > sc_cfg$min_features &
  nFeature_RNA < 8000 &
  pct_mt < sc_cfg$max_mt_pct
)
cat(sprintf("  Cells after QC filter: %d\n", ncol(sobj)))

# ── Normalize — SCTransform v2 (replaces LogNormalize) ───
sobj <- SCTransform(sobj, vst.flavor = "v2", verbose = FALSE)

# ── PCA ───────────────────────────────────────────────────
sobj <- RunPCA(sobj, npcs = sc_cfg$n_pcs, verbose = FALSE)

# ── Batch correction with Harmony ─────────────────────────
# (if multiple samples; here as example)
# sobj <- RunHarmony(sobj, group.by.vars = "orig.ident")

# ── Clustering ────────────────────────────────────────────
sobj <- FindNeighbors(sobj, dims = 1:sc_cfg$n_pcs, verbose = FALSE)
sobj <- FindClusters(sobj, resolution = sc_cfg$resolution, verbose = FALSE)
sobj <- RunUMAP(sobj, dims = 1:sc_cfg$n_pcs, verbose = FALSE)

n_clusters <- length(levels(sobj$seurat_clusters))
cat(sprintf("  Clusters found: %d\n", n_clusters))

# ── UMAP plot ─────────────────────────────────────────────
umap_plot <- DimPlot(sobj, reduction = "umap",
                     label = TRUE, label.size = 4,
                     repel = TRUE, pt.size = 0.5) +
  labs(title = "UMAP — Seurat Clusters") +
  theme_minimal()

ggsave(out_umap, umap_plot, width = 9, height = 7, dpi = 300)

# ── Automated cell type annotation with SingleR ───────────
hpca.se <- HumanPrimaryCellAtlasData()
sce      <- as.SingleCellExperiment(sobj)
pred     <- SingleR(test = sce, ref = hpca.se,
                    assay.type.test = 1,
                    labels = hpca.se$label.main)

sobj$cell_type <- pred$labels

umap_ct <- DimPlot(sobj, reduction = "umap",
                   group.by = "cell_type",
                   label = TRUE, label.size = 3.5,
                   repel = TRUE, pt.size = 0.5) +
  labs(title = "UMAP — Cell type annotation (SingleR)") +
  theme_minimal()

ggsave(gsub("umap_clusters", "umap_celltypes", out_umap),
       umap_ct, width = 11, height = 7, dpi = 300)

# ── Marker genes per cluster ──────────────────────────────
Idents(sobj) <- "seurat_clusters"

markers <- FindAllMarkers(
  sobj,
  only.pos = TRUE,
  min.pct  = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox"
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(top_markers, out_mark, row.names = FALSE)

# ── Dot plot top markers ──────────────────────────────────
top5 <- top_markers %>% group_by(cluster) %>% head(5)

dot_plot <- DotPlot(sobj, features = unique(top5$gene)) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 7)) +
  labs(title = "Top 5 marker genes per cluster")

ggsave(gsub("umap_clusters", "dotplot_markers", out_umap),
       dot_plot, width = 14, height = 6, dpi = 300)

# ── Save Seurat object ────────────────────────────────────
SaveH5Seurat(sobj, filename = out_h5, overwrite = TRUE)

cat("✓ scRNA-Seq analysis complete.\n")
