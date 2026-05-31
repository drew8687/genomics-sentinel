# ================================================================
# Analyse d'expression différentielle — DESeq2
# Auteur  : Driss El Oifi
# Date    : mai 2026
# Dataset : GSE152418 — COVID-19 Severe vs Healthy (Cell, 2020)
# ================================================================
# Approche statistique :
#   Les données RNA-Seq suivent une loi binomiale négative
#   (variance > moyenne, contrairement à Poisson).
#   DESeq2 modélise cette sur-dispersion par gène.
#
#   Modèle : K_ij ~ NB(mu_ij, alpha_i)
#   avec mu_ij = s_j * q_ij  (size factor * expression vraie)
#
#   Test de Wald sur log2FC après shrinkage apeglm
#   Référence : Love et al., Genome Biology, 2014
# ================================================================

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

# ----------------------------------------------------------------
# 1. Chargement des données
# ----------------------------------------------------------------
# count_mat  : matrice de comptage bruts (entiers)
#              lignes = gènes, colonnes = échantillons
# coldata_f  : métadonnées (condition, sexe, etc.)
# ----------------------------------------------------------------
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

message("=== Analyse DESeq2 en cours ===")

# ----------------------------------------------------------------
# 2. Création de l'objet DESeq2
# ----------------------------------------------------------------
# Le design ~ condition indique à DESeq2 quelle variable
# utiliser pour le test. On peut ajouter des covariables :
# ~ sexe + condition  (pour corriger l'effet du sexe)
# ----------------------------------------------------------------
counts  <- read.csv(count_mat, row.names = 1, check.names = FALSE)
coldata <- read.csv(coldata_f, row.names = 1)
coldata$condition <- factor(coldata$condition,
                            levels = c("control", "treated"))

stopifnot(all(colnames(counts) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~ condition
)

# ----------------------------------------------------------------
# 3. Filtrage préliminaire
# ----------------------------------------------------------------
# Règle pratique : garder les gènes avec au moins 10 reads
# dans au moins 3 échantillons.
# Justification : un gène avec counts = 0,0,0,1,0,0 n'apporte
# aucune information statistique exploitable.
# ----------------------------------------------------------------
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message(sprintf("Gènes retenus après filtrage : %d", nrow(dds)))

# ----------------------------------------------------------------
# 4. Estimation des size factors (normalisation)
# ----------------------------------------------------------------
# Problème : certains échantillons ont plus de reads que d'autres.
# Solution DESeq2 : médiane des ratios (pas TPM, pas RPKM).
#
# Pour chaque gène i et échantillon j :
# size_factor_j = médiane_i( K_ij / (prod_j K_ij)^(1/n) )
#
# Cette méthode est robuste aux gènes DE (contrairement à la
# normalisation par somme totale).
# ----------------------------------------------------------------
dds <- DESeq(dds, parallel = FALSE)

message("Size factors calculés :")
print(round(sizeFactors(dds), 3))

# ----------------------------------------------------------------
# 5. Shrinkage apeglm — correction des LFC extrêmes
# ----------------------------------------------------------------
# Problème : les gènes peu exprimés ont des LFC artificiellement
# élevés (peu de reads → grande variance → grands LFC par hasard).
#
# Solution : shrinkage bayésien apeglm
# Les LFC des gènes peu exprimés sont "tirés" vers 0.
# Les gènes fortement exprimés conservent leur LFC réel.
#
# Référence : Zhu et al., Bioinformatics, 2019
# ----------------------------------------------------------------
res_raw <- results(dds,
                   contrast = c("condition", "treated", "control"),
                   alpha    = padj_thr)

res <- lfcShrink(dds,
                 coef = "condition_treated_vs_control",
                 type = "apeglm",
                 quiet = TRUE)

message(sprintf(
  "Gènes significatifs (padj<%.2f, |LFC|>%.1f) : %d",
  padj_thr, lfc_thr,
  sum(res$padj < padj_thr & abs(res$log2FoldChange) > lfc_thr,
      na.rm = TRUE)
))

# ----------------------------------------------------------------
# 6. Export et visualisations
# ----------------------------------------------------------------
res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::mutate(
    significance = dplyr::case_when(
      padj < padj_thr & log2FoldChange >  lfc_thr ~ "Up",
      padj < padj_thr & log2FoldChange < -lfc_thr ~ "Down",
      TRUE ~ "NS"
    )
  ) |>
  dplyr::arrange(padj)

write.csv(res_df, out_de, row.names = FALSE)
message("Résultats exportés : ", out_de)