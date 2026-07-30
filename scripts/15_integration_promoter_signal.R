#!/usr/bin/env Rscript
# ============================================================================
# 15_integration_promoter_signal.R   (system R)   --- Stage E ---
#
# CUT&Tag signal at promoters of the top DEGs (up vs down), WT vs each KO, both
# assays, with a Wilcoxon up-vs-down test. Ported from analysis_V2.Rmd chunks
# integration-tss / integration-top-genes / integration-promoter-signal. Signal
# via the Windows-safe mean_replicate_signal helper.
#
# Inputs:   cache/deseq2_res_list.rds, cache/gene_meta.rds, config bigwig_dir
# Outputs:  results/figures/07_integration/promoter_signal_boxplot_<KO>_vs_WT.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(clusterProfiler)
  library(org.Mm.eg.db); library(ggplot2); library(cowplot); library(ggpubr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures")
fig_root <- "results/figures"
ig <- cfg$integration

res_list  <- readRDS("cache/deseq2_res_list.rds")
gene_meta <- readRDS("cache/gene_meta.rds")
bw_meta   <- build_bw_meta(cfg$paths$bigwig_dir,
                           genotypes = c("WT", "DHX36KO", "FANCJKO", "dKO"), assay = NULL)

gene_map <- clusterProfiler::bitr(unique(gene_meta$gene_name), fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr     <- GenomicFeatures::genes(txdb)
promoter_2kb <- GenomicFeatures::promoters(genes_gr, upstream = ig$promoter_half_width,
                                           downstream = ig$promoter_half_width)

top_deg_promoters <- function(nm, dir) {
  df <- res_list[[nm]]
  df$entrez <- gene_map$ENTREZID[match(df$gene_name, gene_map$SYMBOL)]
  df <- df[!is.na(df$entrez) & !is.na(df$padj) & df$padj < cfg$rnaseq$padj, ]
  df <- if (dir == "up") df[df$log2FoldChange > 0, ][order(-df$log2FoldChange[df$log2FoldChange > 0]), ]
        else            df[df$log2FoldChange < 0, ][order(df$log2FoldChange[df$log2FoldChange < 0]), ]
  ent <- head(df$entrez, ig$top_degs)
  promoter_2kb[names(promoter_2kb) %in% ent]
}

kos <- c("DHX36KO", "FANCJKO", "dKO")
for (ko in kos) {
  nm <- paste0(ko, "_vs_WT")
  proms <- list(up = top_deg_promoters(nm, "up"), down = top_deg_promoters(nm, "down"))
  parts <- list()
  for (cond in c("WT", ko)) for (a in c("G4_BG4", "Rloop_S96")) {
    fp <- bw_meta$filepath[bw_meta$genotype == cond & bw_meta$assay == a]
    if (length(fp) == 0) next
    for (d in c("up", "down")) {
      if (length(proms[[d]]) == 0) next
      sig <- mean_replicate_signal(proms[[d]], fp)
      parts[[paste(cond, a, d, sep = "_")]] <- data.frame(
        signal = sig, condition = cond, assay = a, direction = d)
    }
    evict_bigwigs(fp)
  }
  pdf_df <- dplyr::bind_rows(parts)
  pdf_df <- pdf_df[!is.na(pdf_df$signal), ]
  pdf_df$assay_label <- unname(assay_label_map[pdf_df$assay])
  pdf_df$direction <- factor(pdf_df$direction, levels = c("up", "down"),
                             labels = c("Upregulated", "Downregulated"))
  pdf_df$condition <- factor(pdf_df$condition, levels = c("WT", ko))

  p <- ggplot(pdf_df, aes(direction, log2(signal + 1), fill = direction)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif",
      comparisons = list(c("Upregulated", "Downregulated"))) +
    facet_grid(condition ~ assay_label, scales = "free_y") +
    scale_fill_manual(values = c(Upregulated = "coral", Downregulated = "steelblue")) +
    labs(title = paste0("CUT&Tag signal at DEG promoters (", ko, " vs WT)"),
         subtitle = "Promoter = TSS +/- 2kb", y = "log2(mean signal + 1)", x = "") +
    theme_cowplot() + theme(legend.position = "none",
      strip.background = element_rect(fill = "grey90"))
  save_plot(p, paste0("promoter_signal_boxplot_", ko, "_vs_WT"), "07_integration",
            fig_root, width = 12, height = 8)
}
message("Integration promoter signal done.")
