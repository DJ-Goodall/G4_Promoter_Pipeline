#!/usr/bin/env Rscript
# ============================================================================
# 03_deseq2_plots.R   (system R)   --- Stage A ---
#
# DESeq2 result figures: per-contrast volcano (EnhancedVolcano) + MA plots, union
# -DEG heatmap, and DEG-overlap UpSet. Ported from analysis_V2.Rmd chunks
# de-volcano / de-ma / de-heatmap / de-upset. Gene labels use gene_name; the
# stable key throughout is gene_id.
#
# Inputs:   cache/deseq2_res_list.rds, cache/rld.rds, cache/sample_info.rds
# Outputs:  results/tables/significant_DEG_union.tsv
#           results/figures/02_deseq2/{volcano_<c>,ma_plot_<c>,union_DEG_heatmap,
#             upset_DEGs}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(EnhancedVolcano); library(ggplot2); library(cowplot)
  library(pheatmap); library(RColorBrewer); library(UpSetR); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures", "results/tables")
fig_root <- "results/figures"

res_list    <- readRDS("cache/deseq2_res_list.rds")
rld         <- readRDS("cache/rld.rds")
sample_info <- readRDS("cache/sample_info.rds")
padj_cut <- cfg$rnaseq$padj; lfc_cut <- cfg$rnaseq$lfc
top_n    <- cfg$rnaseq$volcano_top_n %||% 15

# --- 3a/3b. Volcano + MA per contrast ---------------------------------------
for (nm in names(res_list)) {
  df <- res_list[[nm]]
  lab <- df$gene_name; lab[is.na(lab)] <- df$gene_id[is.na(lab)]
  sel <- lab[order(df$padj)][seq_len(min(top_n, nrow(df)))]
  p_v <- EnhancedVolcano::EnhancedVolcano(df, lab = lab,
    x = "log2FoldChange", y = "padj", title = nm,
    subtitle = sprintf("padj < %g, |log2FC| > %g", padj_cut, lfc_cut),
    pCutoff = padj_cut, FCcutoff = lfc_cut, pointSize = 1.5, labSize = 3,
    selectLab = sel, col = c("grey70", "grey70", "grey70", "coral"),
    drawConnectors = TRUE, widthConnectors = 0.5, maxoverlapsConnectors = 20)
  save_plot(p_v, paste0("volcano_", nm), "02_deseq2", fig_root, width = 10, height = 8)

  ma <- df %>% mutate(sig = case_when(
    is.na(padj) ~ "NS",
    padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
    padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down", TRUE ~ "NS"))
  p_ma <- ggplot(ma, aes(log10(baseMean + 1), log2FoldChange, colour = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_manual(values = c(Up = "coral", Down = "steelblue", NS = "grey70")) +
    geom_hline(yintercept = c(-lfc_cut, 0, lfc_cut),
               linetype = c("dashed", "solid", "dashed"),
               colour = c("steelblue", "black", "coral"), linewidth = 0.4) +
    labs(title = paste("MA plot —", nm), x = "log10(mean expression + 1)",
         y = "log2 fold change") + theme_cowplot()
  save_plot(p_ma, paste0("ma_plot_", nm), "02_deseq2", fig_root, width = 10, height = 6)
}

# --- 3d. Union-DEG heatmap (3 vs-WT contrasts) ------------------------------
vs_wt <- names(res_list)[grepl("_vs_WT$", names(res_list))][1:3]
sig_per <- lapply(res_list[vs_wt], function(df)
  df$gene_id[!is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) > lfc_cut])
sig_union <- unique(unlist(sig_per))
message("Union DEGs across 3 vs-WT: ", length(sig_union))

gm <- setNames(res_list[[1]]$gene_name, res_list[[1]]$gene_id)
union_df <- data.frame(gene_id = sig_union, gene_name = gm[sig_union])
for (nm in vs_wt) union_df[[nm]] <- union_df$gene_id %in% sig_per[[nm]]
readr::write_tsv(union_df, "results/tables/significant_DEG_union.tsv")

ann_col <- data.frame(Condition = sample_info$condition, row.names = rownames(sample_info))
ann_colours <- list(Condition = condition_colours)
sig_in_rld <- sig_union[sig_union %in% rownames(SummarizedExperiment::assay(rld))]
if (length(sig_in_rld) > 1) {
  z <- t(scale(t(SummarizedExperiment::assay(rld)[sig_in_rld, ])))
  z[z > 3] <- 3; z[z < -3] <- -3
  hm_fn <- function() pheatmap::pheatmap(z, annotation_col = ann_col,
    annotation_colors = ann_colours, show_rownames = FALSE,
    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    main = paste0("Union DEGs (n=", length(sig_in_rld), ") — z-scored rlog, capped +/-3"))
  save_base_plot(hm_fn, "union_DEG_heatmap", "02_deseq2", fig_root, width = 10, height = 12)
}

# --- 3e. UpSet of DEG overlap -----------------------------------------------
deg_sets <- setNames(lapply(sig_per, as.character), c("DHX36KO", "FANCJKO", "dKO"))
upset_fn <- function() print(UpSetR::upset(UpSetR::fromList(deg_sets),
  order.by = "freq", nsets = 3,
  sets.bar.color = unname(condition_colours[c("DHX36KO", "FANCJKO", "dKO")]),
  main.bar.color = "grey30", text.scale = 1.3))
save_base_plot(upset_fn, "upset_DEGs", "02_deseq2", fig_root, width = 10, height = 6)

message("DESeq2 plots done.")
