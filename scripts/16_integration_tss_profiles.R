#!/usr/bin/env Rscript
# ============================================================================
# 16_integration_tss_profiles.R   (system R)   --- Stage E ---
#
# CUT&Tag metaprofiles around the TSS (+/-5 kb, 100 bins) of top DEGs (up vs
# down), WT vs each KO, both assays. Ported from analysis_V2.Rmd chunk
# integration-tss-profile. Uses the vectorized, strand-aware mean_replicate_
# profile helper (the RMD's row-by-row compute_tss_profile was NOT strand-aware;
# this is the corrected behaviour).
#
# Inputs:   cache/deseq2_res_list.rds, cache/gene_meta.rds, config bigwig_dir
# Outputs:  results/figures/07_integration/{tss_profile_combined_<KO>_vs_WT,
#             tss_profile_{BG4,S96}_<KO>_vs_WT}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(clusterProfiler)
  library(org.Mm.eg.db); library(ggplot2); library(cowplot); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures")
fig_root <- "results/figures"
ig <- cfg$integration
hw <- ig$tss_half_width; nb <- ig$tss_n_bins

res_list  <- readRDS("cache/deseq2_res_list.rds")
gene_meta <- readRDS("cache/gene_meta.rds")
bw_meta   <- build_bw_meta(cfg$paths$bigwig_dir,
                           genotypes = c("WT", "DHX36KO", "FANCJKO", "dKO"), assay = NULL)

gene_map <- clusterProfiler::bitr(unique(gene_meta$gene_name), fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
txdb     <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- GenomicFeatures::genes(txdb)
tss_5kb  <- GenomicFeatures::promoters(genes_gr, upstream = hw, downstream = hw)

top_deg_tss <- function(nm, dir) {
  df <- res_list[[nm]]
  df$entrez <- gene_map$ENTREZID[match(df$gene_name, gene_map$SYMBOL)]
  df <- df[!is.na(df$entrez) & !is.na(df$padj) & df$padj < cfg$rnaseq$padj, ]
  df <- if (dir == "up") df[df$log2FoldChange > 0, ][order(-df$log2FoldChange[df$log2FoldChange > 0]), ]
        else            df[df$log2FoldChange < 0, ][order(df$log2FoldChange[df$log2FoldChange < 0]), ]
  ent <- head(df$entrez, ig$top_degs)
  tss_5kb[names(tss_5kb) %in% ent]
}

kos <- c("DHX36KO", "FANCJKO", "dKO")
for (ko in kos) {
  nm <- paste0(ko, "_vs_WT")
  regs <- list(up = top_deg_tss(nm, "up"), down = top_deg_tss(nm, "down"))
  parts <- list()
  for (cond in c("WT", ko)) for (a in c("G4_BG4", "Rloop_S96")) {
    fp <- bw_meta$filepath[bw_meta$genotype == cond & bw_meta$assay == a]
    if (length(fp) == 0) next
    for (d in c("up", "down")) {
      if (length(regs[[d]]) == 0) next
      prof <- mean_replicate_profile(regs[[d]], fp, n_bins = nb, half_width = hw)
      prof$condition <- cond; prof$assay <- a; prof$direction <- d
      parts[[paste(cond, a, d, sep = "_")]] <- prof
    }
    evict_bigwigs(fp)
  }
  prof_df <- dplyr::bind_rows(parts)
  prof_df$assay_label <- unname(assay_label_map[prof_df$assay])
  prof_df$direction <- factor(prof_df$direction, levels = c("up", "down"),
                              labels = c("Upregulated", "Downregulated"))
  prof_df$condition <- factor(prof_df$condition, levels = c("WT", ko))

  pal <- c(Upregulated = "coral", Downregulated = "steelblue")
  p_comb <- ggplot(prof_df, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_grid(assay_label ~ condition, scales = "free_y") +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    labs(title = paste0("CUT&Tag signal around DEG TSS (", ko, " vs WT)"),
         subtitle = "TSS +/- 5kb, 100 bins", x = "Position relative to TSS (bp)",
         y = "Mean signal", colour = "Gene set", fill = "Gene set") +
    theme_cowplot() + theme(strip.background = element_rect(fill = "grey90"))
  save_plot(p_comb, paste0("tss_profile_combined_", ko, "_vs_WT"), "07_integration",
            fig_root, width = 12, height = 8)

  for (a in c("G4_BG4", "Rloop_S96")) {
    short <- if (a == "G4_BG4") "BG4" else "S96"
    p <- ggplot(dplyr::filter(prof_df, assay == a),
                aes(position, mean, colour = direction, fill = direction)) +
      geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.2, colour = NA) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      facet_wrap(~ condition, nrow = 1) +
      scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
      labs(title = paste(assay_label_map[a], "signal around DEG TSS (", ko, "vs WT)"),
           x = "Position relative to TSS (bp)", y = "Mean signal") + theme_cowplot()
    save_plot(p, paste0("tss_profile_", short, "_", ko, "_vs_WT"), "07_integration",
              fig_root, width = 12, height = 6)
  }
}
message("Integration TSS profiles done.")
