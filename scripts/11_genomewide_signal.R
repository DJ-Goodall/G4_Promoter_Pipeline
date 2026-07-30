#!/usr/bin/env Rscript
# ============================================================================
# 11_genomewide_signal.R   (system R)   --- Stage C ---
#
# Genome-wide CUT&Tag signal QC: 1 kb tiles (subsampled), per-bigWig binned
# signal, then violin, per-assay correlation heatmap, mean-signal bar, and the
# G4-vs-R-loop per-genotype scatter panel. Ported from analysis_V2.Rmd §5
# (cuttag-* chunks). Uses the Windows-safe read_bigwig_rle and evicts each
# bigWig after extraction to bound memory across all 28 files.
#
# Inputs:   config paths.bigwig_dir (28 bigWigs)
# Outputs:  cache/genomewide_signal_mat.rds
#           results/figures/05_cuttag/{violin_signal_<a>,correlation_heatmap_<a>,
#             mean_signal_bar_<a>,G4_vs_Rloop_correlation_panel}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
  library(BSgenome.Mmusculus.UCSC.mm10); library(ggplot2); library(cowplot)
  library(pheatmap); library(RColorBrewer); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures")
fig_root <- "results/figures"
gw <- cfg$genomewide

std_chroms  <- cfg$std_chroms
chrom_sizes <- default_chrom_sizes(std_chroms)
all_geno    <- c("WT", "DHX36KO", "FANCJKO", "dKO", "ERCCWT", "ERCCKO")
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = all_geno, assay = NULL)

sig_cache <- "cache/genomewide_signal_mat.rds"
sig_obj <- cache_or_build(sig_cache, {
  tiles <- GenomicRanges::tileGenome(chrom_sizes, tilewidth = gw$window,
                                     cut.last.tile.in.chrom = TRUE)
  set.seed(gw$subsample_seed)
  n_s <- min(gw$subsample, length(tiles))
  tiles_s <- tiles[sample(length(tiles), n_s)]
  message("Sampled ", length(tiles_s), " of ", length(tiles), " 1kb windows")

  smat <- matrix(NA_real_, nrow = length(tiles_s), ncol = nrow(bw_meta))
  colnames(smat) <- bw_meta$filename
  tile_chr <- as.character(GenomicRanges::seqnames(tiles_s))
  for (i in seq_len(nrow(bw_meta))) {
    message("  [", i, "/", nrow(bw_meta), "] ", bw_meta$filename[i])
    cov <- read_bigwig_rle(bw_meta$filepath[i], chroms = std_chroms, chrom_sizes = chrom_sizes)
    shared <- intersect(names(cov), std_chroms)
    tsub <- tiles_s[tile_chr %in% shared]
    GenomeInfoDb::seqlevels(tsub, pruning.mode = "coarse") <- shared
    GenomeInfoDb::seqlengths(tsub) <- chrom_sizes[shared]
    ba <- GenomicRanges::binnedAverage(tsub, cov[shared], "score")
    smat[tile_chr %in% shared, i] <- ba$score
    evict_bigwigs(bw_meta$filepath[i])
  }
  smat[is.na(smat)] <- 0
  list(signal_mat = smat, bw_meta = bw_meta)
})
signal_mat <- sig_obj$signal_mat
bw_meta    <- sig_obj$bw_meta

assays <- c("G4_BG4", "Rloop_S96")

# --- 5b. Violin (per assay, nonzero subsample) ------------------------------
signal_long <- as.data.frame(signal_mat) %>%
  mutate(window_id = seq_len(nrow(signal_mat))) %>%
  tidyr::pivot_longer(-window_id, names_to = "filename", values_to = "signal") %>%
  dplyr::left_join(bw_meta[, c("filename", "assay", "genotype")], by = "filename")

for (a in assays) {
  pd <- signal_long %>% dplyr::filter(assay == a, signal > 0) %>%
    dplyr::group_by(filename) %>%
    dplyr::slice_sample(n = gw$violin_subsample) %>% dplyr::ungroup()
  p <- ggplot(pd, aes(genotype, log2(signal + 1), fill = genotype)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = condition_colours) +
    labs(title = paste(assay_label_map[a], "— genome-wide signal"),
         y = "log2(signal + 1)", x = "") + theme_cowplot() +
    theme(legend.position = "none")
  save_plot(p, paste0("violin_signal_", a), "05_cuttag", fig_root, width = 12, height = 6)
}

# --- 5c. Per-assay correlation heatmaps -------------------------------------
for (a in assays) {
  files <- bw_meta$filename[bw_meta$assay == a]
  if (length(files) < 2) next
  cor_bw <- cor(signal_mat[, files], method = "pearson", use = "complete.obs")
  ann <- data.frame(Genotype = bw_meta$genotype[match(files, bw_meta$filename)],
                    row.names = files)
  hm_fn <- function() pheatmap::pheatmap(cor_bw, annotation_col = ann,
    annotation_colors = list(Genotype = condition_colours),
    color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
    display_numbers = TRUE, number_format = "%.3f", fontsize_number = 7,
    main = paste(assay_label_map[a], "— sample Pearson correlation"))
  save_base_plot(hm_fn, paste0("correlation_heatmap_", a), "05_cuttag", fig_root,
                 width = 10, height = 8)
}

# --- 5d. Mean signal bar ----------------------------------------------------
sample_means <- data.frame(filename = colnames(signal_mat),
                           mean_signal = colMeans(signal_mat, na.rm = TRUE)) %>%
  dplyr::left_join(bw_meta[, c("filename", "assay", "genotype")], by = "filename")
mean_summary <- sample_means %>% dplyr::group_by(assay, genotype) %>%
  dplyr::summarise(mean_sig = mean(mean_signal), sd_sig = sd(mean_signal),
                   .groups = "drop")
for (a in assays) {
  pd <- mean_summary %>% dplyr::filter(assay == a)
  p <- ggplot(pd, aes(genotype, mean_sig, fill = genotype)) +
    geom_col(colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_sig - sd_sig, ymax = mean_sig + sd_sig), width = 0.2) +
    geom_jitter(data = dplyr::filter(sample_means, assay == a),
                aes(genotype, mean_signal), width = 0.1, size = 2, shape = 21, fill = "white") +
    scale_fill_manual(values = condition_colours) +
    labs(title = paste(assay_label_map[a], "— mean genome-wide signal"),
         y = "Mean signal per 1kb window", x = "") + theme_cowplot() +
    theme(legend.position = "none")
  save_plot(p, paste0("mean_signal_bar_", a), "05_cuttag", fig_root, width = 10, height = 6)
}

# --- 5e. G4 vs R-loop scatter panel -----------------------------------------
geno_mean <- function(a) {
  files <- bw_meta[bw_meta$assay == a, ]
  ms <- lapply(levels(bw_meta$genotype), function(g) {
    gf <- files$filename[files$genotype == g]
    if (length(gf) > 0) rowMeans(signal_mat[, gf, drop = FALSE], na.rm = TRUE) else NULL
  })
  names(ms) <- levels(bw_meta$genotype)
  do.call(cbind, ms[!vapply(ms, is.null, logical(1))])
}
g4m <- geno_mean("G4_BG4"); rlm <- geno_mean("Rloop_S96")
shared_g <- intersect(colnames(g4m), colnames(rlm))
scatter_plots <- lapply(shared_g, function(g) {
  df <- data.frame(G4 = g4m[, g], Rloop = rlm[, g]) %>% dplyr::filter(G4 > 0 | Rloop > 0)
  if (nrow(df) > gw$scatter_subsample) df <- df[sample(nrow(df), gw$scatter_subsample), ]
  r <- cor(log2(df$G4 + 1), log2(df$Rloop + 1), method = "pearson")
  ggplot(df, aes(log2(G4 + 1), log2(Rloop + 1))) +
    geom_point(alpha = 0.03, size = 0.2, colour = "grey30") +
    geom_smooth(method = "loess", colour = "red", se = TRUE, linewidth = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0("r = ", round(r, 3)),
             hjust = 1.2, vjust = 2, size = 5, colour = "red") +
    labs(title = g, x = "BG4 (G4) — log2(signal+1)", y = "S9.6 (R-loop) — log2(signal+1)") +
    theme_cowplot()
})
p_panel <- cowplot::plot_grid(plotlist = scatter_plots, ncol = 2)
save_plot(p_panel, "G4_vs_Rloop_correlation_panel", "05_cuttag", fig_root, width = 12, height = 15)
message("Genome-wide CUT&Tag QC done.")
