#!/usr/bin/env Rscript
# ============================================================================
# 13_region_signal.R   (system R)   --- Stage D ---
#
# Mean CUT&Tag replicate signal per genomic feature (subsampled), per genotype x
# assay, with a boxplot and per-region Wilcoxon tests vs matched controls.
# Ported from extended_V3.Rmd chunks s3-region-signal / s3-region-signal-box.
#
# Inputs:   cache/regions_all.rds, config bigwig_dir
# Outputs:  cache/region_signal_df.rds
#           results/tables/region_signal_stats.tsv
#           results/figures/06_regional/02_region_signal_boxplot.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(dplyr); library(readr)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures", "results/tables")
fig_root <- "results/figures"
rgc <- cfg$regional
N   <- rgc$region_signal_n

bw_meta     <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = extended_genotypes, assay = NULL)
regions_all <- readRDS("cache/regions_all.rds")

region_signal_df <- cache_or_build("cache/region_signal_df.rds", {
  set.seed(rgc$region_signal_seed)
  rg_subs <- lapply(regions_all, function(rg) {
    idx <- if (length(rg) > N) sample(length(rg), N) else seq_along(rg); rg[idx] })
  dplyr::bind_rows(lapply(unique(bw_meta$assay), function(a) {
    dplyr::bind_rows(lapply(extended_genotypes, function(g) {
      bws <- bw_meta$filepath[bw_meta$assay == a & bw_meta$genotype == g]
      if (!length(bws)) return(NULL)
      rows <- dplyr::bind_rows(lapply(names(rg_subs), function(rn)
        data.frame(region = rn, assay = a, genotype = g,
                   signal = mean_replicate_signal(rg_subs[[rn]], bws))))
      evict_bigwigs(bws); rows
    }))
  }))
})
region_signal_df$region      <- factor(region_signal_df$region, levels = region_levels, labels = region_labels)
region_signal_df$genotype    <- factor(region_signal_df$genotype, levels = extended_genotypes)
region_signal_df$assay_label <- unname(assay_label_map[region_signal_df$assay])
region_signal_df <- region_signal_df[!is.na(region_signal_df$signal), ]

p_sig <- ggplot(region_signal_df, aes(genotype, log2(signal + 1), fill = genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, colour = "grey20") +
  scale_fill_manual(values = condition_colours) +
  facet_grid(assay_label ~ region, scales = "free_y") +
  coord_cartesian(ylim = c(0, quantile(log2(region_signal_df$signal + 1), 0.995, na.rm = TRUE))) +
  labs(x = NULL, y = "log2(mean bigwig signal + 1)", title = "CUT&Tag signal per genomic feature",
       subtitle = sprintf("Subsampled to %d regions per class", N)) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_plot(p_sig, "02_region_signal_boxplot", "06_regional", fig_root, width = 14, height = 7)

sig_stats <- dplyr::bind_rows(lapply(list(c("DHX36KO", "WT"), c("FANCJKO", "WT"),
  c("dKO", "WT"), c("ERCCKO", "ERCCWT")), function(pair) {
  ko <- pair[1]; ctrl <- pair[2]
  region_signal_df %>% dplyr::filter(genotype %in% pair) %>%
    dplyr::group_by(assay, assay_label, region) %>%
    dplyr::summarise(comparison = paste0(ko, "_vs_", ctrl),
      wilcox_p = tryCatch(wilcox.test(signal ~ genotype)$p.value, error = function(e) NA_real_),
      median_ctrl = median(signal[genotype == ctrl], na.rm = TRUE),
      median_ko = median(signal[genotype == ko], na.rm = TRUE),
      n_ctrl = sum(genotype == ctrl), n_ko = sum(genotype == ko), .groups = "drop")
}))
readr::write_tsv(sig_stats, "results/tables/region_signal_stats.tsv")
message("Region signal done.")
