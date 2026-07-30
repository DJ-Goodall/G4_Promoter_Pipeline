#!/usr/bin/env Rscript
# ============================================================================
# 14_region_metaprofiles.R   (system R)   --- Stage D ---
#
# Meta-profiles around region centres (+/-2 kb), both a fast subsampled version
# (profile_n) and a larger validation version (profile_full_cap). Ported from
# extended_V3.Rmd chunks s3-region-profile(-full). Uses center_window (== the
# RMD's make_center_window) + the vectorized strand-aware mean_replicate_profile.
#
# Inputs:   cache/regions_all.rds, config bigwig_dir
# Outputs:  cache/{region_profile_df,region_profile_full_df}.rds
#           results/figures/06_regional/{02_region_metaprofile,
#             02_region_metaprofile_full_peaks}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(dplyr)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures")
fig_root <- "results/figures"
rgc  <- cfg$regional
nb   <- rgc$n_bins; half <- rgc$half_width

bw_meta     <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = extended_genotypes, assay = NULL)
regions_all <- readRDS("cache/regions_all.rds")

build_profiles <- function(cap, seed) {
  set.seed(seed)
  rg_wins <- lapply(regions_all, function(rg) {
    idx <- if (length(rg) > cap) sample(length(rg), cap) else seq_along(rg)
    center_window(rg[idx], width = 2 * half) })
  dplyr::bind_rows(lapply(unique(bw_meta$assay), function(a) {
    dplyr::bind_rows(lapply(extended_genotypes, function(g) {
      bws <- bw_meta$filepath[bw_meta$assay == a & bw_meta$genotype == g]
      if (!length(bws)) return(NULL)
      rows <- dplyr::bind_rows(lapply(names(rg_wins), function(rn) {
        prof <- mean_replicate_profile(rg_wins[[rn]], bws, n_bins = nb, half_width = half)
        prof$region <- rn; prof$assay <- a; prof$genotype <- g; prof }))
      evict_bigwigs(bws); rows
    }))
  }))
}
fmt <- function(df) {
  df$region_label <- factor(df$region, levels = region_levels, labels = region_labels)
  df$assay_label  <- unname(assay_label_map[df$assay])
  df$genotype     <- factor(df$genotype, levels = extended_genotypes); df
}
plot_profiles <- function(df, name, subtitle) {
  p <- ggplot(df, aes(position, mean, colour = genotype, fill = genotype)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = condition_colours) +
    scale_fill_manual(values = condition_colours) +
    facet_grid(assay_label ~ region_label, scales = "free_y") +
    labs(x = "Distance to feature centre (bp)", y = "Mean signal",
         title = "Meta-profile around feature centres (+/-2 kb)", subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, name, "06_regional", fig_root, width = 12, height = 7)
}

profile_df <- fmt(cache_or_build("cache/region_profile_df.rds",
                                 build_profiles(rgc$profile_n, rgc$profile_seed)))
plot_profiles(profile_df, "02_region_metaprofile",
  sprintf("Subsampled to %d regions per class, %d bins", rgc$profile_n, nb))

profile_full_df <- fmt(cache_or_build("cache/region_profile_full_df.rds",
                                      build_profiles(rgc$profile_full_cap, rgc$profile_full_seed)))
plot_profiles(profile_full_df, "02_region_metaprofile_full_peaks",
  sprintf("Randomly sampled %d peaks/region (validation vs N=%d), %d bins",
          rgc$profile_full_cap, rgc$profile_n, nb))
message("Region metaprofiles done.")
