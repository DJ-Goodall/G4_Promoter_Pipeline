#!/usr/bin/env Rscript
# ============================================================================
# 24_feature_total_metaprofiles.R   (env: r_g4)   --- Phase 10 ---
#
# SIMPLIFIED companion to rule 14: meta-profiles of TOTAL G4 CUT&Tag signal
# centred on each genomic feature, WITHOUT splitting by G4 topology. One line
# per genotype (WT + 3 KOs), one panel per feature. Gives the overall G4-
# occupancy landscape at each feature class and any KO shift at a glance.
#
# Features (config metaprofile.total_feature_regions):
#   promoter, 5'UTR (per-fragment), 1st intron, enhancer  -> regions_<x>.rds in
#   cache_v2 (via feature.region_files); plus lncRNA -> TSS windows built from
#   cache/lncrna_genes.rds (rule 21).
#
# "Total" = ALL feature instances, with NO G4-peak / topology filter (only the
# seeded max_per_group subsample for runtime). Profiles are centred on the
# feature centre (lncRNA: the gene TSS) and strand-aware (mean_replicate_profile
# flips - strand; enhancers are unstranded).
#
# Two figures: per-feature free y-axis (shape) and one shared, AUTO-DERIVED
# y-axis across all features (absolute-intensity comparison; no clipping).
#
# Inputs:   cache_v2 regions_<feature>.rds, cache/lncrna_genes.rds, G4 bigWigs
# Outputs:  results/tables/feature_total_profiles.csv
#           results/figures/10_metaprofiles/feature_signal_total.{pdf,png}
#           results/figures/10_metaprofiles/feature_signal_total_fixedY.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles", "cache")
fig_root <- "results/figures"

genotypes  <- cfg$genotypes
half_width <- cfg$metaprofile$half_width
n_bins     <- cfg$metaprofile$n_bins
max_group  <- cfg$metaprofile$max_per_group %||% 8000

features <- as.character(unlist(cfg$metaprofile$total_feature_regions %||%
                                c("promoter", "5UTR", "intron1", "enhancer", "lncRNA")))
region_files <- cfg$feature$region_files

# --- Helper: standardise seqlevels + chrom lengths on a feature GRanges ------
std_seqlevels <- function(gr) {
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepSeqlevels(
    gr, intersect(GenomeInfoDb::seqlevels(gr), cfg$std_chroms), pruning.mode = "coarse")
  GenomeInfoDb::seqlengths(gr) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(gr)]
  gr
}

# --- Build one centre-anchored window set per feature (no topology) ----------
set.seed(7)
region_windows <- lapply(features, function(rn) {
  if (rn == "lncRNA") {
    # lncRNA "feature" = lncRNA gene TSS (strand-aware promoters(); midpoint=TSS).
    lnc <- readRDS("cache/lncrna_genes.rds")
    lnc <- std_seqlevels(lnc)
    win <- GenomicRanges::trim(
      GenomicRanges::promoters(lnc, upstream = half_width, downstream = half_width))
  } else {
    rf <- region_files[[rn]]
    if (is.null(rf)) { message("No region file mapped for '", rn, "'; skipping."); return(NULL) }
    gr <- readRDS(file.path(cfg$paths$cache_v2, rf))
    gr <- std_seqlevels(gr)
    if (length(gr) == 0) return(NULL)
    win <- GenomicRanges::trim(center_window(gr, width = 2 * half_width))
  }
  if (length(win) == 0) return(NULL)
  # Cap the feature for tractable runtime (random, seeded) — same rule as 14.
  if (length(win) > max_group) win <- win[sort(sample(seq_along(win), max_group))]
  message(rn, " total windows: ", length(win))
  win
})
names(region_windows) <- features
region_windows <- region_windows[!vapply(region_windows, is.null, logical(1))]
stopifnot(length(region_windows) > 0)

# --- Compute per-genotype, per-feature TOTAL meta-profiles -------------------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)

# Cache key: window params, cap, and feature set (delete cache if any changes).
prof_cache <- sprintf("cache/feature_total_profiles_hw%d_b%d_cap%d_%s.rds",
                      half_width, n_bins, max_group,
                      paste(names(region_windows), collapse = "-"))
profiles <- cache_or_build(prof_cache, {
  dplyr::bind_rows(lapply(genotypes, function(g) {
    bws <- bw_meta$filepath[bw_meta$genotype == g]
    if (length(bws) == 0) return(NULL)
    message("Total feature profiles: ", g, " (", length(bws), " bigWigs)")
    rows <- dplyr::bind_rows(lapply(names(region_windows), function(rn) {
      regions <- region_windows[[rn]]
      if (length(regions) < 5) return(NULL)
      prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
      prof$genotype  <- g
      prof$feature   <- rn
      prof$n_regions <- length(regions)
      prof
    }))
    evict_bigwigs(bws)
    rows
  }))
})
present <- intersect(features, names(region_windows))     # keep config order
profiles$genotype <- factor(profiles$genotype, levels = genotypes)
profiles$feature  <- factor(profiles$feature, levels = present)
readr::write_csv(profiles, "results/tables/feature_total_profiles.csv")

# --- Plot: facet_wrap(~feature), genotype lines -----------------------------
region_labels_map <- c("promoter" = "Promoter", "5UTR" = "5'UTR",
                       "intron1" = "1st intron", "enhancer" = "Enhancer (dELS)",
                       "genebody" = "Gene body", "lncRNA" = "lncRNA (TSS)")
profiles$feature_label <- factor(
  unname(region_labels_map[as.character(profiles$feature)]),
  levels = unname(region_labels_map[present]))

base_layers <- list(
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA),
  geom_line(linewidth = 0.8),
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40"),
  scale_colour_manual(values = condition_colours),
  scale_fill_manual(values = condition_colours),
  labs(x = "Distance to feature centre (bp)", y = "Mean G4 CUT&Tag signal",
       colour = "Genotype", fill = "Genotype"),
  theme_pub(), theme(legend.position = "bottom"))
fig_w <- 2.7 * length(present) + 1

# (1) per-feature free y-axis (shape).
p_free <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_wrap(~ feature_label, nrow = 1, scales = "free_y") +
  labs(title = "Total G4 signal at genomic features (per-feature y-axis)",
       subtitle = sprintf("+/- %d bp, %d bins; all feature instances, no topology split; lines = genotype",
                          half_width, n_bins))
save_plot(p_free, "feature_signal_total", "10_metaprofiles", fig_root,
          width = fig_w, height = 4.2)

# (2) shared, auto-derived y-axis across all features (intensity comparison).
ymax <- max(profiles$mean + profiles$sem, na.rm = TRUE) * 1.05
p_fixed <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_wrap(~ feature_label, nrow = 1) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(title = "Total G4 signal at genomic features (shared y-axis)",
       subtitle = sprintf("+/- %d bp, %d bins; shared auto y-axis 0-%.2g (max signal +SEM x1.05) for cross-feature comparison",
                          half_width, n_bins, ymax))
save_plot(p_fixed, "feature_signal_total_fixedY", "10_metaprofiles", fig_root,
          width = fig_w, height = 4.2)

message("Done. Wrote total feature meta-profiles (free + shared auto-y) for: ",
        paste(present, collapse = ", "))
