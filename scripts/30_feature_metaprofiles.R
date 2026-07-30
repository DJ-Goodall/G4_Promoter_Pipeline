#!/usr/bin/env Rscript
# ============================================================================
# 14_feature_metaprofiles.R   (env: r_g4)   --- Phase 3, point 1 ---
#
# Meta-profiles of G4 CUT&Tag signal centred on NON-promoter genomic features
# (5'UTR, 1st intron, enhancer), split by the topology of the overlapping G4
# peak, across WT and the three KOs. Complements rule 08 (promoter-TSS), which
# is left unchanged.
#
# Each feature region is assigned the topology of its STRONGEST overlapping
# union G4 peak (by max_z, via assign_region_topology()); features with no G4
# peak form the "none" group. Profiles are centred on the feature midpoint and
# strand-aware (mean_replicate_profile flips - strand; enhancers are unstranded).
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           cache_v2 regions_<feature>.rds (config metaprofile.feature_regions
#           keyed into feature.region_files), G4 bigWigs
# Outputs:  results/tables/feature_topology_profiles.csv
#           results/figures/10_metaprofiles/feature_signal_by_topology.{pdf,png}
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
incl_none  <- isTRUE(cfg$metaprofile$include_none_group)
max_group  <- cfg$metaprofile$max_per_group %||% 8000
definite   <- c("parallel", "antiparallel", "hybrid")
keep_groups <- c(definite, if (incl_none) "none")

feature_regions <- as.character(unlist(cfg$metaprofile$feature_regions %||%
                                       c("5UTR", "intron1", "enhancer")))
region_files <- cfg$feature$region_files

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

# --- Build topology-labelled, centre-anchored windows per feature -----------
set.seed(7)
region_windows <- lapply(feature_regions, function(rn) {
  rf <- region_files[[rn]]
  if (is.null(rf)) { message("No region file mapped for '", rn, "'; skipping."); return(NULL) }
  gr <- readRDS(file.path(cfg$paths$cache_v2, rf))
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepSeqlevels(
    gr, intersect(GenomeInfoDb::seqlevels(gr), cfg$std_chroms), pruning.mode = "coarse")
  gr <- assign_region_topology(gr, union)
  gr <- gr[gr$group %in% keep_groups]
  if (length(gr) == 0) return(NULL)
  GenomeInfoDb::seqlengths(gr) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(gr)]
  win <- GenomicRanges::trim(center_window(gr, width = 2 * half_width))
  S4Vectors::mcols(win)$group <- gr$group
  # Cap each topology group for tractable runtime (random, seeded).
  keep_idx <- unlist(lapply(keep_groups, function(g) {
    idx <- which(win$group == g)
    if (length(idx) > max_group) sample(idx, max_group) else idx
  }))
  win <- win[sort(keep_idx)]
  message(rn, " windows per group:"); print(table(win$group))
  win
})
names(region_windows) <- feature_regions
region_windows <- region_windows[!vapply(region_windows, is.null, logical(1))]
stopifnot(length(region_windows) > 0)

# --- Compute per-genotype, per-feature, per-topology meta-profiles ----------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)

# Cache key includes union size, window params AND the feature set so any change
# invalidates it (delete cache/feature_profiles_*.rds if you change max_per_group).
prof_cache <- sprintf("cache/feature_profiles_n%d_hw%d_b%d_%s.rds",
                      length(union), half_width, n_bins,
                      paste(names(region_windows), collapse = "-"))
profiles <- cache_or_build(prof_cache, {
  dplyr::bind_rows(lapply(genotypes, function(g) {
    bws <- bw_meta$filepath[bw_meta$genotype == g]
    if (length(bws) == 0) return(NULL)
    message("Feature profiles: ", g, " (", length(bws), " bigWigs)")
    rows <- dplyr::bind_rows(lapply(names(region_windows), function(rn) {
      win <- region_windows[[rn]]
      dplyr::bind_rows(lapply(keep_groups, function(grp) {
        regions <- win[win$group == grp]
        if (length(regions) < 5) return(NULL)
        prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
        prof$genotype  <- g
        prof$region    <- rn
        prof$group     <- grp
        prof$n_regions <- length(regions)
        prof
      }))
    }))
    evict_bigwigs(bws)
    rows
  }))
})
profiles$genotype <- factor(profiles$genotype, levels = genotypes)
profiles$region   <- factor(profiles$region, levels = feature_regions)
profiles$group    <- factor(profiles$group, levels = keep_groups)
readr::write_csv(profiles, "results/tables/feature_topology_profiles.csv")

# --- Plot: facet_grid(region ~ topology), genotype lines --------------------
region_labels_map <- c("promoter" = "Promoter", "5UTR" = "5'UTR",
                       "intron1" = "1st intron", "enhancer" = "Enhancer (dELS)",
                       "genebody" = "Gene body")
present <- intersect(feature_regions, names(region_windows))   # keep config row order
profiles$region_label <- factor(
  unname(region_labels_map[as.character(profiles$region)]),
  levels = unname(region_labels_map[present]))
group_labels <- c(setNames(topology_labels[match(definite, topology_levels)], definite),
                  none = "No G4 peak")

# Layers shared by both figures (only the facet scaling / title differ).
base_layers <- list(
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA),
  geom_line(linewidth = 0.8),
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40"),
  scale_colour_manual(values = condition_colours),
  scale_fill_manual(values = condition_colours),
  labs(x = "Distance to feature centre (bp)", y = "Mean G4 CUT&Tag signal",
       colour = "Genotype", fill = "Genotype"),
  theme_pub(), theme(legend.position = "bottom"))
fig_h <- 2.2 * length(region_windows) + 1.5

# (1) per-panel free y-axis (each region x topology autoscaled — shows shape).
p_free <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_grid(region_label ~ group, scales = "free_y",
             labeller = labeller(group = group_labels)) +
  labs(title = "Feature-region G4 signal, split by G4 topology",
       subtitle = sprintf("+/- %d bp, %d bins; features grouped by strongest overlapping G4 peak (per-panel y-axis)",
                          half_width, n_bins))
save_plot(p_free, "feature_signal_by_topology", "10_metaprofiles", fig_root,
          width = 11, height = fig_h)

# (2) common y-axis 0..ymax across ALL panels (compare absolute intensities).
ymax <- cfg$metaprofile$feature_fixed_ymax %||% 20
p_fixed <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_grid(region_label ~ group, labeller = labeller(group = group_labels)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(title = "Feature-region G4 signal, split by G4 topology (common y-axis)",
       subtitle = sprintf("+/- %d bp, %d bins; shared y-axis 0-%g for cross-region intensity comparison (lines above %g are clipped)",
                          half_width, n_bins, ymax, ymax))
save_plot(p_fixed, "feature_signal_by_topology_fixedY", "10_metaprofiles", fig_root,
          width = 11, height = fig_h)

message("Done. Wrote feature-region topology meta-profiles (free + fixed y) for: ",
        paste(names(region_windows), collapse = ", "))
