#!/usr/bin/env Rscript
# ============================================================================
# 08_tss_metaprofiles.R   (env: r_g4)   --- Goal G4 ---
#
# Meta-profiles of G4 CUT&Tag signal centred on promoter TSS, split by the
# topology of the overlapping G4 peak, across WT and the three KOs. Shows where
# along the promoter each topology is enriched and how KO changes it.
#
# Each promoter is assigned the topology of the STRONGEST union G4 peak it
# overlaps (by max_z); promoters with no G4 peak form the optional "none" group.
#
# Inputs:   cache/promoters.rds, cache/peaks.rds, results/tables/peak_topology.csv, G4 bigWigs
# Outputs:  results/tables/tss_topology_profiles.csv
#           results/figures/10_metaprofiles/tss_signal_by_topology.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles")
fig_root <- "results/figures"

genotypes  <- cfg$genotypes
half_width <- cfg$metaprofile$half_width
n_bins     <- cfg$metaprofile$n_bins
incl_none  <- isTRUE(cfg$metaprofile$include_none_group)
max_group  <- cfg$metaprofile$max_per_group %||% 8000
definite   <- c("parallel", "antiparallel", "hybrid")

promoters    <- readRDS("cache/promoters.rds")
peaks_obj    <- readRDS("cache/peaks.rds")
union        <- peaks_obj$union
topo         <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

# --- Assign each promoter the topology of its strongest overlapping G4 peak --
# (shared rule; see assign_region_topology() in helpers.R)
promoters <- assign_region_topology(promoters, union)

# Keep only the topology facets we want to display.
keep_groups <- c(definite, if (incl_none) "none")
promoters <- promoters[!is.na(promoters$group) & promoters$group %in% keep_groups]
promoters$group <- factor(promoters$group, levels = keep_groups)

# TSS-centred windows (+/- half_width). Promoters are symmetric (TSS +/- 2kb) so
# their midpoint is the TSS; re-centre and set the requested half-width.
GenomeInfoDb::seqlengths(promoters) <- default_chrom_sizes(cfg$std_chroms)[
  GenomeInfoDb::seqlevels(promoters)]
tss_windows <- center_window(promoters, width = 2 * half_width)
tss_windows <- GenomicRanges::trim(tss_windows)
S4Vectors::mcols(tss_windows)$group <- promoters$group

# Cap group sizes for tractable runtime (random, seeded).
set.seed(7)
keep_idx <- unlist(lapply(keep_groups, function(g) {
  idx <- which(tss_windows$group == g)
  if (length(idx) > max_group) sample(idx, max_group) else idx
}))
tss_windows <- tss_windows[sort(keep_idx)]
message("Promoter windows per group:")
print(table(tss_windows$group))

# --- Compute per-genotype, per-topology meta-profiles ----------------------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)

profiles <- dplyr::bind_rows(lapply(genotypes, function(g) {
  bws <- bw_meta$filepath[bw_meta$genotype == g]
  if (length(bws) == 0) return(NULL)
  message("TSS profiles: ", g, " (", length(bws), " bigWigs)")
  rows <- dplyr::bind_rows(lapply(keep_groups, function(grp) {
    regions <- tss_windows[tss_windows$group == grp]
    if (length(regions) < 5) return(NULL)
    prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
    prof$genotype <- g; prof$group <- grp; prof
  }))
  evict_bigwigs(bws)
  rows
}))
profiles$genotype <- factor(profiles$genotype, levels = genotypes)
profiles$group    <- factor(profiles$group, levels = keep_groups)
readr::write_csv(profiles, "results/tables/tss_topology_profiles.csv")

# --- Plot -------------------------------------------------------------------
group_labels <- c(setNames(topology_labels[match(definite, topology_levels)], definite),
                  none = "No G4 peak")
p <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ group, scales = "free_y", labeller = labeller(group = group_labels)) +
  scale_colour_manual(values = condition_colours) +
  scale_fill_manual(values = condition_colours) +
  labs(x = "Distance to TSS (bp)", y = "Mean G4 CUT&Tag signal",
       colour = "Genotype", fill = "Genotype",
       title = "Promoter G4 signal around TSS, split by G4 topology",
       subtitle = sprintf("+/- %d bp, %d bins; promoters grouped by strongest overlapping G4 peak",
                          half_width, n_bins)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p, "tss_signal_by_topology", "10_metaprofiles", fig_root, width = 11, height = 7)

# --- Background-normalized companion (promoter only) ------------------------
# In KO cells the "No G4 peak" promoters also show elevated TSS-proximal signal
# (general sample background, since the models predict no G4 there). Divide each
# definite-topology bin by the SAME-genotype "none" group bin to remove that
# per-genotype baseline and reveal topology-specific, fold-over-background signal.
if (isTRUE(cfg$metaprofile$normalize_promoter_to_none %||% TRUE) && incl_none &&
    "none" %in% levels(profiles$group)) {
  eps <- 1e-6
  none_df <- profiles %>% dplyr::filter(group == "none") %>%
    dplyr::select(genotype, bin, none_mean = mean, none_sem = sem)
  norm <- profiles %>% dplyr::filter(group %in% definite) %>%
    dplyr::left_join(none_df, by = c("genotype", "bin")) %>%
    dplyr::mutate(
      ratio     = mean / pmax(none_mean, eps),
      # ratio-of-means error propagation (independent relative errors)
      ratio_sem = ratio * sqrt((sem / pmax(mean, eps))^2 +
                               (none_sem / pmax(none_mean, eps))^2))
  norm$group    <- factor(as.character(norm$group), levels = definite)
  norm$genotype <- factor(as.character(norm$genotype), levels = genotypes)
  readr::write_csv(norm, "results/tables/tss_topology_profiles_norm.csv")

  norm_labels <- setNames(topology_labels[match(definite, topology_levels)], definite)
  p_norm <- ggplot(norm, aes(position, ratio, colour = genotype, fill = genotype)) +
    geom_ribbon(aes(ymin = ratio - ratio_sem, ymax = ratio + ratio_sem),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ group, scales = "free_y", labeller = labeller(group = norm_labels)) +
    scale_colour_manual(values = condition_colours) +
    scale_fill_manual(values = condition_colours) +
    labs(x = "Distance to TSS (bp)",
         y = "G4 signal / no-G4-peak signal (matched genotype)",
         colour = "Genotype", fill = "Genotype",
         title = "Promoter G4 signal relative to no-G4-peak background, by topology",
         subtitle = sprintf("Each bin / matched-genotype no-G4-peak signal; +/- %d bp, %d bins",
                            half_width, n_bins)) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p_norm, "tss_signal_by_topology_norm", "10_metaprofiles", fig_root,
            width = 11, height = 5)
  message("Wrote background-normalized promoter meta-profiles.")
}

message("Done. Wrote TSS topology meta-profiles.")
