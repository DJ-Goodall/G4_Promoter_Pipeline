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
# ---------------------------------------------------------------------------
# v1.1 (2026-08-07 audit) -- two changes, both about the "mixed" stratum.
#
# 1. "mixed" is now DISPLAYED, not dropped. Rule 21 no longer back-fills
#    G4ShapePredictor's abstentions with the loop-length heuristic, so at
#    precision 0.80 roughly a fifth of PQS-bearing peaks carry the terminal
#    label "mixed". v1.0's keep_groups filter would have silently DELETED every
#    promoter whose strongest peak was mixed -- they would not have moved to
#    "none", they would have vanished from the denominator. Showing the facet
#    keeps the figure honest about how much of the promoter set is unresolved.
#
# 2. The script is parametrised so the same code renders the precision
#    sensitivity arm (see 21b_topology_sensitivity.R). This is the figure the
#    audit says must be re-derived: v1.0's "hybrid > parallel > antiparallel"
#    ordering rested on 14,429 hybrids of which 87% were heuristic back-fills,
#    and the primary v1.1 arm leaves only 1,879. Reporting it at one precision
#    would be reporting a threshold. Both arms are rendered.
#
# Usage:  Rscript scripts/24_tss_metaprofiles.R [config] [topology_csv] [suffix]
#   topology_csv  default results/tables/peak_topology.csv  (the primary arm)
#   suffix        appended to all four output names, e.g. "_argmax"
# ---------------------------------------------------------------------------
#
# Inputs:   cache/promoters.rds, cache/peaks.rds, results/tables/peak_topology.csv, G4 bigWigs
# Outputs:  results/tables/tss_topology_profiles<suffix>.csv
#           results/tables/tss_topology_profiles_norm<suffix>.csv
#           results/figures/10_metaprofiles/tss_signal_by_topology<suffix>.{pdf,png}
#           results/figures/10_metaprofiles/tss_signal_by_topology_norm<suffix>.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
# v1.1: optional topology table + output suffix, so the precision sensitivity
# arm reuses this script rather than duplicating it. Defaults reproduce v1.0.
topo_path   <- if (length(args) >= 2 && nzchar(args[[2]])) args[[2]] else
  "results/tables/peak_topology.csv"
out_suffix  <- if (length(args) >= 3) args[[3]] else ""
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles")
fig_root <- "results/figures"

genotypes  <- cfg$genotypes
half_width <- cfg$metaprofile$half_width
n_bins     <- cfg$metaprofile$n_bins
incl_none  <- isTRUE(cfg$metaprofile$include_none_group)
max_group  <- cfg$metaprofile$max_per_group %||% 8000
definite   <- topology_definite            # v1.1: shared, was a local copy

message("Topology labels from: ", topo_path,
        if (nzchar(out_suffix)) paste0("  (outputs suffixed '", out_suffix, "')") else "")
promoters    <- readRDS("cache/promoters.rds")
peaks_obj    <- readRDS("cache/peaks.rds")
union        <- peaks_obj$union
topo         <- readr::read_csv(topo_path, show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

# --- Assign each promoter the topology of its strongest overlapping G4 peak --
# (shared rule; see assign_region_topology() in helpers.R)
promoters <- assign_region_topology(promoters, union)

# Keep only the topology facets we want to display.
#
# v1.1: "mixed" is included WHEN PRESENT. Under v1.0 it never was, because rule
# 21 back-filled every abstention. Under v1.1 it is a real stratum, and
# dropping it here would delete those promoters from the figure without them
# appearing in "none" -- a silent change of denominator, which is the class of
# error the audit was about. If G4SP abstained on the strongest peak at a
# promoter, the honest statement is "unresolved", not "absent".
keep_groups <- topology_keep_groups(promoters$group, incl_none)
has_mixed   <- "mixed" %in% keep_groups
n_before <- length(promoters)
promoters <- promoters[!is.na(promoters$group) & promoters$group %in% keep_groups]
promoters$group <- factor(promoters$group, levels = keep_groups)
message(sprintf("Promoters retained: %d / %d (groups: %s)",
                length(promoters), n_before, paste(keep_groups, collapse = ", ")))
if (has_mixed)
  message(sprintf("  of which 'mixed' (G4SP declined to call): %d (%.1f%%)",
                  sum(promoters$group == "mixed"),
                  100 * mean(promoters$group == "mixed")))

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
readr::write_csv(profiles,
                 sprintf("results/tables/tss_topology_profiles%s.csv", out_suffix))

# --- Plot -------------------------------------------------------------------
# Labels drawn from keep_groups so a newly present "mixed" facet is labelled
# from topology_labels automatically rather than showing its raw level name.
group_labels <- topology_group_labels(keep_groups)
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
       subtitle = sprintf(
         "+/- %d bp, %d bins; promoters grouped by strongest overlapping G4 peak\nTopology from %s%s",
         half_width, n_bins, basename(topo_path),
         if (has_mixed) "; 'mixed' = G4ShapePredictor declined to call" else "")) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p, paste0("tss_signal_by_topology", out_suffix), "10_metaprofiles",
          fig_root, width = 11, height = 7)

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
  readr::write_csv(
    norm, sprintf("results/tables/tss_topology_profiles_norm%s.csv", out_suffix))

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
  save_plot(p_norm, paste0("tss_signal_by_topology_norm", out_suffix),
            "10_metaprofiles", fig_root, width = 11, height = 5)
  message("Wrote background-normalized promoter meta-profiles.")
}

message("Done. Wrote TSS topology meta-profiles.")
