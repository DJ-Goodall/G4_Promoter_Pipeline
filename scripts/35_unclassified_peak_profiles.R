#!/usr/bin/env Rscript
# ============================================================================
# 18_unclassified_peak_profiles.R   (env: r_g4)   --- Phase 7 ---
#
# Characterise the G4 peaks that received NO topology call (topology_final ==
# "no_canonical_PQS": pqsfinder found no canonical PQS, so neither G4SP nor the
# loop heuristic could classify them). These ~11.5% of union peaks are dropped
# from the topology-split metaprofiles, so this rule asks, peak-CENTRED:
#   (A) Signal: what does their G4 CUT&Tag profile look like vs classified peaks,
#       across genotypes? (real binding vs weak/noise?)
#   (B) Propensity: are the underlying sequences G4-like at all? DeepG4 prob +
#       window G4Hunter magnitude vs classified peaks and the GC-matched background.
#   (C) Sequence rule: do they share a compositional signature? GC content,
#       G/C 3+ run count on the richer strand (a canonical G4 needs >=4), and
#       |G-skew| -- testing whether they LACK the G-run scaffold (not G4-capable)
#       or merely fail pqsfinder's loop/score rules (non-canonical but G-rich).
#
# Peaks are unstranded -> strand-agnostic centred windows (no flip). Both G-runs
# (GGG+) and C-runs (CCC+, = G4 on the minus strand) are counted.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           results/tables/propensity_metrics.csv, results/tables/deepg4_scores.csv,
#           cache/peak_seqs_201bp.fa, cache/background_seqs.fa, G4 bigWigs
# Outputs:  results/tables/unclassified_peak_profiles.csv
#           results/tables/unclassified_peak_seqstats.csv
#           results/figures/15_unclassified/unclassified_peak_metaprofile.{pdf,png}
#           results/figures/15_unclassified/unclassified_peak_propensity.{pdf,png}
#           results/figures/15_unclassified/unclassified_peak_composition.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(Biostrings); library(dplyr); library(readr); library(ggplot2); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/15_unclassified")
fig_root <- "results/figures"

genotypes  <- cfg$genotypes
half_width <- cfg$metaprofile$half_width %||% 2000
n_bins     <- cfg$metaprofile$n_bins %||% 80
max_group  <- cfg$metaprofile$max_per_group %||% 8000
unclassified_lvl <- "no_canonical_PQS"

# Two-way class: the no-PQS bucket vs everything with a topology.
class2_levels <- c("unclassified", "classified")
class2_labels <- c(unclassified = "No canonical PQS",
                   classified   = "Classified (P/AP/H)")
group_levels  <- c("No canonical PQS", "Classified (P/AP/H)", "GC-matched background")
group_cols    <- c("No canonical PQS" = "#E41A1C",
                   "Classified (P/AP/H)" = "#377EB8",
                   "GC-matched background" = "grey55")

# --- Peaks + topology -> 2-way class ----------------------------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
union$class2 <- ifelse(union$topology == unclassified_lvl, "unclassified", "classified")
GenomeInfoDb::seqlevelsStyle(union) <- "UCSC"
union <- GenomeInfoDb::keepSeqlevels(
  union, intersect(GenomeInfoDb::seqlevels(union), cfg$std_chroms), pruning.mode = "coarse")
GenomicRanges::strand(union) <- "*"   # peaks are unstranded -> symmetric windows
message("Peak classes: "); print(table(union$class2))

# ============================================================================
# (A) Peak-centred meta-profile: unclassified vs classified, across genotypes
# ============================================================================
GenomeInfoDb::seqlengths(union) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(union)]
win <- GenomicRanges::trim(center_window(union, width = 2 * half_width))
S4Vectors::mcols(win)$class2 <- union$class2
set.seed(7)
keep <- unlist(lapply(class2_levels, function(g) {
  idx <- which(win$class2 == g)
  if (length(idx) > max_group) sample(idx, max_group) else idx
}))
win <- win[sort(keep)]
message("Windows profiled per class (capped at ", max_group, "):"); print(table(win$class2))

bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)
prof_cache <- sprintf("cache/unclassified_peak_profiles_n%d_hw%d_b%d.rds",
                      length(union), half_width, n_bins)
profiles <- cache_or_build(prof_cache, {
  dplyr::bind_rows(lapply(genotypes, function(g) {
    bws <- bw_meta$filepath[bw_meta$genotype == g]
    if (length(bws) == 0) return(NULL)
    message("Peak-centred profiles: ", g, " (", length(bws), " bigWigs)")
    rows <- dplyr::bind_rows(lapply(class2_levels, function(cl) {
      regions <- win[win$class2 == cl]
      if (length(regions) < 5) return(NULL)
      prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
      prof$genotype <- g; prof$class2 <- cl; prof$n_peaks <- length(regions)
      prof
    }))
    evict_bigwigs(bws)
    rows
  }))
})
profiles$genotype     <- factor(profiles$genotype, levels = genotypes)
profiles$class_label  <- factor(unname(class2_labels[profiles$class2]),
                                levels = unname(class2_labels[class2_levels]))
readr::write_csv(profiles, "results/tables/unclassified_peak_profiles.csv")

n_lab <- profiles %>% dplyr::distinct(class_label, genotype, n_peaks) %>%
  dplyr::group_by(class_label) %>% dplyr::summarise(n = n_peaks[1], .groups = "drop") %>%
  dplyr::mutate(label = sprintf("n = %d peaks", n))
p_meta <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.4, size = 3, colour = "grey20") +
  facet_wrap(~ class_label, nrow = 1) +
  scale_colour_manual(values = condition_colours, name = "Genotype") +
  scale_fill_manual(values = condition_colours, name = "Genotype") +
  labs(x = "Distance to peak centre (bp)", y = "Mean G4 CUT&Tag signal",
       title = "G4 CUT&Tag signal at peaks: unclassified (no canonical PQS) vs classified",
       subtitle = sprintf("+/- %d bp, %d bins; peak-centred, strand-agnostic; auto y-axis",
                          half_width, n_bins)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_meta, "unclassified_peak_metaprofile", "15_unclassified", fig_root,
          width = 11, height = 4.5)

# ============================================================================
# (B + C) Sequence propensity and composition
# ============================================================================
peak_seqs <- Biostrings::readDNAStringSet("cache/peak_seqs_201bp.fa")
bg_seqs   <- Biostrings::readDNAStringSet("cache/background_seqs.fa")

# Per-sequence composition: GC, G/C 3+ run count (richer strand), |G-skew|.
seq_compose <- function(dss) {
  s  <- as.character(dss)
  lf <- Biostrings::letterFrequency(dss, c("G", "C"))
  G  <- lf[, "G"]; C <- lf[, "C"]; w <- Biostrings::width(dss)
  data.frame(
    gc      = (G + C) / w,
    g4_runs = pmax(stringr::str_count(s, "G{3,}"), stringr::str_count(s, "C{3,}")),
    g_skew  = (G - C) / pmax(G + C, 1),
    stringsAsFactors = FALSE)
}
peak_comp <- seq_compose(peak_seqs); peak_comp$peak_id <- names(peak_seqs)
bg_comp   <- seq_compose(bg_seqs)

peak_comp$class2 <- union$class2[match(peak_comp$peak_id, names(union))]
peak_comp <- peak_comp[!is.na(peak_comp$class2), ]
peak_comp$group <- unname(class2_labels[peak_comp$class2])
bg_comp$group   <- "GC-matched background"

comp <- dplyr::bind_rows(
  peak_comp[, c("group", "gc", "g4_runs", "g_skew")],
  bg_comp[,   c("group", "gc", "g4_runs", "g_skew")])
comp$group <- factor(comp$group, levels = group_levels)

# Per-group summary table (the "rule": % structurally G4-capable = >=4 G/C runs).
seqstats <- comp %>% dplyr::group_by(group) %>% dplyr::summarise(
  n               = dplyr::n(),
  gc_med          = round(median(gc), 3),
  g4_runs_med     = median(g4_runs),
  pct_ge4_runs    = round(100 * mean(g4_runs >= 4), 1),
  abs_gskew_med   = round(median(abs(g_skew)), 3),
  .groups = "drop")
readr::write_csv(seqstats, "results/tables/unclassified_peak_seqstats.csv")
message("Sequence composition by group:"); print(as.data.frame(seqstats))

# (C) composition figure -- GC, G/C runs, |G-skew| by group
comp_long <- dplyr::bind_rows(
  data.frame(group = comp$group, metric = "GC content",                 value = comp$gc),
  data.frame(group = comp$group, metric = "G/C 3+ runs (richer strand)", value = comp$g4_runs),
  data.frame(group = comp$group, metric = "|G-skew|",                    value = abs(comp$g_skew)))
comp_long$metric <- factor(comp_long$metric,
  levels = c("GC content", "G/C 3+ runs (richer strand)", "|G-skew|"))
p_comp <- ggplot(comp_long, aes(group, value, fill = group)) +
  geom_violin(scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Value",
       title = "Sequence composition: do unclassified peaks share a rule?",
       subtitle = "A canonical G4 needs >=4 runs of >=3 G (or C, minus strand); see pct_ge4_runs in seqstats CSV") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_comp, "unclassified_peak_composition", "15_unclassified", fig_root,
          width = 11, height = 4.5)

# (B) propensity figure -- DeepG4 prob + window G4Hunter magnitude by group
prop <- readr::read_csv("results/tables/propensity_metrics.csv", show_col_types = FALSE)
prop$class2 <- union$class2[match(prop$peak_id, names(union))]
prop <- prop[!is.na(prop$class2), ]
prop$group <- unname(class2_labels[prop$class2])
dg  <- readr::read_csv("results/tables/deepg4_scores.csv", show_col_types = FALSE)
bg_dg <- dg$deepg4_prob[dg$set == "background"]
bg_g4h <- abs(g4hunter_score(bg_seqs))   # same algorithm rule 09 used for peaks

prop_long <- dplyr::bind_rows(
  data.frame(group = prop$group,                 metric = "DeepG4 formation prob", value = prop$deepg4_prob),
  data.frame(group = "GC-matched background",    metric = "DeepG4 formation prob", value = bg_dg),
  data.frame(group = prop$group,                 metric = "|G4Hunter| (window)",   value = abs(prop$g4hunter_win)),
  data.frame(group = "GC-matched background",    metric = "|G4Hunter| (window)",   value = bg_g4h))
prop_long$group  <- factor(prop_long$group, levels = group_levels)
prop_long$metric <- factor(prop_long$metric, levels = c("DeepG4 formation prob", "|G4Hunter| (window)"))
p_prop <- ggplot(prop_long, aes(group, value, fill = group)) +
  geom_violin(scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Value",
       title = "Are unclassified peaks G4-like by sequence?",
       subtitle = "DeepG4 + window G4Hunter magnitude vs classified peaks and GC-matched background") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_prop, "unclassified_peak_propensity", "15_unclassified", fig_root,
          width = 9, height = 4.5)

message("Done. Unclassified-peak profiles + sequence characterisation.")
