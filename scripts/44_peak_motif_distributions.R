#!/usr/bin/env Rscript
# ============================================================================
# 29_peak_motif_distributions.R   (env: r_g4)   --- Phase 13, Analysis 1 ---
#
# How big are G4 peaks, how many distinct G4 motifs do they contain, and does
# counting EVERY motif (rule 27) change the topology picture vs the single-top-PQS
# peak-level call (rule 05)? This quantifies the size mismatch that motivates the
# whole all-PQS reframe: a wide peak can hold several G4s of different topologies.
#
# Inputs:   cache/peaks.rds (union widths), cache/motifs_all.rds (+ $topology),
#           results/tables/motif_topology.csv, results/tables/peak_topology.csv
# Outputs:  results/tables/peak_motif_summary.csv, motif_count_distribution.csv,
#           topology_composition_compare.csv
#           results/figures/20_motif_distributions/{peak_width_distribution,
#             motifs_per_peak,motifs_vs_width,topology_composition_motif_vs_peak}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(S4Vectors); library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/20_motif_distributions")
fig_root <- "results/figures"
definite <- c("parallel", "antiparallel", "hybrid")

union     <- readRDS("cache/peaks.rds")$union
motifs_gr <- readRDS("cache/motifs_all.rds")
motif_topo <- as.character(S4Vectors::mcols(motifs_gr)$topology)
motif_peak <- as.character(S4Vectors::mcols(motifs_gr)$peak_id)
peak_topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)

# --- Per-peak width + motif count -------------------------------------------
n_motif_tab <- table(motif_peak)
peak_df <- data.frame(
  peak_id  = names(union),
  width    = as.integer(GenomicRanges::width(union)),
  n_motifs = as.integer(n_motif_tab[match(names(union), names(n_motif_tab))]),
  stringsAsFactors = FALSE)
peak_df$n_motifs[is.na(peak_df$n_motifs)] <- 0L
peak_df$has_pqs <- peak_df$n_motifs > 0
peak_df$motif_bucket <- cut(peak_df$n_motifs, breaks = c(-1, 0, 1, 2, 3, Inf),
                            labels = c("0", "1", "2", "3", "4+"))
readr::write_csv(peak_df, "results/tables/peak_motif_summary.csv")

count_dist <- peak_df %>% count(n_motifs, name = "n_peaks")
readr::write_csv(count_dist, "results/tables/motif_count_distribution.csv")

# --- (1) Peak-width distribution (overall + by motif-count bucket) -----------
p_w <- ggplot(peak_df, aes(width)) +
  geom_histogram(bins = 60, fill = "#377EB8", colour = "white", linewidth = 0.1) +
  scale_x_log10() +
  labs(x = "Union peak width (bp, log10)", y = "Peaks",
       title = "G4 peak width distribution",
       subtitle = sprintf("n = %d union peaks; median %d bp, mean %.0f bp, max %d bp",
                          nrow(peak_df), median(peak_df$width),
                          mean(peak_df$width), max(peak_df$width))) +
  theme_pub()
save_plot(p_w, "peak_width_distribution", "20_motif_distributions", fig_root, width = 7, height = 5)

p_wb <- ggplot(peak_df %>% filter(n_motifs > 0), aes(width, fill = motif_bucket)) +
  geom_histogram(bins = 50, colour = "white", linewidth = 0.05) +
  facet_wrap(~ motif_bucket, scales = "free_y", nrow = 1) +
  scale_x_log10() + scale_fill_brewer(palette = "YlOrRd", guide = "none") +
  labs(x = "Union peak width (bp, log10)", y = "Peaks",
       title = "Peak width by number of distinct G4 motifs",
       subtitle = "Wider peaks carry more motifs (collapsed, score>=20)") +
  theme_pub()
save_plot(p_wb, "peak_width_by_motif_count", "20_motif_distributions", fig_root, width = 12, height = 4)

# --- (2) Motifs per peak -----------------------------------------------------
mc <- peak_df %>% mutate(nm = pmin(n_motifs, 10)) %>% count(nm, name = "n_peaks")
p_m <- ggplot(mc, aes(factor(nm), n_peaks)) +
  geom_col(fill = "#1B9E77", colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n_peaks), vjust = -0.3, size = 2.8) +
  scale_x_discrete(labels = function(x) ifelse(x == "10", "10+", x)) +
  labs(x = "Distinct G4 motifs in peak", y = "Peaks",
       title = "Number of G4 motifs per peak",
       subtitle = sprintf("%.1f%% of peaks have >=2 motifs; mean %.2f motifs/peak (PQS+ peaks)",
                          100 * mean(peak_df$n_motifs >= 2),
                          mean(peak_df$n_motifs[peak_df$n_motifs > 0]))) +
  theme_pub()
save_plot(p_m, "motifs_per_peak", "20_motif_distributions", fig_root, width = 7, height = 5)

# --- (3) Motif count vs peak width ------------------------------------------
p_mw <- ggplot(peak_df %>% filter(n_motifs > 0), aes(width, n_motifs)) +
  geom_hex(bins = 40) + scale_x_log10() +
  scale_fill_viridis_c(trans = "log10", name = "Peaks") +
  labs(x = "Union peak width (bp, log10)", y = "Distinct G4 motifs",
       title = "More motifs in wider peaks",
       subtitle = "Each wide peak is currently labelled by ONE motif's topology (rules 05-26)") +
  theme_pub()
save_plot(p_mw, "motifs_vs_width", "20_motif_distributions", fig_root, width = 7, height = 5)

# --- (4) Topology composition: motif-level vs peak-level --------------------
motif_comp <- data.frame(level = "Per motif (all PQS)",
                         topology = motif_topo, stringsAsFactors = FALSE) %>%
  filter(!is.na(topology)) %>% count(level, topology, name = "n")
peak_comp <- data.frame(level = "Per peak (top PQS)",
                        topology = as.character(peak_topo$topology_final),
                        stringsAsFactors = FALSE) %>%
  filter(!is.na(topology)) %>% count(level, topology, name = "n")
comp <- bind_rows(motif_comp, peak_comp) %>%
  group_by(level) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
comp$topology <- factor(comp$topology,
                        levels = c(definite, "no_canonical_PQS", "unclassified"))
readr::write_csv(comp, "results/tables/topology_composition_compare.csv")

topo_pal <- c(topology_palette, unclassified = "grey50")
p_c <- ggplot(comp, aes(level, pct, fill = topology)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = topo_pal, name = "Topology",
                    labels = c(setNames(topology_labels, topology_levels),
                               unclassified = "Unclassified")) +
  labs(x = NULL, y = "Share of calls (%)",
       title = "G4 topology composition: motif-level vs peak-level",
       subtitle = "Does counting every motif shift the parallel-dominated peak-level picture?") +
  theme_pub()
save_plot(p_c, "topology_composition_motif_vs_peak", "20_motif_distributions",
          fig_root, width = 7, height = 6)

message("Done. Analysis 1 (peak/motif distributions).")
print(comp)
