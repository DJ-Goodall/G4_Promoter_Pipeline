#!/usr/bin/env Rscript
# ============================================================================
# 31_deg_promoter_motif_metaprofiles.R  (env: r_g4)  --- Phase 13, Analysis 3 ---
#
# The rule-15 DEG promoter metaprofile, re-done at MOTIF resolution. Instead of one
# topology per promoter, every G4 motif within a DEG promoter's TSS +/- 2 kb is a
# point with its OWN topology, signed distance to TSS, and per-genotype signal
# (reused from rule 30's cache). Split by DEG direction (up=red/down=blue) and G4
# topology, per KO, with a matched-WT control column (mirrors rule 15 layout):
#   [WT (DHX36KO DEGs)] [DHX36KO] [WT (dKO DEGs)] [dKO]
#
# Inputs:   cache/promoter_motif_signal.rds (rule 30), results/tables/deseq2_KO_vs_WT.csv
# Outputs:  results/tables/deg_promoter_motif_profiles.csv
#           results/figures/10_metaprofiles/{deg_motif_topology_grid,
#             deg_motif_pooled,deg_motif_pooled_unfiltered,
#             deg_motif_density}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles")
fig_root <- "results/figures"

ref_geno   <- cfg$ref_genotype
mcfg       <- cfg$motif_analysis
hw         <- mcfg$promoter_half_width %||% 2000
n_bins     <- mcfg$n_bins %||% cfg$metaprofile$n_bins %||% 80
grid_kos   <- as.character(unlist(cfg$metaprofile$deg_grid_kos %||% c("DHX36KO", "dKO")))
deg_padj   <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc    <- cfg$metaprofile$deg_lfc %||% 0.5
# NB: metaprofile.deg_fixed_ymax (28) is deliberately NOT used here. It caps rule
# 31's window-AVERAGED promoter profile; this script plots a per-motif POINT
# statistic (motif centre +/- signal_half_width) whose values reach ~395, so the
# same cap clipped most of the grid into flat-topped noise. Rows autoscale instead.
definite   <- topology_definite            # v1.1: shared, was a local copy
dir_cols   <- c(up = "#E41A1C", down = "#377EB8")
dir_labels <- c(up = "Upregulated", down = "Downregulated")

pm <- readRDS("cache/promoter_motif_signal.rds")          # rule 45 cache
# v1.1 (2026-08-11): definite-only dropped 23.5% of promoter motifs once rule 43
# stopped back-filling the G4SP abstentions -- and dropped them from the
# denominator too, so the figure could not show that they existed. "mixed" is now
# its own row, as in rules 24/30/33. incl_none = FALSE: every motif IS a PQS, so
# there is no "no G4 peak" class here.
keep_groups <- topology_keep_groups(pm$topology, incl_none = FALSE)
n_pre <- nrow(pm)
pm <- pm[pm$topology %in% keep_groups, ]
message(sprintf(
  "Motifs profiled: %d / %d (%.1f%% excluded as unclassified)",
  nrow(pm), n_pre, 100 * (1 - nrow(pm) / max(n_pre, 1))))
message("Motifs per topology:"); print(table(pm$topology))
pm$symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = as.character(pm$prom_entrez),
                                   column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
pm <- pm[!is.na(pm$symbol), ]

breaks  <- seq(-hw, hw, length.out = n_bins + 1)
centres <- (breaks[-1] + breaks[-(n_bins + 1)]) / 2
pm$bin  <- cut(pm$dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
pm <- pm[!is.na(pm$bin), ]; pm$position <- centres[pm$bin]

de <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)

all_rows <- list(); dens_rows <- list()
for (k in grid_kos) {
  dek <- de %>% dplyr::filter(ko == k, !is.na(padj), !is.na(log2FC))
  dek$direction <- dplyr::case_when(
    dek$padj < deg_padj & dek$log2FC >  deg_lfc ~ "up",
    dek$padj < deg_padj & dek$log2FC < -deg_lfc ~ "down",
    TRUE ~ NA_character_)
  dek <- dek %>% dplyr::filter(!is.na(direction)) %>% dplyr::distinct(gene_name, .keep_all = TRUE)
  dir_map <- setNames(dek$direction, dek$gene_name)

  pmk <- pm; pmk$direction <- unname(dir_map[pmk$symbol]); pmk <- pmk[!is.na(pmk$direction), ]
  if (nrow(pmk) == 0) { message("No DEG-promoter motifs for ", k); next }
  message(sprintf("DEG-promoter motifs: %s; %d up / %d down motifs",
                  k, sum(pmk$direction == "up"), sum(pmk$direction == "down")))

  # motif density per (topology, direction, bin) -- genotype-independent
  dens_rows[[k]] <- pmk %>% count(topology, direction, bin, position, name = "n_motifs") %>%
    mutate(ko = k)

  # signal: matched WT control column then the KO's own column
  for (sg in c(ref_geno, k)) {
    col <- if (sg == ref_geno) sprintf("WT (%s DEGs)", k) else k
    val <- pmk[[paste0("signal_", sg)]]
    all_rows[[length(all_rows) + 1]] <- pmk %>%
      mutate(column = col, ko = k, signal = val) %>%
      group_by(column, ko, topology, direction, bin, position) %>%
      summarise(n_motifs = dplyr::n(),
                mean = mean(signal, na.rm = TRUE),
                sem  = sd(signal, na.rm = TRUE) / sqrt(sum(!is.na(signal))),
                .groups = "drop")
  }
}
profiles <- bind_rows(all_rows)
if (nrow(profiles) == 0) stop("No DEG-promoter motif profiles computed.")

col_levels <- unlist(lapply(grid_kos, function(g) c(sprintf("WT (%s DEGs)", g), g)))
profiles$column    <- factor(profiles$column, levels = col_levels)
profiles$topology  <- factor(profiles$topology, levels = keep_groups)
profiles$direction <- factor(profiles$direction, levels = c("up", "down"))
readr::write_csv(profiles, "results/tables/deg_promoter_motif_profiles.csv")

topo_labs <- topology_group_labels(keep_groups)
sub_base  <- sprintf("motifs within TSS +/- %d bp; significant DEGs (padj<%g & |log2FC|>%g); WT-matched vs KO",
                     hw, deg_padj, deg_lfc)

# Total motifs behind each panel. Each motif contributes to exactly ONE of the
# n_bins bins (unlike rule 31, where every promoter contributes to every bin), so
# a panel's n divided by n_bins is roughly its per-bin support -- which is the
# number that decides whether a wiggle means anything. Annotate it.
strat_n <- profiles %>%
  group_by(column, topology, direction) %>%
  summarise(n = sum(n_motifs), .groups = "drop")

# --- (A) topology x column grid (y autoscaled per topology row) -------------
p_grid <- ggplot(profiles, aes(position, mean, colour = direction, fill = direction)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(data = panel_n_labels(strat_n, c("topology", "column"), "n"),
            aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 2.8, colour = "grey20") +
  facet_grid(topology ~ column, scales = "free_y", drop = FALSE,
             labeller = labeller(topology = topo_labs)) +
  scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  labs(x = "Distance to TSS (bp)", y = "Mean per-motif G4 signal",
       title = "DEG promoter G4 signal at individual motifs, by direction and topology",
       subtitle = sprintf(paste("%s; y autoscaled per topology row (magnitudes differ ~10x);",
                                "n = motifs behind each panel, spread over %d bins"),
                          sub_base, n_bins)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_grid, "deg_motif_topology_grid", "10_metaprofiles", fig_root, width = 13, height = 8)

# --- (B) pooled, in two arms exactly as rule 31 does it ---------------------
#   "pooled"     -> the three definite classes only (weighted average of the grid
#                   rows), so the number stays comparable with rule 31's pooled panel.
#   "unfiltered" -> every motif, mixed included. Without this arm the only pooled
#                   motif figure would still discard 23.5% of the data silently,
#                   which is the defect class the 2026-08-07 audit raised.
# na.rm throughout: a bin holding a single motif has sd = NA, and without it one
# such bin in any topology nulled that bin's pooled mean and ribbon entirely.
build_pooled <- function(df, fname, title, subtitle) {
  pooled <- df %>%
    group_by(column, ko, direction, bin, position) %>%
    summarise(mean = weighted.mean(mean, n_motifs, na.rm = TRUE),
              n_motifs = sum(n_motifs),
              sem = sqrt(sum((sem * n_motifs)^2, na.rm = TRUE)) / sum(n_motifs),
              .groups = "drop")
  pooled$column <- factor(pooled$column, levels = col_levels)
  panel_n <- pooled %>% group_by(column, direction) %>%
    summarise(n = sum(n_motifs), .groups = "drop")
  p <- ggplot(pooled, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = panel_n_labels(panel_n, "column", "n"),
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey20") +
    facet_wrap(~ column, nrow = 1) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = "Distance to TSS (bp)", y = "Mean per-motif G4 signal",
         title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, "10_metaprofiles", fig_root, width = 13, height = 4)
}
build_pooled(profiles %>% dplyr::filter(topology %in% definite), "deg_motif_pooled",
             "DEG promoter G4 signal at motifs by direction (G4 topologies pooled)",
             sprintf("%s; parallel+antiparallel+hybrid pooled (mixed excluded); auto y-axis", sub_base))
build_pooled(profiles, "deg_motif_pooled_unfiltered",
             "DEG promoter G4 signal at motifs by direction (all motifs - no topology filter)",
             sprintf("%s; ALL motifs incl. mixed/unresolved; auto y-axis", sub_base))

# --- (C) motif density by topology x direction, per KO ----------------------
dens <- bind_rows(dens_rows)
dens$topology  <- factor(dens$topology, levels = keep_groups)
dens$direction <- factor(dens$direction, levels = c("up", "down"))
dens$ko        <- factor(dens$ko, levels = grid_kos)
p_den <- ggplot(dens, aes(position, n_motifs, colour = topology)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(ko ~ direction, labeller = labeller(direction = dir_labels)) +
  scale_colour_manual(values = topology_palette, labels = topo_labs, name = "Topology") +
  labs(x = "Distance to TSS (bp)", y = "G4 motifs (count per bin)",
       title = "DEG promoter G4 motif density by topology and DEG direction",
       subtitle = sprintf("motifs within TSS +/- %d bp, %d bins", hw, n_bins)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_den, "deg_motif_density", "10_metaprofiles", fig_root, width = 9, height = 6)

message("Done. Analysis 3 (DEG promoter motif metaprofiles).")
