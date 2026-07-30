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
#             deg_motif_pooled,deg_motif_density}.{pdf,png}
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
deg_ymax   <- cfg$metaprofile$deg_fixed_ymax %||% 28
definite   <- c("parallel", "antiparallel", "hybrid")
dir_cols   <- c(up = "#E41A1C", down = "#377EB8")
dir_labels <- c(up = "Upregulated", down = "Downregulated")

pm <- readRDS("cache/promoter_motif_signal.rds")          # rule 30 cache
pm <- pm[pm$topology %in% definite, ]
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
profiles$topology  <- factor(profiles$topology, levels = definite)
profiles$direction <- factor(profiles$direction, levels = c("up", "down"))
readr::write_csv(profiles, "results/tables/deg_promoter_motif_profiles.csv")

topo_labs <- setNames(topology_labels[match(definite, topology_levels)], definite)
sub_base  <- sprintf("motifs within TSS +/- %d bp; significant DEGs (padj<%g & |log2FC|>%g); WT-matched vs KO",
                     hw, deg_padj, deg_lfc)

# --- (A) topology x column grid (shared y 0..deg_ymax) ----------------------
p_grid <- ggplot(profiles, aes(position, mean, colour = direction, fill = direction)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(topology ~ column, labeller = labeller(topology = topo_labs)) +
  coord_cartesian(ylim = c(0, deg_ymax)) +
  scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  labs(x = "Distance to TSS (bp)", y = "Mean per-motif G4 signal",
       title = "DEG promoter G4 signal at individual motifs, by direction and topology",
       subtitle = sprintf("%s; shared y 0-%g", sub_base, deg_ymax)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_grid, "deg_motif_topology_grid", "10_metaprofiles", fig_root, width = 13, height = 8)

# --- (B) pooled (G4+ topologies merged, auto y) -----------------------------
pooled <- profiles %>%
  group_by(column, ko, direction, bin, position) %>%
  summarise(mean = weighted.mean(mean, n_motifs), n_motifs = sum(n_motifs),
            sem = sqrt(sum((sem * n_motifs)^2)) / sum(n_motifs), .groups = "drop")
pooled$column <- factor(pooled$column, levels = col_levels)
p_pool <- ggplot(pooled, aes(position, mean, colour = direction, fill = direction)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ column, nrow = 1) +
  scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  labs(x = "Distance to TSS (bp)", y = "Mean per-motif G4 signal",
       title = "DEG promoter G4 signal at motifs by direction (G4 topologies pooled)",
       subtitle = sprintf("%s; auto y-axis", sub_base)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_pool, "deg_motif_pooled", "10_metaprofiles", fig_root, width = 13, height = 4)

# --- (C) motif density by topology x direction, per KO ----------------------
dens <- bind_rows(dens_rows)
dens$topology  <- factor(dens$topology, levels = definite)
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
