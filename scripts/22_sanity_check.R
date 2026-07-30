#!/usr/bin/env Rscript
# ============================================================================
# 06_sanity_check.R   (env: r_g4)   --- Goal G1 (DISCRIMINATION ONLY) ---
#
# DeepG4 is a *classifier* of G4-formation propensity, not an occupancy
# predictor, so correlating it with continuous CUT&Tag peak strength is weak by
# design (Phase-1 finding: Spearman rho ~ 0.13). This module therefore keeps only
# the question DeepG4 is built to answer:
#   (a) Discrimination: AUROC of DeepG4 separating peaks from GC-matched background.
#   (b) Peak vs background score distributions.
#   (c) DeepG4 probability by assigned topology class.
# The cross-validation of DeepG4 against pqsfinder/G4Hunter, and DeepG4 vs
# helicase-dependence, now live in rule 09 (propensity_metrics).
#
# Inputs:   results/tables/deepg4_scores.csv, peak_catalog.csv, peak_topology.csv
# Outputs:  results/tables/sanity_metrics.csv
#           results/figures/08_topology/{deepg4_roc,deepg4_score_distributions,
#                                      deepg4_by_topology}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(pROC)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/08_topology")
fig_root <- "results/figures"

scores   <- readr::read_csv("results/tables/deepg4_scores.csv", show_col_types = FALSE)
catalog  <- readr::read_csv("results/tables/peak_catalog.csv", show_col_types = FALSE)
topo     <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)

stopifnot(nrow(scores) > 0)
peak_scores <- scores %>% filter(set == "peak") %>%
  left_join(catalog, by = "peak_id") %>%
  left_join(select(topo, peak_id, topology_final), by = "peak_id")

# Supplementary only (expected weak — kept in the metrics table, not plotted):
rho_maxz <- suppressWarnings(cor(peak_scores$deepg4_prob, peak_scores$max_z,
                                 method = "spearman", use = "complete.obs"))
rho_sig  <- suppressWarnings(cor(peak_scores$deepg4_prob, peak_scores$mean_signal,
                                 method = "spearman", use = "complete.obs"))

# --- (a) AUROC peaks vs background -----------------------------------------
auroc <- NA_real_
if (any(scores$set == "background")) {
  roc_obj <- pROC::roc(response = scores$set, predictor = scores$deepg4_prob,
                       levels = c("background", "peak"), direction = "<", quiet = TRUE)
  auroc <- as.numeric(pROC::auc(roc_obj))
  roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                       tpr = roc_obj$sensitivities)
  p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_path(colour = "#D95F02", linewidth = 1) +
    annotate("text", x = 0.6, y = 0.1, size = 5, colour = "#D95F02",
             label = sprintf("AUROC = %.3f", auroc)) +
    coord_equal() +
    labs(x = "False positive rate", y = "True positive rate",
         title = "DeepG4 discriminates G4 peaks from GC-matched background") +
    theme_pub()
  save_plot(p_roc, "deepg4_roc", "08_topology", fig_root, width = 6, height = 6)

  # --- (b) Score distributions, peaks vs background ------------------------
  p_dist <- ggplot(scores, aes(deepg4_prob, fill = set, colour = set)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_fill_manual(values = c(peak = "#1B9E77", background = "grey60")) +
    scale_colour_manual(values = c(peak = "#1B9E77", background = "grey50")) +
    labs(x = "DeepG4 G4-formation probability", y = "Density", fill = NULL, colour = NULL,
         title = "DeepG4 score: G4 peaks vs GC-matched background") +
    theme_pub()
  save_plot(p_dist, "deepg4_score_distributions", "08_topology", fig_root, width = 7, height = 4.5)
}

# --- (c) DeepG4 probability by topology ------------------------------------
if ("topology_final" %in% colnames(peak_scores)) {
  topo_df <- peak_scores %>%
    filter(!is.na(topology_final)) %>%
    mutate(topology_final = factor(topology_final, levels = topology_levels))
  p_topo <- ggplot(topo_df, aes(topology_final, deepg4_prob, fill = topology_final)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, colour = "grey25") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = topology_palette, guide = "none") +
    scale_x_discrete(labels = setNames(topology_labels, topology_levels)) +
    labs(x = NULL, y = "DeepG4 G4-formation probability",
         title = "DeepG4 probability by assigned topology") +
    theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_plot(p_topo, "deepg4_by_topology", "08_topology", fig_root, width = 7.5, height = 5)
}

# --- Metrics table ----------------------------------------------------------
metrics <- data.frame(
  metric = c("auroc_peak_vs_background", "median_prob_peak", "median_prob_background",
             "n_peaks", "n_background",
             "spearman_prob_vs_maxz_SUPPL", "spearman_prob_vs_meansignal_SUPPL"),
  value = c(auroc,
            median(scores$deepg4_prob[scores$set == "peak"], na.rm = TRUE),
            median(scores$deepg4_prob[scores$set == "background"], na.rm = TRUE),
            sum(scores$set == "peak"), sum(scores$set == "background"),
            rho_maxz, rho_sig)
)
readr::write_csv(metrics, "results/tables/sanity_metrics.csv")
message("Done. AUROC = ", round(auroc, 3),
        " | median prob peak/bg = ",
        round(median(scores$deepg4_prob[scores$set == "peak"], na.rm = TRUE), 3), "/",
        round(median(scores$deepg4_prob[scores$set == "background"], na.rm = TRUE), 3))
print(metrics)
