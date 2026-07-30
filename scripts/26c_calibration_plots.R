#!/usr/bin/env Rscript
# ============================================================================
# 10c_calibration_plots.R   (env: r_g4)   --- topology calibration (Q3) ---
#
# Visualise the calibration panel:
#   (1) precision-threshold sweep  -> is the parallel fraction threshold-stable
#       (real biology) or does it collapse (artifact)?
#   (2) per-class G4SP probability distributions.
#   (3) known-topology control benchmark confusion (reported vs predicted).
#
# Inputs:   results/tables/topology_precision_sweep.csv,
#           results/tables/g4sp_topology.csv,
#           results/tables/g4sp_control_benchmark.csv
# Outputs:  results/figures/08_topology/{calibration_precision_sweep,
#           calibration_proba_distributions,calibration_control_benchmark}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures/08_topology")
fig_root <- "results/figures"

class_levels <- c("parallel", "antiparallel", "hybrid", "mixed")
class_pal <- c(parallel = "#1B9E77", antiparallel = "#D95F02",
               hybrid = "#7570B3", mixed = "grey60")

# --- (1) Precision sweep ----------------------------------------------------
sweep <- readr::read_csv("results/tables/topology_precision_sweep.csv",
                         show_col_types = FALSE)
sweep_long <- sweep %>%
  select(precision, starts_with("frac_")) %>%
  pivot_longer(-precision, names_to = "class", values_to = "frac") %>%
  mutate(class = factor(sub("frac_", "", class), levels = class_levels),
         precision = factor(precision, levels = sweep$precision))
p_sweep <- ggplot(sweep_long, aes(precision, frac, colour = class, group = class)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = class_pal, name = "Topology") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = "G4SP precision setting", y = "Fraction of peaks",
       title = "Topology fractions vs G4SP precision threshold",
       subtitle = "Stable parallel dominance = real; collapse with precision = threshold artifact") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_sweep, "calibration_precision_sweep", "08_topology", fig_root,
          width = 8, height = 5)

# --- (2) Per-class probability distributions --------------------------------
g4sp <- readr::read_csv("results/tables/g4sp_topology.csv", show_col_types = FALSE)
proba_long <- g4sp %>%
  filter(!is.na(p_parallel)) %>%
  select(peak_id, p_parallel, p_antiparallel, p_hybrid) %>%
  pivot_longer(-peak_id, names_to = "class", values_to = "prob") %>%
  mutate(class = factor(sub("p_", "", class),
                        levels = c("parallel", "antiparallel", "hybrid")))
p_proba <- ggplot(proba_long, aes(prob, fill = class, colour = class)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = class_pal, guide = "none") +
  scale_colour_manual(values = class_pal, name = "Class prob") +
  labs(x = "G4SP class probability", y = "Density",
       title = "G4SP per-class probability distributions") +
  theme_pub()
save_plot(p_proba, "calibration_proba_distributions", "08_topology", fig_root,
          width = 7.5, height = 4.5)

# --- (3) Control benchmark confusion ---------------------------------------
bench <- readr::read_csv("results/tables/g4sp_control_benchmark.csv",
                         show_col_types = FALSE)
if (nrow(bench) > 0 && "predicted_topology" %in% colnames(bench)) {
  lv <- c("parallel", "antiparallel", "hybrid", "mixed")
  conf <- bench %>%
    mutate(reported_topology = factor(reported_topology, levels = lv),
           predicted_topology = factor(predicted_topology, levels = lv)) %>%
    count(reported_topology, predicted_topology, name = "n", .drop = FALSE)
  acc <- mean(bench$correct, na.rm = TRUE)
  p_bench <- ggplot(conf, aes(predicted_topology, reported_topology, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 4) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
    labs(x = "G4SP predicted", y = "Reported (literature)",
         title = "G4SP calibration on known-topology controls",
         subtitle = sprintf("Accuracy = %.0f%% (n = %d)", 100 * acc, nrow(bench))) +
    theme_pub()
  save_plot(p_bench, "calibration_control_benchmark", "08_topology", fig_root,
            width = 6.5, height = 5)
  message("Control benchmark accuracy: ", round(100 * acc, 1), "%")
} else {
  message("No control benchmark rows; skipping confusion plot.")
}

message("Done. Calibration figures written to results/figures/08_topology/.")
