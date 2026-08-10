#!/usr/bin/env Rscript
# ============================================================================
# 21b_topology_sensitivity.R   (env: r_g4)   --- NEW IN v1.1 ---
#
# The confidence-stratified sensitivity arm for topology assignment.
#
# Why this exists
# ---------------
# v1.1 stopped back-filling G4ShapePredictor's abstentions with the pqsfinder
# loop-length heuristic (see 21_assign_topology.R and 26d, which measures the
# heuristic at 2/10 on the literature controls against G4SP's 8/10). That is
# the right call, but it has a cost: at the configured g4sp.precision = 0.80
# it moves 29,225 of 130,863 PQS-bearing union peaks into a terminal "mixed"
# stratum, and it shrinks the hybrid stratum from 14,429 to 1,879. The
# "hybrid > parallel > antiparallel" TSS ordering then has to be re-derived on
# an order of magnitude fewer hybrids.
#
# So the size of the hybrid stratum is now entirely a function of one config
# knob. Reporting a result that depends that heavily on a threshold, without
# showing what the threshold does, would repeat the mistake the audit found.
# This script provides the other end of the range: G4SP's plain ARGMAX call,
# where there is no abstention at all (n_mixed = 0) and 11,062 hybrids are
# retained.
#
# Neither arm is "the truth". The primary (0.80) says only what G4SP is
# confident about; the sensitivity arm (argmax) says what G4SP would say if
# forced to choose. A conclusion that holds in both is robust to the knob; one
# that holds in only one is a statement about the knob.
#
# This is CHEAP: 20_g4sp_topology.py already writes the per-class
# probabilities (p_parallel, p_antiparallel, p_hybrid) for every peak, so the
# argmax arm is a row-wise max in R. No model reload, no conda call, no
# G4ShapePredictor re-inference.
#
# Inputs:   results/tables/g4sp_topology.csv   (per-class probabilities)
#           results/tables/peak_topology.csv   (the primary arm, from rule 21)
# Outputs:  results/tables/peak_topology_argmax.csv
#             Same schema as peak_topology.csv, so any downstream script can
#             consume it unmodified (see rule tss_metaprofiles_argmax).
#           results/tables/topology_sensitivity_precision.csv
#             Composition under each arm, plus the primary x argmax cross-tab:
#             where the abstentions actually land when forced to choose.
#           results/figures/08_topology/topology_sensitivity_precision.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/08_topology")
fig_root <- "results/figures"

primary_precision <- cfg$g4sp$precision %||% "default"
sens_precision    <- cfg$g4sp$sensitivity_precision %||% "default"

# Only the argmax arm is computable here from probabilities alone. A numeric
# sensitivity_precision would need the per-class threshold .pkl, which lives on
# the python side (26a_threshold_sweep.py already sweeps it) -- so refuse
# rather than silently produce the wrong thing.
if (!(is.character(sens_precision) &&
      tolower(trimws(sens_precision)) %in% c("default", "argmax", "none"))) {
  stop("g4sp.sensitivity_precision must be 'default' (plain argmax). ",
       "Got: ", sens_precision, ". For numeric precisions the per-class ",
       "thresholds live in the G4SP .pkl -- use 26a_threshold_sweep.py, which ",
       "already sweeps calibration.precision_grid.")
}

g4sp <- readr::read_csv("results/tables/g4sp_topology.csv", show_col_types = FALSE)
prim <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)

prob_cols <- c(parallel = "p_parallel", antiparallel = "p_antiparallel",
               hybrid = "p_hybrid")
missing <- setdiff(unname(prob_cols), names(g4sp))
if (length(missing) > 0)
  stop("g4sp_topology.csv is missing the per-class probability columns (",
       paste(missing, collapse = ", "), "). 20_g4sp_topology.py writes these; ",
       "re-run rule g4sp_topology.")

# --- Plain argmax over the three class probabilities ------------------------
# This reproduces G4SP's predict_G4 rule with thresholds disabled, which is
# exactly what precision: default does in 20_g4sp_topology.py.
P <- as.matrix(g4sp[, unname(prob_cols)])
has_prob <- rowSums(is.na(P)) == 0
argmax_class <- rep(NA_character_, nrow(P))
argmax_class[has_prob] <- names(prob_cols)[max.col(P[has_prob, , drop = FALSE],
                                                   ties.method = "first")]
argmax_conf <- rep(NA_real_, nrow(P))
argmax_conf[has_prob] <- apply(P[has_prob, , drop = FALSE], 1, max)

# A peak with no PQS has no sequence to classify and stays no_canonical_PQS in
# both arms -- the argmax arm changes G4SP's abstentions, not its inputs.
no_pqs <- prim$topology_final[match(g4sp$peak_id, prim$peak_id)] == "no_canonical_PQS"
no_pqs[is.na(no_pqs)] <- TRUE

sens <- tibble::tibble(
  peak_id         = g4sp$peak_id,
  topology_final  = ifelse(no_pqs | is.na(argmax_class),
                           "no_canonical_PQS", argmax_class),
  topology_g4sp   = argmax_class,
  # Carried through from the primary table so the two arms stay joinable on a
  # common schema; the heuristic is never used to assign a label in either.
  topology_heuristic = prim$topology_heuristic[match(g4sp$peak_id, prim$peak_id)],
  topology_source = ifelse(no_pqs | is.na(argmax_class), "none", "g4sp_argmax"),
  g4sp_confidence = argmax_conf)
sens$agree <- !is.na(sens$topology_heuristic) &
  sens$topology_final %in% topology_definite &
  sens$topology_final == sens$topology_heuristic

sens$topology_final <- factor(sens$topology_final, levels = topology_levels)
readr::write_csv(sens, "results/tables/peak_topology_argmax.csv")

# --- Composition under each arm + the cross-tab -----------------------------
comp <- bind_rows(
  prim %>% count(topology_final, name = "n") %>%
    mutate(arm = "primary", precision = as.character(primary_precision)),
  sens %>% count(topology_final, name = "n") %>%
    mutate(arm = "sensitivity", precision = "default (argmax)")) %>%
  group_by(arm) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
  select(arm, precision, topology_final, n, frac)

# Where do the primary arm's abstentions go when G4SP is forced to choose?
# This is the number that says whether "mixed" is a coherent stratum or just
# the low-confidence tail of the three real classes.
xt <- prim %>%
  select(peak_id, primary = topology_final) %>%
  inner_join(sens %>% select(peak_id, argmax = topology_final), by = "peak_id") %>%
  count(primary, argmax, name = "n")

readr::write_csv(comp, "results/tables/topology_sensitivity_precision.csv")
readr::write_csv(xt, "results/tables/topology_sensitivity_crosstab.csv")

mixed_fate <- xt %>% filter(primary == "mixed") %>%
  mutate(frac = n / sum(n)) %>% arrange(desc(n))
if (nrow(mixed_fate) > 0) {
  message("Fate of the primary arm's 'mixed' peaks under plain argmax:")
  print(as.data.frame(mixed_fate))
}

# --- Figure: the two arms side by side --------------------------------------
p_comp <- ggplot(comp, aes(topology_final, n, fill = topology_final)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  facet_wrap(~ paste0(arm, ": precision ", precision), nrow = 1) +
  scale_fill_manual(values = topology_palette, guide = "none") +
  scale_x_discrete(labels = setNames(topology_labels, topology_levels)) +
  labs(x = NULL, y = "Union peak count",
       title = "Topology composition is a function of the G4SP precision threshold",
       subtitle = paste0(
         "Primary keeps G4SP's abstentions as 'mixed'; the sensitivity arm ",
         "forces a call (argmax, no abstention).\n",
         "A conclusion that holds in both arms is robust to this knob.")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_comp, "topology_sensitivity_precision", "08_topology", fig_root,
          width = 11, height = 5.5)

message("Done. Composition by arm:")
print(as.data.frame(comp))
