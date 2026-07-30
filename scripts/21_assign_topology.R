#!/usr/bin/env Rscript
# ============================================================================
# 05_assign_topology.R   (env: r_g4)
#
# Combine G4ShapePredictor topology calls with an independent pqsfinder
# loop-length heuristic into one per-peak topology label, and report concordance.
#
# Inputs:   results/tables/peak_pqs.csv       (loop lengths -> heuristic)
#           results/tables/g4sp_topology.csv  (G4SP class)
# Outputs:  results/tables/peak_topology.csv
#           results/tables/topology_concordance.csv
#           results/figures/08_topology/{topology_counts,g4sp_vs_heuristic_concordance}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/08_topology")
fig_root <- "results/figures"

par_max  <- cfg$heuristic$parallel_max_loop
anti_min <- cfg$heuristic$antiparallel_min_long_loops

pqs  <- readr::read_csv("results/tables/peak_pqs.csv", show_col_types = FALSE) %>%
  mutate(has_pqs = as.logical(has_pqs))
g4sp <- readr::read_csv("results/tables/g4sp_topology.csv", show_col_types = FALSE)

# --- Heuristic topology from loop lengths ----------------------------------
heuristic_topology <- function(has_pqs, ll1, ll2, ll3) {
  if (isFALSE(has_pqs) || is.na(has_pqs)) return(NA_character_)
  loops <- c(ll1, ll2, ll3)
  if (anyNA(loops)) return(NA_character_)
  long <- sum(loops > par_max)
  if (long == 0) "parallel" else if (long >= anti_min) "antiparallel" else "hybrid"
}

pqs <- pqs %>%
  rowwise() %>%
  mutate(topology_heuristic = heuristic_topology(has_pqs, ll1, ll2, ll3)) %>%
  ungroup()

# --- Merge + final label ----------------------------------------------------
definite <- c("parallel", "antiparallel", "hybrid")
merged <- pqs %>%
  select(peak_id, has_pqs, topology_heuristic) %>%
  left_join(select(g4sp, peak_id, topology_g4sp = g4sp_class, g4sp_confidence),
            by = "peak_id") %>%
  mutate(
    g4sp_definite = topology_g4sp %in% definite,
    no_pqs = is.na(has_pqs) | !has_pqs,
    topology_final = case_when(
      no_pqs                          ~ "no_canonical_PQS",
      g4sp_definite                   ~ topology_g4sp,
      !is.na(topology_heuristic)      ~ topology_heuristic,
      TRUE                            ~ "no_canonical_PQS"
    ),
    topology_source = case_when(
      no_pqs                     ~ "none",
      g4sp_definite              ~ "g4sp",
      !is.na(topology_heuristic) ~ "heuristic",
      TRUE                       ~ "none"
    ),
    agree = g4sp_definite & !is.na(topology_heuristic) &
            (topology_g4sp == topology_heuristic)
  )

merged$topology_final <- factor(merged$topology_final, levels = topology_levels)
readr::write_csv(
  merged %>% select(peak_id, topology_final, topology_g4sp, topology_heuristic,
                    topology_source, g4sp_confidence, agree),
  "results/tables/peak_topology.csv")

# --- Concordance (only where both definite) --------------------------------
conc <- merged %>% filter(g4sp_definite, !is.na(topology_heuristic))
if (nrow(conc) > 0) {
  conf <- conc %>% count(topology_g4sp, topology_heuristic, name = "n")
  readr::write_csv(conf, "results/tables/topology_concordance.csv")

  # Cohen's kappa (manual; avoids extra dependency)
  lv <- definite
  tab <- table(factor(conc$topology_g4sp, lv), factor(conc$topology_heuristic, lv))
  po <- sum(diag(tab)) / sum(tab)
  pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  kappa <- (po - pe) / (1 - pe)
  message(sprintf("G4SP vs heuristic: agreement %.1f%% (Cohen's kappa = %.3f, n = %d)",
                  100 * po, kappa, nrow(conc)))

  p_conc <- ggplot(conf, aes(topology_heuristic, topology_g4sp, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 4) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
    labs(x = "pqsfinder loop heuristic", y = "G4ShapePredictor",
         title = "Topology call concordance",
         subtitle = sprintf("Cohen's kappa = %.3f (n = %d definite both)", kappa, nrow(conc))) +
    theme_pub()
  save_plot(p_conc, "g4sp_vs_heuristic_concordance", "08_topology", fig_root,
            width = 6.5, height = 5.5)
} else {
  message("No peaks with definite calls from BOTH methods; skipping concordance ",
          "(G4SP likely unavailable -> final labels are heuristic-based).")
  readr::write_csv(data.frame(), "results/tables/topology_concordance.csv")
}

# --- Overall topology composition ------------------------------------------
counts <- merged %>% count(topology_final, name = "n")
p_counts <- ggplot(counts, aes(topology_final, n, fill = topology_final)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = topology_palette, guide = "none") +
  scale_x_discrete(labels = setNames(topology_labels, topology_levels)) +
  labs(x = NULL, y = "Union peak count",
       title = "G4 topology composition (all union peaks)",
       subtitle = "Final label = G4SP where definite, else pqsfinder loop heuristic") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_counts, "topology_counts", "08_topology", fig_root, width = 7, height = 5)

message("Done. Final topology breakdown:")
print(counts)
