#!/usr/bin/env Rscript
# ============================================================================
# 05_assign_topology.R   (env: r_g4)
#
# Assign one per-peak topology label from G4ShapePredictor, keeping G4SP's
# abstentions ("mixed") as their own stratum, and report concordance with an
# independent pqsfinder loop-length heuristic.
#
# ---------------------------------------------------------------------------
# v1.1 CHANGE (2026-08-07 audit) -- this is the change that moves a conclusion.
#
# v1.0 back-filled every peak G4SP declined to classify with the loop-length
# heuristic:
#     topology_final = case_when(
#       no_pqs                     ~ "no_canonical_PQS",
#       g4sp_definite              ~ topology_g4sp,
#       !is.na(topology_heuristic) ~ topology_heuristic,   # <- removed in v1.1
#       TRUE                       ~ "no_canonical_PQS")
#
# At the configured g4sp.precision = 0.80 that fallback fired for 29,225 of
# 130,863 union peaks with a PQS (22.3%), and because G4SP abstains
# preferentially on the non-parallel classes it supplied the MAJORITY of two
# of the three strata (from topology_composition_compare.csv vs the 0.80 row
# of topology_precision_sweep.csv):
#
#     stratum        final n   G4SP n   heuristic n   % heuristic
#     parallel       104,336   95,660         8,676         8.3%
#     antiparallel    12,098    4,099         7,999        66.1%
#     hybrid          14,429    1,879        12,550        87.0%
#
# The two methods do not measure the same thing. Of the heuristic's 41,122
# hybrid calls only 741 (1.8%) are G4SP-hybrid; of its 25,539 antiparallel
# calls only 191 (0.7%) are G4SP-antiparallel (topology_concordance.csv).
# Cohen's kappa is -0.0098, though note kappa is depressed here by the very
# different marginals (G4SP 94.1% parallel vs heuristic 34.4%) -- the
# cross-classification rates above are the more honest summary and are what
# v1.1 reports.
#
# So the v1.0 "hybrid > parallel > antiparallel" TSS ordering was largely a
# statement about "PQS with exactly one loop longer than 3 nt". v1.1:
#   - keeps "mixed" as a terminal stratum (no back-fill);
#   - persists kappa AND the cross-classification rates to a table, rather
#     than a message() into a git-ignored log and a figure subtitle;
#   - aggregates the topology_source column, which v1.0 recorded but never
#     surfaced anywhere.
# See docs/audit-2026-08/REVIEW.md section 2.1 and scripts/26b_control_benchmark.py, which
# now scores the heuristic on the same literature controls as G4SP.
# ---------------------------------------------------------------------------
#
# Inputs:   results/tables/peak_pqs.csv       (loop lengths -> heuristic)
#           results/tables/g4sp_topology.csv  (G4SP class)
# Outputs:  results/tables/peak_topology.csv
#           results/tables/topology_concordance.csv
#           results/tables/topology_concordance_stats.csv   [v1.1]
#           results/tables/topology_source_composition.csv  [v1.1]
#           results/figures/08_topology/{topology_counts,topology_source_composition,
#             g4sp_vs_heuristic_concordance}.{pdf,png}
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
    # v1.1: NO heuristic back-fill. A peak with a PQS that G4SP declined to
    # classify is labelled "mixed" and stays there. The heuristic call is still
    # computed and retained in its own column for the concordance analysis, but
    # it never becomes the reported topology.
    topology_final = case_when(
      no_pqs        ~ "no_canonical_PQS",
      g4sp_definite ~ topology_g4sp,
      TRUE          ~ "mixed"
    ),
    topology_source = case_when(
      no_pqs        ~ "none",
      g4sp_definite ~ "g4sp",
      TRUE          ~ "g4sp_abstained"
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

  # --- v1.1: PERSIST the concordance statistics ----------------------------
  # In v1.0 kappa existed only in a message() to a git-ignored logs/ dir and in
  # a figure subtitle under a git-ignored results/figures/ -- i.e. in no
  # committed artefact at all. It is now a table.
  #
  # Note on interpretation: kappa is depressed here by the very different
  # marginals (G4SP calls 94% parallel, the heuristic 34%), and with marginals
  # that divergent kappa is a known-unreliable summary -- quoting "worse than
  # chance" invites a methodological rebuttal. The marginal-free
  # cross-classification rates below say the same thing unarguably: when the
  # heuristic says hybrid, G4SP almost never agrees.
  ppv <- function(cls) {
    denom <- sum(tab[, cls])
    if (denom == 0) NA_real_ else tab[cls, cls] / denom
  }
  conc_stats <- data.frame(
    metric = c("n_both_definite", "observed_agreement", "expected_agreement",
               "cohens_kappa",
               "frac_heuristic_parallel_confirmed_by_g4sp",
               "frac_heuristic_antiparallel_confirmed_by_g4sp",
               "frac_heuristic_hybrid_confirmed_by_g4sp",
               "g4sp_marginal_parallel", "heuristic_marginal_parallel"),
    value = c(nrow(conc), po, pe, kappa,
              ppv("parallel"), ppv("antiparallel"), ppv("hybrid"),
              sum(tab["parallel", ]) / sum(tab), sum(tab[, "parallel"]) / sum(tab)),
    note = c("peaks where both methods gave a definite call",
             "diagonal / total",
             "sum(row_i * col_i) / total^2",
             "(po - pe) / (1 - pe); depressed by marginal mismatch, see header",
             "of all heuristic 'parallel' calls, fraction G4SP also calls parallel",
             "of all heuristic 'antiparallel' calls, fraction G4SP agrees",
             "of all heuristic 'hybrid' calls, fraction G4SP agrees",
             "G4SP class balance", "heuristic class balance"),
    stringsAsFactors = FALSE)
  readr::write_csv(conc_stats, "results/tables/topology_concordance_stats.csv")
  message(sprintf(
    "  heuristic calls confirmed by G4SP: parallel %.1f%%, antiparallel %.1f%%, hybrid %.1f%%",
    100 * ppv("parallel"), 100 * ppv("antiparallel"), 100 * ppv("hybrid")))

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
       subtitle = paste0("Final label = G4ShapePredictor at precision ",
                         cfg$g4sp$precision,
                         "; 'mixed' = G4SP declined to call (no back-fill)")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_counts, "topology_counts", "08_topology", fig_root, width = 7, height = 5)

# --- v1.1: surface topology_source ------------------------------------------
# The column existed in v1.0 (it was even written to peak_topology.csv) but was
# never counted, plotted or reported anywhere, so a reader had no way to see
# how much of each stratum was a model call. It is now a table and a figure.
src_comp <- merged %>%
  count(topology_final, topology_source, name = "n") %>%
  group_by(topology_final) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
readr::write_csv(src_comp, "results/tables/topology_source_composition.csv")

p_src <- ggplot(src_comp, aes(topology_final, n, fill = topology_source)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  scale_x_discrete(labels = setNames(topology_labels, topology_levels)) +
  scale_fill_manual(values = c(g4sp = "#1B9E77", g4sp_abstained = "#D95F02",
                               none = "grey70")) +
  labs(x = NULL, y = "Union peak count", fill = "call source",
       title = "Provenance of every topology call",
       subtitle = "v1.1 makes no heuristic substitutions; abstentions stay 'mixed'") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_src, "topology_source_composition", "08_topology", fig_root,
          width = 7.5, height = 5)

message("Done. Final topology breakdown:")
print(counts)
message("Call provenance:")
print(merged %>% count(topology_source, name = "n"))
