#!/usr/bin/env Rscript
# ============================================================================
# 09_propensity_metrics.R   (env: r_g4)   --- DeepG4 reframe (Q2) ---
#
# DeepG4-vs-peak-strength is uninformative (DeepG4 is a saturating classifier).
# Instead:
#   (a) Cross-validate three INDEPENDENT sequence G4-propensity metrics on our
#       data: DeepG4 probability, pqsfinder score, G4Hunter score. High agreement
#       validates all three; this replaces the strength scatter.
#   (b) Propensity vs helicase-dependence: are peaks LOST upon KO different in
#       intrinsic G4 propensity than stable/gained peaks? (Do the helicases act
#       on the strongest predicted G4s?)
#
# Also writes results/tables/propensity_metrics.csv (peak_id + all metrics +
# topology), reused by rules 11 (feature) and 12 (R-loop).
#
# Inputs:   results/tables/deepg4_scores.csv, peak_pqs.csv, peak_topology.csv,
#           topology_differential.csv, cache/peak_seqs_201bp.fa
# Outputs:  results/tables/propensity_metrics.csv
#           results/tables/propensity_crossvalidation.csv
#           results/tables/propensity_by_dependence.csv
#           results/figures/11_propensity/{propensity_crossvalidation,
#                                          propensity_by_dependence}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/11_propensity")
fig_root <- "results/figures"

# --- Assemble the per-peak propensity table --------------------------------
deepg4 <- readr::read_csv("results/tables/deepg4_scores.csv", show_col_types = FALSE) %>%
  filter(set == "peak") %>% select(peak_id, deepg4_prob)
pqs <- readr::read_csv("results/tables/peak_pqs.csv", show_col_types = FALSE) %>%
  mutate(has_pqs = as.logical(has_pqs)) %>%
  select(peak_id, has_pqs, pqs_seq, pqsfinder_score = score)
topo <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE) %>%
  select(peak_id, topology_final)

# G4Hunter on the 201 bp DeepG4 window (fair comparison) and on the PQS.
win_seqs <- Biostrings::readDNAStringSet("cache/peak_seqs_201bp.fa")
g4h_win <- data.frame(peak_id = names(win_seqs),
                      g4hunter_win = g4hunter_score(win_seqs),
                      stringsAsFactors = FALSE)
pqs$g4hunter_pqs <- ifelse(pqs$has_pqs & !is.na(pqs$pqs_seq),
                           g4hunter_score(pqs$pqs_seq), NA_real_)

prop <- topo %>%
  left_join(deepg4, by = "peak_id") %>%
  left_join(select(pqs, peak_id, pqsfinder_score, g4hunter_pqs), by = "peak_id") %>%
  left_join(g4h_win, by = "peak_id")
readr::write_csv(prop, "results/tables/propensity_metrics.csv")

# --- (a) Cross-validation among the three metrics --------------------------
metric_cols <- c(deepg4_prob = "DeepG4 prob",
                 pqsfinder_score = "pqsfinder score",
                 g4hunter_win = "G4Hunter (201bp)",
                 g4hunter_pqs = "G4Hunter (PQS)")
M <- as.matrix(prop[, names(metric_cols)])
cmat <- suppressWarnings(cor(M, method = "spearman", use = "pairwise.complete.obs"))
cv_long <- as.data.frame(as.table(cmat)) %>%
  setNames(c("metric_a", "metric_b", "spearman_rho"))
readr::write_csv(cv_long, "results/tables/propensity_crossvalidation.csv")

cv_long2 <- cv_long %>%
  mutate(metric_a = factor(metric_a, levels = names(metric_cols), labels = metric_cols),
         metric_b = factor(metric_b, levels = names(metric_cols), labels = metric_cols))
p_heat <- ggplot(cv_long2, aes(metric_a, metric_b, fill = spearman_rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 3.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  labs(x = NULL, y = NULL, title = "Cross-validation of G4-propensity metrics",
       subtitle = "Three independent sequence measures over the same peaks") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Representative scatter: DeepG4 vs G4Hunter (on the PQS — the window-G4Hunter is
# diluted by 201 bp of flank and does not track motif propensity; see heatmap).
sc_df <- prop %>% filter(is.finite(deepg4_prob), is.finite(g4hunter_pqs))
rho_dg <- suppressWarnings(cor(sc_df$deepg4_prob, sc_df$g4hunter_pqs,
                               method = "spearman", use = "complete.obs"))
p_sc <- ggplot(sc_df, aes(g4hunter_pqs, deepg4_prob)) +
  geom_point(alpha = 0.08, size = 0.4, colour = "grey30") +
  # GAM (not loess): at full scale (~130k peaks) loess se=TRUE needs ~25 GB and
  # silently drops the trend line; a penalised spline scales fine and keeps the CI.
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              colour = "#1B9E77", se = TRUE, linewidth = 1) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 4,
           colour = "#1B9E77", label = sprintf("Spearman rho = %.3f", rho_dg)) +
  labs(x = "G4Hunter score (PQS)", y = "DeepG4 probability",
       title = "DeepG4 vs G4Hunter (PQS)") +
  theme_pub()
save_plot(cowplot::plot_grid(p_heat, p_sc, nrow = 1, rel_widths = c(1, 1)),
          "propensity_crossvalidation", "11_propensity", fig_root, width = 12, height = 5)

# --- (b) Propensity vs helicase-dependence ---------------------------------
dz <- readr::read_csv("results/tables/topology_differential.csv", show_col_types = FALSE) %>%
  select(peak_id, ko, dz_class) %>%
  filter(dz_class %in% c("gained", "stable", "lost"))
dep <- dz %>% left_join(prop, by = "peak_id") %>%
  mutate(dz_class = factor(dz_class, levels = c("gained", "stable", "lost")),
         ko = factor(ko, levels = cfg$genotypes))

dep_long <- dep %>%
  pivot_longer(c(deepg4_prob, pqsfinder_score, g4hunter_pqs),
               names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric,
                         levels = c("deepg4_prob", "pqsfinder_score", "g4hunter_pqs"),
                         labels = c("DeepG4 prob", "pqsfinder score", "G4Hunter (PQS)")))

# Kruskal-Wallis: does each metric differ across gained/stable/lost, per KO?
kw <- dep_long %>%
  group_by(ko, metric) %>%
  summarise(kw_stat = tryCatch(kruskal.test(value ~ dz_class)$statistic, error = function(e) NA_real_),
            kw_p    = tryCatch(kruskal.test(value ~ dz_class)$p.value,   error = function(e) NA_real_),
            n = dplyr::n(), .groups = "drop")
readr::write_csv(kw, "results/tables/propensity_by_dependence.csv")

p_dep <- ggplot(dep_long, aes(dz_class, value, fill = dz_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.7) +
  facet_grid(metric ~ ko, scales = "free_y") +
  scale_fill_manual(values = c(gained = "#E41A1C", stable = "grey75", lost = "#377EB8"),
                    guide = "none") +
  labs(x = NULL, y = "Sequence G4-propensity",
       title = "G4 propensity by helicase-dependence (KO vs WT)",
       subtitle = "Do helicases act on the strongest predicted G4s? (Kruskal-Wallis per panel)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_dep, "propensity_by_dependence", "11_propensity", fig_root, width = 11, height = 8)

message("Done. Cross-validation rho (DeepG4 vs G4Hunter win) = ", round(rho_dg, 3))
print(kw)
