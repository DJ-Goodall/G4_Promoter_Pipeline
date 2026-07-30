#!/usr/bin/env Rscript
# ============================================================================
# 11_feature_partition.R   (env: r_g4)   --- Q1c ---
#
# Where do G4 topologies and high-propensity G4s sit in the genome?
#   (a) Topology x genomic feature (promoter / 5'UTR / intron1 / gene-body /
#       enhancer / intergenic): composition + chi-square + standardized residuals.
#   (b) Sequence propensity (DeepG4 / pqsfinder / G4Hunter) by feature.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           results/tables/propensity_metrics.csv,
#           cache_v2 feature region RDS (config: feature.region_files)
# Outputs:  results/tables/topology_by_feature.csv,
#           results/tables/topology_feature_tests.csv,
#           results/tables/propensity_by_feature.csv,
#           results/figures/12_feature/{topology_by_feature,propensity_by_feature}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/12_feature")
fig_root <- "results/figures"

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
prop  <- readr::read_csv("results/tables/propensity_metrics.csv", show_col_types = FALSE)
union$topology <- factor(topo$topology_final[match(names(union), topo$peak_id)],
                         levels = topology_levels)

# --- Feature annotation (priority order from config) ------------------------
region_files <- cfg$feature$region_files
region_list <- lapply(names(region_files), function(nm) {
  gr <- readRDS(file.path(cfg$paths$cache_v2, region_files[[nm]]))
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr
})
names(region_list) <- names(region_files)
union <- annotate_peak_regions(union, region_list, intergenic_label = "intergenic")

df <- data.frame(peak_id = names(union),
                 region = as.character(S4Vectors::mcols(union)$region),
                 topology = as.character(union$topology),
                 stringsAsFactors = FALSE) %>%
  left_join(prop %>% select(peak_id, deepg4_prob, pqsfinder_score, g4hunter_pqs),
            by = "peak_id")
feature_levels <- c(names(region_files), "intergenic")
df$region   <- factor(df$region, levels = feature_levels)
df$topology <- factor(df$topology, levels = topology_levels)

# --- (a) Topology x feature -------------------------------------------------
comp <- df %>% count(region, topology, name = "n") %>%
  group_by(region) %>% mutate(frac = n / sum(n)) %>% ungroup()
readr::write_csv(comp, "results/tables/topology_by_feature.csv")

definite <- c("parallel", "antiparallel", "hybrid")
ct <- table(factor(df$topology[df$topology %in% definite], levels = definite),
            df$region[df$topology %in% definite])
chisq <- suppressWarnings(chisq.test(ct))
stdres <- as.data.frame(as.table(chisq$stdres)) %>%
  setNames(c("topology", "region", "std_residual"))
stdres$chisq_p <- chisq$p.value
readr::write_csv(stdres, "results/tables/topology_feature_tests.csv")

p_comp <- ggplot(comp, aes(region, frac, fill = topology)) +
  geom_col(colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = topology_palette,
                    labels = setNames(topology_labels, topology_levels), name = "Topology") +
  labs(x = NULL, y = "Fraction of peaks",
       title = "G4 topology composition by genomic feature") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_comp, "topology_by_feature", "12_feature", fig_root, width = 9, height = 5)

# --- (b) Propensity by feature ----------------------------------------------
prop_long <- df %>%
  pivot_longer(c(deepg4_prob, pqsfinder_score, g4hunter_pqs),
               names_to = "metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(metric = factor(metric,
                         levels = c("deepg4_prob", "pqsfinder_score", "g4hunter_pqs"),
                         labels = c("DeepG4 prob", "pqsfinder score", "G4Hunter (PQS)")))
prop_summ <- prop_long %>% group_by(region, metric) %>%
  summarise(median = median(value, na.rm = TRUE), mean = mean(value, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")
readr::write_csv(prop_summ, "results/tables/propensity_by_feature.csv")

p_prop <- ggplot(prop_long, aes(region, value, fill = region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "Sequence G4-propensity",
       title = "G4 propensity by genomic feature") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1),
                      legend.position = "none")
save_plot(p_prop, "propensity_by_feature", "12_feature", fig_root, width = 8, height = 9)

message("Done. Topology x feature (chi-square p = ", signif(chisq$p.value, 3), ")")
print(comp)
