#!/usr/bin/env Rscript
# ============================================================================
# 12_rloop_coupling.R   (env: r_g4)   --- Q1b ---
#
# Do G4s (by topology / DeepG4 propensity) couple to R-loops, and does helicase
# loss change R-loops specifically at certain G4 topologies?
#   (a) Fraction of G4 peaks overlapping an S9.6 R-loop peak, by topology and by
#       DeepG4-propensity quartile.
#   (b) R-loop CUT&Tag signal at G4 peaks (WT) by topology.
#   (c) Differential R-loop signal at G4 peaks KO-vs-WT (voom), by topology.
#
# Inputs:   cache/peaks.rds, peak_topology.csv, propensity_metrics.csv,
#           cache_v2 R-loop union peaks (config rloop.union_rds), R-loop bigWigs
# Outputs:  results/tables/g4_rloop_overlap.csv,
#           results/tables/rloop_signal_by_topology.csv,
#           results/tables/rloop_differential_at_g4.csv,
#           results/figures/13_rloop/{rloop_fraction_by_topology,
#           rloop_signal_by_topology,rloop_differential_at_g4}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/13_rloop")
fig_root <- "results/figures"

genotypes <- cfg$genotypes
ref_geno  <- cfg$ref_genotype
kos       <- setdiff(genotypes, ref_geno)
definite  <- c("parallel", "antiparallel", "hybrid")
topo_lab  <- setNames(topology_labels, topology_levels)

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
prop  <- readr::read_csv("results/tables/propensity_metrics.csv", show_col_types = FALSE)
union$topology <- factor(topo$topology_final[match(names(union), topo$peak_id)],
                         levels = topology_levels)

# --- (a) Overlap with R-loop peaks -----------------------------------------
rloop <- readRDS(file.path(cfg$paths$cache_v2, cfg$rloop$union_rds))
GenomeInfoDb::seqlevelsStyle(rloop) <- "UCSC"
union$rloop_pos <- IRanges::overlapsAny(union, rloop, ignore.strand = TRUE)

ov_df <- data.frame(peak_id = names(union),
                    topology = as.character(union$topology),
                    rloop_pos = union$rloop_pos, stringsAsFactors = FALSE) %>%
  left_join(prop %>% select(peak_id, deepg4_prob), by = "peak_id")

# by topology
by_topo <- ov_df %>% filter(topology %in% c(definite, "no_canonical_PQS")) %>%
  group_by(group_type = "topology", group = topology) %>%
  summarise(n = dplyr::n(), n_rloop_pos = sum(rloop_pos),
            frac_rloop_pos = mean(rloop_pos), .groups = "drop")
# by DeepG4 propensity quartile
nbin <- cfg$rloop$propensity_bins %||% 4
ov_df$prop_bin <- tryCatch(
  cut(ov_df$deepg4_prob, breaks = quantile(ov_df$deepg4_prob, probs = seq(0, 1, length.out = nbin + 1),
      na.rm = TRUE), include.lowest = TRUE, labels = paste0("Q", seq_len(nbin))),
  error = function(e) factor(rep(NA, nrow(ov_df))))
by_prop <- ov_df %>% filter(!is.na(prop_bin)) %>%
  group_by(group_type = "deepg4_quartile", group = as.character(prop_bin)) %>%
  summarise(n = dplyr::n(), n_rloop_pos = sum(rloop_pos),
            frac_rloop_pos = mean(rloop_pos), .groups = "drop")
overlap_summary <- dplyr::bind_rows(by_topo, by_prop)
readr::write_csv(overlap_summary, "results/tables/g4_rloop_overlap.csv")

# chi-square: R-loop+ vs topology (definite)
ct <- table(ov_df$topology[ov_df$topology %in% definite],
            ov_df$rloop_pos[ov_df$topology %in% definite])
chisq_p <- suppressWarnings(tryCatch(chisq.test(ct)$p.value, error = function(e) NA_real_))

p_frac <- ggplot(by_topo %>% mutate(group = factor(group, levels = topology_levels)),
                 aes(group, frac_rloop_pos, fill = group)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac_rloop_pos)), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = topology_palette, guide = "none") +
  scale_x_discrete(labels = topo_lab) +
  labs(x = NULL, y = "Fraction overlapping an R-loop peak",
       title = "G4-R-loop co-occupancy by topology",
       subtitle = sprintf("chi-square p = %s", signif(chisq_p, 3))) +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_frac, "rloop_fraction_by_topology", "13_rloop", fig_root, width = 7.5, height = 5)

# --- (b) R-loop signal at G4 peaks (WT) by topology ------------------------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$rloop$assay)
wt_files <- bw_meta$filepath[bw_meta$genotype == ref_geno]
rloop_sig <- mean_replicate_signal(union, wt_files)
evict_bigwigs(wt_files)
sig_df <- data.frame(peak_id = names(union),
                     topology = factor(as.character(union$topology), levels = topology_levels),
                     rloop_signal_WT = rloop_sig, stringsAsFactors = FALSE)
sig_summ <- sig_df %>% filter(topology %in% c(definite, "no_canonical_PQS")) %>%
  group_by(topology) %>%
  summarise(median = median(rloop_signal_WT, na.rm = TRUE),
            mean = mean(rloop_signal_WT, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
readr::write_csv(sig_summ, "results/tables/rloop_signal_by_topology.csv")

p_sig <- ggplot(sig_df %>% filter(topology %in% c(definite, "no_canonical_PQS"),
                                  is.finite(rloop_signal_WT)),
                aes(topology, log1p(rloop_signal_WT), fill = topology)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = topology_palette, guide = "none") +
  scale_x_discrete(labels = topo_lab) +
  labs(x = NULL, y = "WT R-loop signal at G4 peak  [log1p]",
       title = "R-loop occupancy at G4 peaks by topology (WT)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_sig, "rloop_signal_by_topology", "13_rloop", fig_root, width = 7.5, height = 5)

# --- (c) Differential R-loop signal at G4 peaks (KO vs WT) by topology ------
dz_rloop <- dplyr::bind_rows(lapply(kos, function(k) {
  sub <- bw_meta[bw_meta$genotype %in% c(ref_geno, k), , drop = FALSE]
  message("Differential R-loop at G4 peaks: ", k, " vs ", ref_geno, " (", nrow(sub), " bigWigs)")
  dz  <- peak_voom_lfc(union, sub, ref_genotype = ref_geno)
  cls <- classify_peak_dz(dz, padj_thresh = cfg$differential$voom_padj,
                          lfc_thresh = cfg$differential$voom_lfc)
  evict_bigwigs(sub$filepath)
  data.frame(peak_id = names(union), ko = k, dz_class = as.character(cls),
             topology = as.character(union$topology), stringsAsFactors = FALSE)
}))
dz_summary <- dz_rloop %>%
  filter(topology %in% definite, dz_class %in% c("gained", "stable", "lost")) %>%
  mutate(topology = factor(topology, levels = definite),
         dz_class = factor(dz_class, levels = c("gained", "stable", "lost")),
         ko = factor(ko, levels = kos)) %>%
  count(ko, topology, dz_class, name = "n") %>%
  group_by(ko, topology) %>% mutate(frac = n / sum(n)) %>% ungroup()
readr::write_csv(dz_summary, "results/tables/rloop_differential_at_g4.csv")

p_dz <- ggplot(dz_summary, aes(topology, frac, fill = dz_class)) +
  geom_col(colour = "white", linewidth = 0.2) +
  facet_wrap(~ ko) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = topo_lab) +
  scale_fill_manual(values = c(gained = "#E41A1C", stable = "grey75", lost = "#377EB8"),
                    name = "R-loop vs WT") +
  labs(x = NULL, y = "Fraction of G4 peaks",
       title = "Differential R-loop signal at G4 peaks by topology (KO vs WT)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_dz, "rloop_differential_at_g4", "13_rloop", fig_root, width = 10, height = 5)

message("Done. R-loop+ fraction by topology:")
print(by_topo)
