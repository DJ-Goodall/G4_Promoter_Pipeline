#!/usr/bin/env Rscript
# ============================================================================
# 07_gain_loss.R   (env: r_g4)   --- Goal G3 ---
#
# Does helicase KO gain or lose particular G4 topologies?
#   (a) Composition: topology breakdown of per-genotype peak sets (counts/fractions)
#       + chi-square (topology x genotype) and per-topology KO-vs-WT proportion tests.
#   (b) Differential signal: limma-voom LFC of per-peak G4 signal KO-vs-WT on the
#       union peak set, classified gained/stable/lost, summarised per topology.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv, G4 bigWigs
# Outputs:  results/tables/topology_composition{,_tests}.csv
#           results/tables/topology_differential{,_summary}.csv
#           results/figures/09_gain_loss/*.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(S4Vectors); library(dplyr); library(tidyr)
  library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/09_gain_loss")
fig_root <- "results/figures"

genotypes <- cfg$genotypes
ref_geno  <- cfg$ref_genotype
kos       <- setdiff(genotypes, ref_geno)
definite  <- c("parallel", "antiparallel", "hybrid")
topo_lab  <- setNames(topology_labels, topology_levels)

peaks_obj    <- readRDS("cache/peaks.rds")
union        <- peaks_obj$union
per_genotype <- peaks_obj$per_genotype
topo         <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- factor(topo$topology_final[match(names(union), topo$peak_id)],
                         levels = topology_levels)

# ===========================================================================
# (a) Per-genotype topology composition
# ===========================================================================
# Each per-genotype peak inherits the topology of the union peak it falls in.
comp <- dplyr::bind_rows(lapply(genotypes, function(g) {
  pg <- per_genotype[[g]]
  if (length(pg) == 0) return(NULL)
  hit <- GenomicRanges::findOverlaps(pg, union, select = "first", ignore.strand = TRUE)
  data.frame(genotype = g,
             topology = as.character(union$topology[hit]),
             stringsAsFactors = FALSE)
}))
comp$topology <- factor(comp$topology, levels = topology_levels)
comp$genotype <- factor(comp$genotype, levels = genotypes)

comp_counts <- comp %>% count(genotype, topology, name = "n") %>%
  group_by(genotype) %>% mutate(frac = n / sum(n)) %>% ungroup()
readr::write_csv(comp_counts, "results/tables/topology_composition.csv")

p_stack <- ggplot(comp_counts, aes(genotype, n, fill = topology)) +
  geom_col(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = topology_palette, labels = topo_lab, name = "Topology") +
  labs(x = NULL, y = "Peak count", title = "G4 topology composition per genotype") +
  theme_pub()
save_plot(p_stack, "composition_stacked", "09_gain_loss", fig_root, width = 8, height = 5)

p_frac <- ggplot(comp_counts, aes(genotype, frac, fill = topology)) +
  geom_col(colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = topology_palette, labels = topo_lab, name = "Topology") +
  labs(x = NULL, y = "Fraction of peaks",
       title = "G4 topology fractions per genotype") +
  theme_pub()
save_plot(p_frac, "composition_fraction", "09_gain_loss", fig_root, width = 8, height = 5)

# Statistical tests on composition (definite topologies only)
comp_def <- comp %>% filter(topology %in% definite) %>%
  mutate(topology = factor(topology, levels = definite))
ct <- table(comp_def$topology, comp_def$genotype)
chisq <- suppressWarnings(chisq.test(ct))
test_rows <- list(data.frame(test = "chisq_topology_x_genotype",
                             topology = NA, comparison = NA,
                             statistic = unname(chisq$statistic),
                             p_value = chisq$p.value))
# Per-topology KO-vs-WT proportion test (is topology t enriched/depleted in KO?)
for (k in kos) for (t in definite) {
  sub <- comp_def %>% filter(genotype %in% c(ref_geno, k))
  tab <- table(factor(sub$topology == t, c(TRUE, FALSE)),
               factor(sub$genotype, c(k, ref_geno)))
  pt <- suppressWarnings(tryCatch(prop.test(tab[1, ], colSums(tab))$p.value,
                                  error = function(e) NA_real_))
  test_rows[[length(test_rows) + 1]] <- data.frame(
    test = "prop_test_topology_KO_vs_WT", topology = t,
    comparison = paste0(k, "_vs_", ref_geno),
    statistic = NA_real_, p_value = pt)
}
readr::write_csv(dplyr::bind_rows(test_rows), "results/tables/topology_composition_tests.csv")

# ===========================================================================
# (b) Differential peak signal (gained / stable / lost) by topology
# ===========================================================================
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)
chrom_sizes <- default_chrom_sizes(cfg$std_chroms)

dz_all <- dplyr::bind_rows(lapply(kos, function(k) {
  sub <- bw_meta[bw_meta$genotype %in% c(ref_geno, k), , drop = FALSE]
  message("Differential G4 signal: ", k, " vs ", ref_geno,
          " (", nrow(sub), " bigWigs)")
  dz <- peak_voom_lfc(union, sub, ref_genotype = ref_geno)
  cls <- classify_peak_dz(dz, padj_thresh = cfg$differential$voom_padj,
                          lfc_thresh = cfg$differential$voom_lfc)
  evict_bigwigs(sub$filepath)
  data.frame(peak_id = names(union), ko = k,
             logFC = dz$logFC, adj.P.Val = dz$adj.P.Val,
             dz_class = as.character(cls),
             topology = as.character(union$topology),
             stringsAsFactors = FALSE)
}))
dz_all$topology <- factor(dz_all$topology, levels = topology_levels)
dz_all$dz_class <- factor(dz_all$dz_class, levels = c("gained", "stable", "lost", "ambiguous"))
dz_all$ko <- factor(dz_all$ko, levels = kos)
readr::write_csv(dz_all, "results/tables/topology_differential.csv")

# Summary: gained / lost fraction within each topology x KO
dz_summary <- dz_all %>%
  filter(topology %in% definite, dz_class %in% c("gained", "stable", "lost")) %>%
  mutate(topology = factor(topology, levels = definite),
         dz_class = factor(dz_class, levels = c("gained", "stable", "lost"))) %>%
  count(ko, topology, dz_class, name = "n") %>%
  group_by(ko, topology) %>% mutate(frac = n / sum(n)) %>% ungroup()
readr::write_csv(dz_summary, "results/tables/topology_differential_summary.csv")

# Plot: gained/stable/lost fraction per topology, faceted by KO
p_dz <- ggplot(dz_summary, aes(topology, frac, fill = dz_class)) +
  geom_col(colour = "white", linewidth = 0.2) +
  facet_wrap(~ ko) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = topo_lab) +
  scale_fill_manual(values = c(gained = "#E41A1C", stable = "grey75", lost = "#377EB8"),
                    name = "vs WT") +
  labs(x = NULL, y = "Fraction of peaks",
       title = "Differential G4 signal by topology (KO vs WT)",
       subtitle = sprintf("limma-voom on union peaks; |log2FC|>%.2f & adj.P<%.2g",
                          cfg$differential$voom_lfc, cfg$differential$voom_padj)) +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_dz, "differential_class_fraction", "09_gain_loss", fig_root,
          width = 10, height = 5)

# Plot: absolute gained vs lost counts per topology x KO
gl <- dz_summary %>% filter(dz_class %in% c("gained", "lost"))
p_gl <- ggplot(gl, aes(topology, n, fill = dz_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ ko, scales = "free_y") +
  scale_x_discrete(labels = topo_lab) +
  scale_fill_manual(values = c(gained = "#E41A1C", lost = "#377EB8"), name = "vs WT") +
  labs(x = NULL, y = "Peak count",
       title = "Gained vs lost G4 peaks by topology (KO vs WT)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_gl, "differential_gainloss_counts", "09_gain_loss", fig_root,
          width = 10, height = 5)

# Per-KO chi-square: is dz_class associated with topology?
chisq_rows <- dplyr::bind_rows(lapply(kos, function(k) {
  d <- dz_all %>% filter(ko == k, topology %in% definite,
                         dz_class %in% c("gained", "stable", "lost"))
  tb <- table(droplevels(d$topology), droplevels(d$dz_class))
  ch <- suppressWarnings(tryCatch(chisq.test(tb), error = function(e) NULL))
  data.frame(ko = k,
             statistic = if (is.null(ch)) NA_real_ else unname(ch$statistic),
             p_value   = if (is.null(ch)) NA_real_ else ch$p.value)
}))
readr::write_csv(chisq_rows, "results/tables/topology_differential_tests.csv")

# ===========================================================================
# (c) DHX36-parallel headline test: is signal LOSS (or gain) topology-specific?
#     Per KO, Fisher 2x2 (topology t vs rest) x (direction vs rest) -> odds ratio.
#     Hypothesis: OR(lost | parallel) > 1 in DHX36KO & dKO (DHX36 is the canonical
#     parallel-G4 resolvase), ~1 in FANCJKO (topology-agnostic).
# ===========================================================================
enrich_rows <- list()
for (k in kos) {
  d <- dz_all %>% filter(ko == k, topology %in% definite,
                         dz_class %in% c("gained", "stable", "lost"))
  for (dir in c("lost", "gained")) {
    for (t in definite) {
      in_topo <- d$topology == t
      is_dir  <- d$dz_class == dir
      a <- sum(in_topo & is_dir);   b <- sum(in_topo & !is_dir)
      cc <- sum(!in_topo & is_dir); dd <- sum(!in_topo & !is_dir)
      ft <- tryCatch(fisher.test(matrix(c(a, b, cc, dd), nrow = 2, byrow = TRUE)),
                     error = function(e) NULL)
      enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
        ko = k, topology = t, direction = dir,
        n_dir_in_topo = a, n_topo = a + b,
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        ci_low  = if (is.null(ft)) NA_real_ else ft$conf.int[1],
        ci_high = if (is.null(ft)) NA_real_ else ft$conf.int[2],
        p_value = if (is.null(ft)) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE)
    }
  }
}
enrich <- dplyr::bind_rows(enrich_rows)
readr::write_csv(enrich, "results/tables/topology_lost_enrichment.csv")

clamp <- function(x) pmin(pmax(x, 1e-2), 1e2)
lost_enr <- enrich %>% filter(direction == "lost") %>%
  mutate(topology = factor(topology, levels = definite),
         ko = factor(ko, levels = kos))
p_enr <- ggplot(lost_enr, aes(topology, clamp(odds_ratio), fill = topology)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = clamp(ci_low), ymax = clamp(ci_high)), width = 0.2) +
  facet_wrap(~ ko) +
  scale_y_log10() +
  scale_fill_manual(values = topology_palette, guide = "none") +
  scale_x_discrete(labels = topo_lab) +
  labs(x = NULL, y = "Odds ratio for 'lost' (log scale)",
       title = "Is G4 signal loss topology-specific? (KO vs WT)",
       subtitle = "OR>1 = 'lost' peaks enriched for this topology (OR<1 = protected from loss)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p_enr, "lost_enrichment_by_topology", "09_gain_loss", fig_root, width = 10, height = 5)

message("Done. Differential class x topology chi-square p-values:")
print(chisq_rows)
message("Lost-enrichment odds ratios (parallel) per KO:")
print(enrich %>% filter(direction == "lost", topology == "parallel") %>%
        select(ko, odds_ratio, p_value))
