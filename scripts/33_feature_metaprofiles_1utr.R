#!/usr/bin/env Rscript
# ============================================================================
# 17_feature_metaprofiles_1utr.R   (env: r_g4)   --- Phase 4, point 3 ---
#
# Reproduces the rule-14 feature-region G4 meta-profile (promoter / 5'UTR / 1st
# intron / enhancer x topology, all 4 genotypes) but with the 5'UTR collapsed to
# ONE TSS-proximal fragment per gene (collapse_to_tss_proximal). This tests the
# hypothesis that the strong 5'UTR peak in rule 14 is inflated by the ~2.2
# fragments/gene in regions_5UTR.rds.
#
# EFFICIENT: recomputes profiles ONLY for the collapsed 5'UTR, then reuses rule
# 14's promoter/intron1/enhancer rows verbatim from feature_topology_profiles.csv
# (exact apples-to-apples; the only thing that changes between the two figures is
# the 5'UTR definition). The collapsed-5'UTR profiles are cached under a DISTINCT
# key so they never collide with rule 14's cache.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           cache_v2 regions_5UTR.rds, results/tables/feature_topology_profiles.csv
#           (rule 14, for the reused feature rows), G4 bigWigs
# Outputs:  results/tables/feature_topology_profiles_1utr.csv
#           results/figures/10_metaprofiles/feature_signal_by_topology_1utr.{pdf,png}
#           results/figures/10_metaprofiles/feature_signal_by_topology_1utr_fixedY.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles", "cache")
fig_root <- "results/figures"

genotypes  <- cfg$genotypes
half_width <- cfg$metaprofile$half_width
n_bins     <- cfg$metaprofile$n_bins
incl_none  <- isTRUE(cfg$metaprofile$include_none_group)
max_group  <- cfg$metaprofile$max_per_group %||% 8000
definite   <- c("parallel", "antiparallel", "hybrid")
keep_groups <- c(definite, if (incl_none) "none")
feature_regions <- as.character(unlist(cfg$metaprofile$feature_regions %||%
                                       c("promoter", "5UTR", "intron1", "enhancer")))

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

# --- Collapsed-5'UTR windows (one TSS-proximal fragment per gene) -----------
set.seed(7)
utr5 <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$`5UTR`))
GenomeInfoDb::seqlevelsStyle(utr5) <- "UCSC"
utr5 <- GenomeInfoDb::keepSeqlevels(
  utr5, intersect(GenomeInfoDb::seqlevels(utr5), cfg$std_chroms), pruning.mode = "coarse")
utr5 <- collapse_to_tss_proximal(utr5, "gene_id")
utr5 <- assign_region_topology(utr5, union)
utr5 <- utr5[utr5$group %in% keep_groups]
stopifnot(length(utr5) > 0)
GenomeInfoDb::seqlengths(utr5) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(utr5)]
win <- GenomicRanges::trim(center_window(utr5, width = 2 * half_width))
S4Vectors::mcols(win)$group <- utr5$group
keep_idx <- unlist(lapply(keep_groups, function(g) {
  idx <- which(win$group == g)
  if (length(idx) > max_group) sample(idx, max_group) else idx
}))
win <- win[sort(keep_idx)]
message("Collapsed 5'UTR windows per group:"); print(table(win$group))

# --- Per-genotype profiles for the collapsed 5'UTR (distinct cache key) ------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)
prof_cache <- sprintf("cache/feature_profiles_5utr1pg_n%d_hw%d_b%d.rds",
                      length(union), half_width, n_bins)
utr_profiles <- cache_or_build(prof_cache, {
  dplyr::bind_rows(lapply(genotypes, function(g) {
    bws <- bw_meta$filepath[bw_meta$genotype == g]
    if (length(bws) == 0) return(NULL)
    message("Collapsed-5'UTR profiles: ", g, " (", length(bws), " bigWigs)")
    rows <- dplyr::bind_rows(lapply(keep_groups, function(grp) {
      regions <- win[win$group == grp]
      if (length(regions) < 5) return(NULL)
      prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
      prof$genotype  <- g
      prof$region    <- "5UTR"
      prof$group     <- grp
      prof$n_regions <- length(regions)
      prof
    }))
    evict_bigwigs(bws)
    rows
  }))
})

# --- Reuse rule-14 rows for promoter/intron1/enhancer, swap in collapsed 5'UTR
csv14 <- "results/tables/feature_topology_profiles.csv"
if (!file.exists(csv14))
  stop("rule-14 output missing — run `feature_metaprofiles` first: ", csv14)
other <- readr::read_csv(csv14, show_col_types = FALSE) %>% dplyr::filter(region != "5UTR")
reused <- setdiff(feature_regions, "5UTR")
missing <- setdiff(reused, unique(other$region))
if (length(missing) > 0)
  stop("rule-14 CSV lacks reused feature rows: ", paste(missing, collapse = ", "))

profiles <- dplyr::bind_rows(other, utr_profiles)
profiles$genotype <- factor(profiles$genotype, levels = genotypes)
profiles$region   <- factor(profiles$region, levels = feature_regions)
profiles$group    <- factor(profiles$group, levels = keep_groups)
readr::write_csv(profiles, "results/tables/feature_topology_profiles_1utr.csv")

# --- Plot (identical to rule 14; only the 5'UTR row differs) -----------------
region_labels_map <- c("promoter" = "Promoter", "5UTR" = "5'UTR (1 per gene)",
                       "intron1" = "1st intron", "enhancer" = "Enhancer (dELS)",
                       "genebody" = "Gene body")
present <- intersect(feature_regions, unique(as.character(profiles$region)))
profiles$region_label <- factor(
  unname(region_labels_map[as.character(profiles$region)]),
  levels = unname(region_labels_map[present]))
group_labels <- c(setNames(topology_labels[match(definite, topology_levels)], definite),
                  none = "No G4 peak")

base_layers <- list(
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA),
  geom_line(linewidth = 0.8),
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40"),
  scale_colour_manual(values = condition_colours),
  scale_fill_manual(values = condition_colours),
  labs(x = "Distance to feature centre (bp)", y = "Mean G4 CUT&Tag signal",
       colour = "Genotype", fill = "Genotype"),
  theme_pub(), theme(legend.position = "bottom"))
fig_h <- 2.2 * length(present) + 1.5

p_free <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_grid(region_label ~ group, scales = "free_y",
             labeller = labeller(group = group_labels)) +
  labs(title = "Feature-region G4 signal by topology — 5'UTR collapsed to 1 per gene",
       subtitle = sprintf("+/- %d bp, %d bins; per-panel y-axis (promoter/intron1/enhancer reused from rule 14)",
                          half_width, n_bins))
save_plot(p_free, "feature_signal_by_topology_1utr", "10_metaprofiles", fig_root,
          width = 11, height = fig_h)

ymax <- cfg$metaprofile$feature_fixed_ymax %||% 30
p_fixed <- ggplot(profiles, aes(position, mean, colour = genotype, fill = genotype)) +
  base_layers +
  facet_grid(region_label ~ group, labeller = labeller(group = group_labels)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(title = "Feature-region G4 signal by topology — 5'UTR collapsed to 1 per gene (common y-axis)",
       subtitle = sprintf("+/- %d bp, %d bins; shared y-axis 0-%g (compare the 5'UTR row against rule 14's per-fragment plot)",
                          half_width, n_bins, ymax))
save_plot(p_fixed, "feature_signal_by_topology_1utr_fixedY", "10_metaprofiles", fig_root,
          width = 11, height = fig_h)

message("Done. 1-UTR-per-gene feature meta-profiles (free + fixed y). Collapsed 5'UTR genes: ",
        length(unique(names(utr5))))
