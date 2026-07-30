#!/usr/bin/env Rscript
# ============================================================================
# 15_deg_promoter_metaprofiles.R   (env: r_g4)   --- Phase 3, point 2 ---
#
# Promoter-TSS meta-profiles of G4 CUT&Tag signal, split by RNA-seq DEG
# DIRECTION (up vs down, KO-vs-WT), zoomed on the TSS. up = red, down = blue.
# FANCJKO is dropped (too few DEGs -> wide CI, sparse topologies). Each remaining
# KO (DHX36KO, dKO; config metaprofile.deg_grid_kos) is shown as a MATCHED PAIR
# of columns: the WT G4 bigWigs at that KO's DEG promoters (matched control)
# immediately BEFORE the KO's own G4 bigWigs at the SAME promoters. Column order:
#   [WT (DHX36KO DEGs)] [DHX36KO] [WT (dKO DEGs)] [dKO]
# Both views share a fixed 0..deg_fixed_ymax y-axis (cross-panel comparison):
#   (A) facet_grid(topology ~ column) -- direction within each G4 topology.
#   (B) facet_wrap(~ column), topologies pooled -- high-power direction comparison.
#
# DEGs = all significant (padj < deg_padj AND |log2FC| > deg_lfc). Each promoter
# is mapped Entrez -> SYMBOL and joined to DESeq2 by gene_name (same join as
# 13b_promoter_g4_expression.R), and labelled with the topology of its strongest
# overlapping union G4 peak (assign_region_topology()).
#
# Inputs:   results/tables/deseq2_KO_vs_WT.csv (rule 13a),
#           cache_v2 regions_promoter.rds (Entrez names), cache/peaks.rds,
#           results/tables/peak_topology.csv, G4 bigWigs
# Outputs:  results/tables/deg_promoter_profiles.csv
#           results/figures/10_metaprofiles/deg_promoter_topology_grid.{pdf,png}
#           results/figures/10_metaprofiles/deg_promoter_pooled.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles")
fig_root <- "results/figures"

kos        <- setdiff(cfg$genotypes, cfg$ref_genotype)
ref_geno   <- cfg$ref_genotype
# KOs shown with a matched-WT control column (FANCJKO dropped by default).
grid_kos   <- as.character(unlist(cfg$metaprofile$deg_grid_kos %||% setdiff(kos, "FANCJKO")))
half_width <- cfg$metaprofile$deg_half_width %||% 2000
n_bins     <- cfg$metaprofile$n_bins
deg_padj   <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc    <- cfg$metaprofile$deg_lfc %||% 0.5
deg_ymax   <- cfg$metaprofile$deg_fixed_ymax %||% 28    # grid y-cap (pooled plot uses auto y)
max_group  <- cfg$metaprofile$max_per_group %||% 8000
definite   <- c("parallel", "antiparallel", "hybrid")
dir_cols   <- c(up = "#E41A1C", down = "#377EB8")
dir_labels <- c(up = "Upregulated", down = "Downregulated")

# --- Promoters: Entrez -> SYMBOL, topology of strongest G4 peak, TSS window --
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

prom <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$promoter))
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), cfg$std_chroms), pruning.mode = "coarse")
prom_entrez <- names(prom)
if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)
prom$symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_entrez,
                                     column = "SYMBOL", keytype = "ENTREZID",
                                     multiVals = "first")
prom <- assign_region_topology(prom, union)            # adds mcols$group

# Zoomed TSS window (promoter midpoint = TSS; strand carried for - flip).
GenomeInfoDb::seqlengths(prom) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(prom)]
tss <- GenomicRanges::trim(center_window(prom, width = 2 * half_width))
S4Vectors::mcols(tss)$symbol <- prom$symbol
S4Vectors::mcols(tss)$group  <- prom$group
# Keep ALL promoters with a valid symbol (any G4 group, incl. none/mixed/NA) so the
# "unfiltered" set below is a true no-G4-filter pool; the definite/pooled/grid subsets
# still filter on group downstream, so those outputs are unchanged.
tss <- tss[!is.na(tss$symbol)]

de <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = cfg$assay)

# Helper: profile a region set with given bigWigs; NULL if too few regions.
# Promoter subsampling is done ONCE per (group, direction) upstream so the
# matched-WT and KO columns are computed on IDENTICAL promoter sets.
profile_or_null <- function(regions, bws, ko, signal_geno, tp, dir) {
  if (length(regions) < 5) return(NULL)
  prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
  prof$ko <- ko; prof$signal_geno <- signal_geno
  prof$topology <- tp; prof$direction <- dir; prof$n_promoters <- length(regions)
  prof
}

set.seed(7)
wt_bws <- bw_meta$filepath[bw_meta$genotype == ref_geno]   # WT G4 bigWigs (matched control)
all_rows <- list()
for (k in grid_kos) {
  dek <- de %>% dplyr::filter(ko == k, !is.na(padj), !is.na(log2FC))
  dek$direction <- dplyr::case_when(
    dek$padj < deg_padj & dek$log2FC >  deg_lfc ~ "up",
    dek$padj < deg_padj & dek$log2FC < -deg_lfc ~ "down",
    TRUE ~ NA_character_)
  dek <- dek %>% dplyr::filter(!is.na(direction)) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)
  dir_map <- setNames(dek$direction, dek$gene_name)

  tk <- tss
  S4Vectors::mcols(tk)$direction <- unname(dir_map[tk$symbol])
  tk <- tk[!is.na(tk$direction)]
  ko_bws <- bw_meta$filepath[bw_meta$genotype == k]
  if (length(tk) == 0 || length(ko_bws) == 0 || length(wt_bws) == 0) {
    message("Skipping ", k, " (", length(tk), " DEG promoters, ",
            length(ko_bws), " KO bigWigs, ", length(wt_bws), " WT bigWigs)")
    next
  }
  message("DEG promoter profiles: ", k, "; ",
          sum(tk$direction == "up"), " up / ", sum(tk$direction == "down"),
          " down promoters; matched WT + KO signal")

  # Fix the promoter set per (topology|pooled|unfiltered, direction) ONCE, then profile
  # both WT and the KO over the identical regions (true matched control).
  #   "pooled"     -> G4+ (definite) topologies only: a true weighted average of the grid
  #                   rows; the no-G4-peak ("none") group excluded so it cannot dilute.
  #   "unfiltered" -> ALL DEG promoters (any group incl. none/mixed): no G4 filter at all.
  # "unfiltered" is appended LAST so the RNG draw order for the existing definite/pooled
  # sets is untouched (defensive; no set trips max_group here anyway).
  region_sets <- list()
  for (tp in c("pooled", definite, "unfiltered")) for (dir in c("up", "down")) {
    regs <- if (tp == "unfiltered") tk[tk$direction == dir]
            else if (tp == "pooled") tk[tk$direction == dir & tk$group %in% definite]
            else tk[tk$direction == dir & tk$group == tp]
    if (length(regs) > max_group) regs <- regs[sample(length(regs), max_group)]
    region_sets[[paste(tp, dir, sep = "@@")]] <- regs
  }

  # Column order in the figures: WT-matched control THEN the KO itself.
  for (sg in c(ref_geno, k)) {
    bws <- if (sg == ref_geno) wt_bws else ko_bws
    for (key in names(region_sets)) {
      parts <- strsplit(key, "@@", fixed = TRUE)[[1]]
      all_rows <- c(all_rows, list(
        profile_or_null(region_sets[[key]], bws, k, sg, parts[1], parts[2])))
    }
    if (sg != ref_geno) evict_bigwigs(ko_bws)   # keep WT cached across KOs
  }
}
evict_bigwigs(wt_bws)

all_rows <- all_rows[!vapply(all_rows, is.null, logical(1))]
profiles <- dplyr::bind_rows(all_rows)
if (nrow(profiles) == 0) stop("No DEG-promoter profiles computed (check DESeq2 output / joins).")

# Column = matched-WT control ("WT (<KO> DEGs)") or the KO's own signal (<KO>),
# ordered WT-then-KO per grid KO -> [WT (DHX36KO DEGs)] [DHX36KO] [WT (dKO DEGs)] [dKO]
col_levels <- unlist(lapply(grid_kos, function(g) c(sprintf("WT (%s DEGs)", g), g)))
profiles$column <- ifelse(profiles$signal_geno == ref_geno,
                          sprintf("WT (%s DEGs)", profiles$ko),
                          as.character(profiles$signal_geno))
profiles$column    <- factor(profiles$column, levels = col_levels)
profiles$ko        <- factor(profiles$ko, levels = grid_kos)
profiles$direction <- factor(profiles$direction, levels = c("up", "down"))
readr::write_csv(profiles, "results/tables/deg_promoter_profiles.csv")

sub_base <- sprintf(paste0("+/- %d bp, %d bins; significant DEGs (padj<%g & |log2FC|>%g); ",
                           "WT-matched control vs KO G4 signal"),
                    half_width, n_bins, deg_padj, deg_lfc)

# --- (B) pooled figures (auto y-axis; per-panel up/down n annotation) --------
# Per-column "up N / down M" label (matched WT and KO columns share identical sets).
panel_counts <- function(df) {
  up <- df %>% dplyr::filter(direction == "up") %>%
    dplyr::distinct(column, n_promoters) %>% dplyr::rename(n_up = n_promoters)
  dn <- df %>% dplyr::filter(direction == "down") %>%
    dplyr::distinct(column, n_promoters) %>% dplyr::rename(n_down = n_promoters)
  m <- merge(up, dn, by = "column", all = TRUE)
  m$label <- sprintf("up %s / down %s",
                     ifelse(is.na(m$n_up), "0", as.character(m$n_up)),
                     ifelse(is.na(m$n_down), "0", as.character(m$n_down)))
  m$column <- factor(m$column, levels = levels(df$column))
  m
}
build_pooled <- function(topo_key, fname, title, subtitle) {
  df <- profiles %>% dplyr::filter(topology == topo_key)
  p <- ggplot(df, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = panel_counts(df), aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey20") +
    facet_wrap(~ column, nrow = 1) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = "Distance to TSS (bp)", y = "Mean G4 CUT&Tag signal",
         title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, "10_metaprofiles", fig_root, width = 13, height = 4)
}
build_pooled("pooled", "deg_promoter_pooled",
             "Promoter G4 signal at DEGs by direction (G4+ topologies pooled)",
             sprintf("%s; parallel+antiparallel+hybrid pooled (no-G4-peak excluded); auto y-axis", sub_base))
build_pooled("unfiltered", "deg_promoter_pooled_unfiltered",
             "Promoter G4 signal at DEGs by direction (all DEGs - no G4 filter)",
             sprintf("%s; ALL DEG promoters incl. no-G4-peak and mixed; auto y-axis", sub_base))

# --- (A) topology x column grid ---------------------------------------------
grid <- profiles %>% dplyr::filter(topology %in% definite)
grid$topology <- factor(grid$topology, levels = definite)
topo_labs <- setNames(topology_labels[match(definite, topology_levels)], definite)
p_grid <- ggplot(grid, aes(position, mean, colour = direction, fill = direction)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(topology ~ column, labeller = labeller(topology = topo_labs)) +
  coord_cartesian(ylim = c(0, deg_ymax)) +
  scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  labs(x = "Distance to TSS (bp)", y = "Mean G4 CUT&Tag signal",
       title = "Promoter G4 signal at DEGs by direction and G4 topology",
       subtitle = sprintf("%s; shared y-axis 0-%g", sub_base, deg_ymax)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_grid, "deg_promoter_topology_grid", "10_metaprofiles", fig_root,
          width = 13, height = 8)

# --- (C) standalone gene-count summary (how many promoters feed each panel) --
filter_levels <- c("unfiltered", "pooled", "parallel", "antiparallel", "hybrid")
filter_labels <- c(unfiltered = "All DEGs", pooled = "G4+ (pooled)", parallel = "Parallel",
                   antiparallel = "Antiparallel", hybrid = "Hybrid")
counts <- profiles %>%
  dplyr::distinct(ko, topology, direction, n_promoters) %>%
  dplyr::rename(filter = topology, n = n_promoters) %>%
  dplyr::filter(filter %in% filter_levels)
counts$filter <- factor(counts$filter, levels = filter_levels)
readr::write_csv(counts, "results/tables/deg_promoter_counts.csv")

p_counts <- ggplot(counts, aes(filter, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.8) +
  facet_wrap(~ ko, nrow = 1) +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_x_discrete(labels = filter_labels) +
  labs(x = NULL, y = "DEG promoters (n)",
       title = "DEG promoters per metaprofile panel",
       subtitle = "All DEGs = no-G4 filter; G4+ (pooled) = parallel+antiparallel+hybrid") +
  theme_pub() + theme(legend.position = "bottom",
                      axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_counts, "deg_promoter_counts", "10_metaprofiles", fig_root, width = 10, height = 5)

message("Done. DEG-promoter meta-profiles (WT-matched + KO). Promoters per (ko, topology, direction):")
print(profiles %>% dplyr::distinct(ko, topology, direction, n_promoters) %>%
        dplyr::arrange(ko, topology, direction))
