#!/usr/bin/env Rscript
# ============================================================================
# 16_deg_5utr_metaprofiles.R   (env: r_g4)   --- Phase 4, point 2 ---
#
# 5'UTR analog of rule 15 (DEG-direction promoter profiles): G4 CUT&Tag signal
# centred on 5'UTRs, split by RNA-seq DEG DIRECTION (up=red / down=blue) and by
# G4 topology, with a MATCHED-WT control column before each KO. Column order:
#   [WT (DHX36KO DEGs)] [DHX36KO] [WT (dKO DEGs)] [dKO]   (FANCJKO dropped)
#
# Run for TWO region units so they can be compared side-by-side:
#   - "fragment": every 5'UTR exon piece as-is (regions_5UTR.rds; ~42.5k, multi
#                 fragments per gene -> may inflate the centred peak).
#   - "gene":     one TSS-proximal 5'UTR fragment per gene (collapse_to_tss_proximal;
#                 ~19.5k, one window per gene).
#
# Each 5'UTR is mapped mcols$gene_id (Entrez) -> SYMBOL (org.Mm.eg.db) and joined
# to DESeq2 by gene_name (same join as 13b/15), and labelled with the topology of
# its strongest overlapping union G4 peak (assign_region_topology()).
#
# Inputs:   results/tables/deseq2_KO_vs_WT.csv (rule 13a), cache/peaks.rds,
#           results/tables/peak_topology.csv, cache_v2 regions_5UTR.rds, G4 bigWigs
# Outputs:  results/tables/deg_5utr_profiles.csv  (adds `unit` column)
#           results/figures/10_metaprofiles/deg_5utr_{fragment,gene}_topology_grid.{pdf,png}
#           results/figures/10_metaprofiles/deg_5utr_{fragment,gene}_pooled.{pdf,png}
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

kos          <- setdiff(cfg$genotypes, cfg$ref_genotype)
ref_geno     <- cfg$ref_genotype
grid_kos     <- as.character(unlist(cfg$metaprofile$deg_grid_kos %||% setdiff(kos, "FANCJKO")))
half_width   <- cfg$metaprofile$deg_half_width %||% 2000
n_bins       <- cfg$metaprofile$n_bins
deg_padj     <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc      <- cfg$metaprofile$deg_lfc %||% 0.5
grid_ymax    <- cfg$metaprofile$deg_5utr_grid_ymax %||% 40   # grid y-cap (pooled uses auto y)
max_group    <- cfg$metaprofile$max_per_group %||% 8000
min_regions  <- max(1L, as.integer(cfg$metaprofile$min_regions_per_panel %||% 1))
definite     <- topology_definite          # v1.1: shared, was a local copy
dir_cols     <- c(up = "#E41A1C", down = "#377EB8")
dir_labels   <- c(up = "Upregulated", down = "Downregulated")

# --- Union G4 peaks with topology ------------------------------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])

# v1.1: grid rows follow the labels actually present. "mixed" is G4SP's abstention
# class and terminal now that rule 21 no longer back-fills it; a hard-coded
# `definite` computed it in no region set, so those 5'UTRs left the figure without
# reaching any denominator. Measured 2026-08-11: 19-69 per (unit, ko, direction),
# against 0-1 antiparallel and 0-8 hybrid. "none" is excluded -- the no-G4-peak
# 5'UTRs are already the difference between the "pooled" and "unfiltered" panels.
keep_groups <- topology_keep_groups(union$topology, incl_none = FALSE)

# --- Build topology/symbol-labelled, centre-anchored 5'UTR windows per unit --
utr5_raw <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$`5UTR`))
GenomeInfoDb::seqlevelsStyle(utr5_raw) <- "UCSC"
utr5_raw <- GenomeInfoDb::keepSeqlevels(
  utr5_raw, intersect(GenomeInfoDb::seqlevels(utr5_raw), cfg$std_chroms),
  pruning.mode = "coarse")

build_unit_windows <- function(gr) {
  entrez <- as.character(S4Vectors::mcols(gr)$gene_id)
  gr$symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = entrez,
                                     column = "SYMBOL", keytype = "ENTREZID",
                                     multiVals = "first")
  gr <- assign_region_topology(gr, union)                 # topology on the small UTR piece
  GenomeInfoDb::seqlengths(gr) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(gr)]
  win <- GenomicRanges::trim(center_window(gr, width = 2 * half_width))
  S4Vectors::mcols(win)$symbol <- gr$symbol
  S4Vectors::mcols(win)$group  <- gr$group
  # Keep ALL 5'UTRs with a valid symbol (any G4 group incl. none/mixed/NA) so the
  # "unfiltered" set is a true no-G4-filter pool; definite/pooled/grid subsets still
  # filter on group downstream, so those outputs are unchanged.
  win <- win[!is.na(win$symbol)]
  win
}

unit_windows <- list(
  fragment = build_unit_windows(utr5_raw),
  gene     = build_unit_windows(collapse_to_tss_proximal(utr5_raw, "gene_id")))
for (u in names(unit_windows))
  message("5'UTR unit '", u, "': ", length(unit_windows[[u]]), " windows; group table:")
print(lapply(unit_windows, function(w) table(w$group)))

de      <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = cfg$assay)
wt_bws  <- bw_meta$filepath[bw_meta$genotype == ref_geno]

# Profile a region set with given bigWigs; NULL if too few regions. (Subsampling
# is done ONCE per (unit, group, direction) so matched-WT and KO use identical sets.)
profile_or_null <- function(regions, bws, unit, ko, signal_geno, tp, dir) {
  if (length(regions) < min_regions) return(NULL)
  prof <- mean_replicate_profile(regions, bws, n_bins = n_bins, half_width = half_width)
  prof$unit <- unit; prof$ko <- ko; prof$signal_geno <- signal_geno
  prof$topology <- tp; prof$direction <- dir; prof$n_5utr <- length(regions)
  prof
}

set.seed(7)
all_rows <- list(); count_rows <- list()
for (k in grid_kos) {
  dek <- de %>% dplyr::filter(ko == k, !is.na(padj), !is.na(log2FC))
  dek$direction <- dplyr::case_when(
    dek$padj < deg_padj & dek$log2FC >  deg_lfc ~ "up",
    dek$padj < deg_padj & dek$log2FC < -deg_lfc ~ "down",
    TRUE ~ NA_character_)
  dek <- dek %>% dplyr::filter(!is.na(direction)) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)
  dir_map <- setNames(dek$direction, dek$gene_name)

  ko_bws <- bw_meta$filepath[bw_meta$genotype == k]
  if (length(ko_bws) == 0 || length(wt_bws) == 0) {
    message("Skipping ", k, " (", length(ko_bws), " KO bigWigs, ", length(wt_bws), " WT bigWigs)")
    next
  }

  # Fix the 5'UTR set per (unit, topology|pooled|unfiltered, direction) ONCE, then
  # profile both WT and the KO over identical regions (true matched control).
  #   "pooled"     -> G4+ (definite) topologies only: weighted average of grid rows
  #                   (no-G4-peak "none" excluded so it cannot dilute).
  #   "unfiltered" -> ALL DEG 5'UTRs (any group incl. none/mixed): no G4 filter.
  # "unfiltered" and the extra keep_groups strata ("mixed") are appended LAST so
  # existing sets' RNG draw order is untouched.
  region_sets <- list()   # key "unit@@tp@@dir"
  for (u in names(unit_windows)) {
    wu <- unit_windows[[u]]
    S4Vectors::mcols(wu)$direction <- unname(dir_map[wu$symbol])
    wu <- wu[!is.na(wu$direction)]
    for (tp in c("pooled", definite, "unfiltered", setdiff(keep_groups, definite)))
      for (dir in c("up", "down")) {
        regs <- if (tp == "unfiltered") wu[wu$direction == dir]
                else if (tp == "pooled") wu[wu$direction == dir & wu$group %in% definite]
                else wu[wu$direction == dir & wu$group == tp]
        if (length(regs) > max_group) regs <- regs[sample(length(regs), max_group)]
        region_sets[[paste(u, tp, dir, sep = "@@")]] <- regs
        # True n recorded BEFORE the min_regions guard, so a stratum too small to
        # profile still reports its size instead of vanishing from the counts table.
        count_rows[[length(count_rows) + 1]] <-
          data.frame(unit = u, ko = k, filter = tp, direction = dir, n = length(regs),
                     stringsAsFactors = FALSE)
      }
    message("DEG 5'UTR (", u, ") ", k, ": ",
            sum(wu$direction == "up"), " up / ", sum(wu$direction == "down"), " down 5'UTRs")
  }

  # Column order: WT-matched control THEN the KO; both units profiled while the
  # genotype's bigWigs are cached (load each genotype once).
  for (sg in c(ref_geno, k)) {
    bws <- if (sg == ref_geno) wt_bws else ko_bws
    for (key in names(region_sets)) {
      parts <- strsplit(key, "@@", fixed = TRUE)[[1]]
      all_rows <- c(all_rows, list(
        profile_or_null(region_sets[[key]], bws, parts[1], k, sg, parts[2], parts[3])))
    }
    if (sg != ref_geno) evict_bigwigs(ko_bws)   # keep WT cached across KOs
  }
}
evict_bigwigs(wt_bws)

all_rows <- all_rows[!vapply(all_rows, is.null, logical(1))]
profiles <- dplyr::bind_rows(all_rows)
if (nrow(profiles) == 0) stop("No DEG-5'UTR profiles computed (check DESeq2 output / joins).")

col_levels <- unlist(lapply(grid_kos, function(g) c(sprintf("WT (%s DEGs)", g), g)))
profiles$column <- ifelse(profiles$signal_geno == ref_geno,
                          sprintf("WT (%s DEGs)", profiles$ko),
                          as.character(profiles$signal_geno))
profiles$column    <- factor(profiles$column, levels = col_levels)
profiles$ko        <- factor(profiles$ko, levels = grid_kos)
profiles$direction <- factor(profiles$direction, levels = c("up", "down"))
profiles$unit      <- factor(profiles$unit, levels = c("fragment", "gene"))
readr::write_csv(profiles, "results/tables/deg_5utr_profiles.csv")

# --- Stratum sizes, including the strata too small to profile ---------------
# Built from the region sets, NOT from `profiles`: a stratum below min_regions is
# never profiled, so deriving the counts from the profiles would have deleted the
# very rows a reader needs in order to see why a facet is empty. Defined here
# rather than beside its figure because the per-unit grid annotation reads it.
filter_levels <- c("unfiltered", "pooled", keep_groups)
filter_labels <- c(unfiltered = "All DEGs", pooled = "G4+ (pooled)",
                   topology_group_labels(keep_groups))
counts <- dplyr::bind_rows(count_rows) %>%
  dplyr::filter(filter %in% filter_levels)
counts$filter    <- factor(counts$filter, levels = filter_levels)
counts$unit      <- factor(counts$unit, levels = c("fragment", "gene"))
counts$ko        <- factor(counts$ko, levels = grid_kos)
counts$direction <- factor(counts$direction, levels = c("up", "down"))
readr::write_csv(counts, "results/tables/deg_5utr_counts.csv")

# The matched-WT column and the KO column are profiled over IDENTICAL 5'UTR sets,
# so one count serves both; expand each ko into its two grid columns.
grid_n <- counts %>% dplyr::filter(filter %in% keep_groups) %>%
  dplyr::rename(topology = filter)
grid_n <- dplyr::bind_rows(
  grid_n %>% dplyr::mutate(column = sprintf("WT (%s DEGs)", ko)),
  grid_n %>% dplyr::mutate(column = as.character(ko)))
grid_n$column   <- factor(grid_n$column, levels = col_levels)
grid_n$topology <- factor(as.character(grid_n$topology), levels = keep_groups)

topo_labs <- topology_group_labels(keep_groups)
unit_titles <- c(fragment = "all 5'UTR fragments",
                 gene     = "one TSS-proximal 5'UTR per gene")
sub_base <- sprintf(paste0("+/- %d bp, %d bins; significant DEGs (padj<%g & |log2FC|>%g); ",
                           "WT-matched control vs KO G4 signal"),
                    half_width, n_bins, deg_padj, deg_lfc)

# Per-column "up N / down M" label (matched WT and KO columns share identical sets).
# panel_n_labels() lives in helpers.R -- the grid below needs the same label keyed
# on (topology, column), so a local copy would have been duplicated.
build_pooled <- function(df, fname, title, subtitle) {
  p <- ggplot(df, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = panel_n_labels(df, "column", "n_5utr"),
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey20") +
    facet_wrap(~ column, nrow = 1) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = "Distance to 5'UTR centre (bp)", y = "Mean G4 CUT&Tag signal",
         title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, "10_metaprofiles", fig_root, width = 13, height = 4)
}

# --- Per-unit figures: grid (topology x column) + pooled (G4+ & unfiltered) --
for (u in levels(profiles$unit)) {
  pu <- profiles %>% dplyr::filter(unit == u)

  # drop = FALSE keeps a stratum's row visible even when it holds no 5'UTRs at all,
  # so "no data" reads differently from "no signal"; the n annotation says which.
  grid <- pu %>% dplyr::filter(topology %in% keep_groups)
  grid$topology <- factor(grid$topology, levels = keep_groups)
  p_grid <- ggplot(grid, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = panel_n_labels(grid_n %>% dplyr::filter(unit == u),
                                    c("topology", "column"), "n"),
              aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 2.8, colour = "grey20") +
    facet_grid(topology ~ column, drop = FALSE,
               labeller = labeller(topology = topo_labs)) +
    coord_cartesian(ylim = c(0, grid_ymax)) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = "Distance to 5'UTR centre (bp)", y = "Mean G4 CUT&Tag signal",
         title = sprintf("5'UTR G4 signal at DEGs by direction and G4 topology (%s)", unit_titles[u]),
         subtitle = sprintf("%s; shared y-axis 0-%g (lines above %g are clipped); n per panel",
                            sub_base, grid_ymax, grid_ymax)) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p_grid, sprintf("deg_5utr_%s_topology_grid", u), "10_metaprofiles",
            fig_root, width = 13, height = 8)

  build_pooled(pu %>% dplyr::filter(topology == "pooled"),
               sprintf("deg_5utr_%s_pooled", u),
               sprintf("5'UTR G4 signal at DEGs by direction (G4+ topologies pooled, %s)", unit_titles[u]),
               sprintf("%s; parallel+antiparallel+hybrid pooled (no-G4-peak excluded); auto y-axis", sub_base))
  build_pooled(pu %>% dplyr::filter(topology == "unfiltered"),
               sprintf("deg_5utr_%s_pooled_unfiltered", u),
               sprintf("5'UTR G4 signal at DEGs by direction (all DEGs - no G4 filter, %s)", unit_titles[u]),
               sprintf("%s; ALL DEG 5'UTRs incl. no-G4-peak and mixed; auto y-axis", sub_base))
}

# --- Standalone region-count summary (how many 5'UTRs feed each panel) ------
p_counts <- ggplot(counts, aes(filter, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5) +
  facet_grid(unit ~ ko, labeller = labeller(unit = unit_titles)) +
  scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
  scale_x_discrete(labels = filter_labels) +
  labs(x = NULL, y = "DEG 5'UTRs (n)",
       title = "DEG 5'UTRs per metaprofile panel",
       subtitle = paste("All DEGs = no-G4 filter; G4+ (pooled) = parallel+antiparallel+hybrid;",
                        "counted before the min-n profiling guard")) +
  theme_pub() + theme(legend.position = "bottom",
                      axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_counts, "deg_5utr_counts", "10_metaprofiles", fig_root, width = 11, height = 7)

message("Done. DEG-5'UTR meta-profiles (fragment + gene). 5'UTRs per (unit, ko, filter, direction):")
print(counts %>% dplyr::arrange(unit, ko, filter, direction))
skipped <- counts %>% dplyr::filter(n < min_regions)
if (nrow(skipped) > 0) {
  message("Strata below min_regions_per_panel = ", min_regions,
          " (facet shown empty, no curve drawn):")
  print(skipped)
}
