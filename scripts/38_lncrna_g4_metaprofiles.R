#!/usr/bin/env Rscript
# ============================================================================
# 23_lncrna_g4_metaprofiles.R   (env: system R 4.5.1)   --- Phase 9 ---
#
# G4 CUT&Tag signal at lncRNA TSS (GENCODE lncRNA genes from rule 21), two views:
#
#  View 1 -- GAINED / LOST.  Each lncRNA TSS inherits the gained/stable/lost class
#    (rule-07 topology_differential.csv dz_class) of its strongest overlapping
#    union G4 peak. For every KO we profile WT *and* KO G4 signal over the same
#    lncRNA-TSS sets -> facet_grid(dz_class ~ ko), lines = genotype. "lost" facets
#    => the KO line drops below WT; "gained" => rises.
#
#  View 2 -- DEG up / down (mirror rule 15).  lncRNA TSS joined to DESeq2 by
#    unversioned ENSEMBL gene_id (recovers lncRNAs the rule-15 symbol-join drops);
#    classify up/down per KO (grid_kos: DHX36KO, dKO); matched WT-vs-KO columns
#    [WT (<KO> DEGs)] [<KO>]; topologies pooled (lncRNA G4 overlap is sparse).
#    up = #E41A1C, down = #377EB8; auto y-axis; per-column up/down n annotation.
#
# Efficiency: lncRNA TSS windows are built ONCE; each genotype's G4 bigWigs are
# loaded once and profiled over every region set that needs that genotype's
# signal (12 bigWig loads total), then evicted.
#
# Inputs:   cache/lncrna_genes.rds (rule 21), cache/peaks.rds,
#           results/tables/peak_topology.csv, results/tables/topology_differential.csv
#           (rule 07), results/tables/deseq2_KO_vs_WT.csv (rule 13a), G4 bigWigs
# Outputs:  results/tables/lncrna_g4_{gainloss,deg}_profiles.csv, lncrna_g4_counts.csv
#           results/figures/17_lncrna/{lncrna_g4_gainloss,lncrna_g4_deg_pooled,
#                                      lncrna_g4_counts}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2); library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/17_lncrna")
fig_root <- "results/figures"

ref_geno   <- cfg$ref_genotype
kos        <- setdiff(cfg$genotypes, ref_geno)
grid_kos   <- as.character(unlist(cfg$metaprofile$deg_grid_kos %||% setdiff(kos, "FANCJKO")))
n_bins     <- cfg$metaprofile$n_bins
win_half   <- cfg$metaprofile$deg_half_width %||% cfg$metaprofile$half_width %||% 2000
deg_padj   <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc    <- cfg$metaprofile$deg_lfc  %||% 0.5
max_group  <- cfg$metaprofile$max_per_group %||% 8000
dz_levels  <- c("gained", "stable", "lost")
dir_cols   <- c(up = "#E41A1C", down = "#377EB8")
dir_labels <- c(up = "Upregulated", down = "Downregulated")
strip_ver  <- function(x) sub("\\..*$", "", as.character(x))

# --- lncRNA TSS windows (strand-aware; mean_replicate_profile flips '-') -----
lnc <- readRDS("cache/lncrna_genes.rds")
GenomeInfoDb::seqlevelsStyle(lnc) <- "UCSC"
lnc <- GenomeInfoDb::keepSeqlevels(
  lnc, intersect(GenomeInfoDb::seqlevels(lnc), cfg$std_chroms), pruning.mode = "coarse")
GenomeInfoDb::seqlengths(lnc) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(lnc)]
tss <- GenomicRanges::trim(GenomicRanges::promoters(lnc, upstream = win_half, downstream = win_half))
tss$gene_id <- strip_ver(names(lnc))
message(sprintf("lncRNA TSS windows: %d (+/- %d bp)", length(tss), win_half))

# --- Strongest overlapping union peak per lncRNA TSS (by max_z) --------------
union <- readRDS("cache/peaks.rds")$union
mz <- S4Vectors::mcols(union)$max_z
if (is.null(mz)) mz <- rep(0, length(union))
mz[is.na(mz)] <- -Inf
sp_idx <- rep(NA_integer_, length(tss))
ov <- GenomicRanges::findOverlaps(tss, union, ignore.strand = TRUE)
if (length(ov) > 0) {
  qh <- S4Vectors::queryHits(ov); sh <- S4Vectors::subjectHits(ov)
  o  <- order(qh, -mz[sh]); qo <- qh[o]; so <- sh[o]
  first <- !duplicated(qo); sp_idx[qo[first]] <- so[first]
}
tss$peak_id <- names(union)[sp_idx]   # NA where no overlapping G4 peak
message(sprintf("lncRNA TSS with an overlapping G4 peak: %d / %d",
                sum(!is.na(tss$peak_id)), length(tss)))

dz  <- readr::read_csv("results/tables/topology_differential.csv", show_col_types = FALSE)
de  <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
de$gene_id <- strip_ver(de$gene_id)
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = cfg$assay)

# ===========================================================================
# Build region sets (fixed ONCE), tagged with which genotype signals they need.
# ===========================================================================
set.seed(7)
cap <- function(g) if (length(g) > max_group) g[sample(length(g), max_group)] else g
region_sets <- list()   # key -> list(regs, view, ko, class, signal_genos)

# View 1: gained/stable/lost (per KO), profiled with WT + that KO.
for (k in kos) {
  dzk <- dz %>% dplyr::filter(ko == k)
  dz_map <- setNames(as.character(dzk$dz_class), dzk$peak_id)
  cls <- unname(dz_map[tss$peak_id])
  for (d in dz_levels) {
    regs <- cap(tss[!is.na(cls) & cls == d])
    if (length(regs) < 5) { message("  skip gainloss ", k, "/", d, " (", length(regs), ")"); next }
    region_sets[[paste("gainloss", k, d, sep = "@@")]] <-
      list(regs = regs, view = "gainloss", ko = k, class = d,
           signal_genos = c(ref_geno, k))
  }
}

# View 2: DEG up/down (grid KOs), matched WT + KO; lncRNA pooled (no topology).
for (k in grid_kos) {
  dek <- de %>% dplyr::filter(ko == k, !is.na(padj), !is.na(log2FC))
  dek$direction <- dplyr::case_when(
    dek$padj < deg_padj & dek$log2FC >  deg_lfc ~ "up",
    dek$padj < deg_padj & dek$log2FC < -deg_lfc ~ "down",
    TRUE ~ NA_character_)
  dir_map <- setNames(dek$direction, dek$gene_id)
  cls <- unname(dir_map[tss$gene_id])
  for (d in c("up", "down")) {
    regs <- cap(tss[!is.na(cls) & cls == d])
    if (length(regs) < 5) { message("  skip deg ", k, "/", d, " (", length(regs), ")"); next }
    region_sets[[paste("deg", k, d, sep = "@@")]] <-
      list(regs = regs, view = "deg", ko = k, class = d,
           signal_genos = c(ref_geno, k))
  }
}
# View 0: GLOBAL -- all lncRNA TSS (no G4/DEG/topology filter), every genotype's
# signal over the identical set. Not capped: "all lncRNA" is the point, and
# compute_profile_matrix is vectorised so the full ~13-14k set is cheap.
region_sets[["global@@all@@all"]] <- list(
  regs = tss, view = "global", ko = "all", class = "all",
  signal_genos = c(ref_geno, kos))

if (length(region_sets) == 0) stop("No lncRNA region sets passed the >=5 guard (check joins).")

# ===========================================================================
# Profile: load each genotype's bigWigs once, cover all sets needing it, evict.
# ===========================================================================
profile_set <- function(regs, bws, view, ko, class, signal_geno) {
  prof <- mean_replicate_profile(regs, bws, n_bins = n_bins, half_width = win_half)
  prof$view <- view; prof$ko <- ko; prof$class <- class
  prof$signal_geno <- signal_geno; prof$n <- length(regs)
  prof
}
all_rows <- list()
for (g in c(ref_geno, kos)) {
  bws <- bw_meta$filepath[bw_meta$genotype == g]
  if (length(bws) == 0) next
  jobs <- Filter(function(s) g %in% s$signal_genos, region_sets)
  if (length(jobs) == 0) next
  message("Profiling ", g, " G4 signal over ", length(jobs), " lncRNA region sets")
  for (s in jobs)
    all_rows <- c(all_rows, list(profile_set(s$regs, bws, s$view, s$ko, s$class, g)))
  evict_bigwigs(bws)
}
profiles <- dplyr::bind_rows(all_rows)

# ===========================================================================
# View 0 figure + table: GLOBAL G4 signal at all lncRNA TSS, lines = genotype
# ===========================================================================
v0 <- profiles %>% dplyr::filter(view == "global")
v0$signal_geno <- factor(v0$signal_geno, levels = c(ref_geno, kos))
readr::write_csv(
  v0 %>% dplyr::select(genotype = signal_geno, bin, position, mean, sem, n),
  "results/tables/lncrna_g4_global_profiles.csv")

p_v0 <- ggplot(v0, aes(position, mean, colour = signal_geno, fill = signal_geno)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = condition_colours, name = "Genotype") +
  scale_fill_manual(values = condition_colours, name = "Genotype") +
  labs(x = "Distance to lncRNA TSS (bp)", y = "Mean G4 CUT&Tag signal",
       title = "Global G4 signal at all lncRNA TSS",
       subtitle = sprintf("all %d GENCODE lncRNA TSS +/- %d bp; lines = genotype",
                          length(tss), win_half)) +
  theme_pub() + theme(legend.position = "right")
save_plot(p_v0, "lncrna_g4_global_tss", "17_lncrna", fig_root, width = 7, height = 5)

# ===========================================================================
# View 1 figure + table: gained/stable/lost, WT vs KO lines
# ===========================================================================
v1 <- profiles %>% dplyr::filter(view == "gainloss")
v1$class <- factor(v1$class, levels = dz_levels)
v1$ko    <- factor(v1$ko, levels = kos)
v1$signal_geno <- factor(v1$signal_geno, levels = c(ref_geno, kos))
readr::write_csv(
  v1 %>% dplyr::select(ko, dz_class = class, genotype = signal_geno, bin, position, mean, sem, n),
  "results/tables/lncrna_g4_gainloss_profiles.csv")

p_v1 <- ggplot(v1, aes(position, mean, colour = signal_geno, fill = signal_geno)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_grid(class ~ ko) +
  scale_colour_manual(values = condition_colours, name = "Genotype") +
  scale_fill_manual(values = condition_colours, name = "Genotype") +
  labs(x = "Distance to lncRNA TSS (bp)", y = "Mean G4 CUT&Tag signal",
       title = "G4 signal at lncRNA TSS by differential class (KO vs WT)",
       subtitle = sprintf("lncRNA TSS +/- %d bp; dz_class = gained/stable/lost of strongest overlapping G4 peak (|log2FC|>%.2f, adj.P<%.2g)",
                          win_half, cfg$differential$voom_lfc, cfg$differential$voom_padj)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_v1, "lncrna_g4_gainloss", "17_lncrna", fig_root, width = 11, height = 8)

# ===========================================================================
# View 2 figure + table: DEG up/down, matched WT-vs-KO columns (auto y)
# ===========================================================================
v2 <- profiles %>% dplyr::filter(view == "deg")
counts_out <- NULL
if (nrow(v2) > 0) {
  v2$column <- ifelse(v2$signal_geno == ref_geno,
                      sprintf("WT (%s DEGs)", v2$ko), as.character(v2$signal_geno))
  col_levels <- unlist(lapply(grid_kos, function(g) c(sprintf("WT (%s DEGs)", g), g)))
  v2$column    <- factor(v2$column, levels = col_levels)
  v2$direction <- factor(v2$class, levels = c("up", "down"))
  readr::write_csv(
    v2 %>% dplyr::select(ko, column, direction, bin, position, mean, sem, n),
    "results/tables/lncrna_g4_deg_profiles.csv")

  panel_counts <- function(df) {
    up <- df %>% dplyr::filter(direction == "up") %>%
      dplyr::distinct(column, n) %>% dplyr::rename(n_up = n)
    dn <- df %>% dplyr::filter(direction == "down") %>%
      dplyr::distinct(column, n) %>% dplyr::rename(n_down = n)
    m <- merge(up, dn, by = "column", all = TRUE)
    m$label <- sprintf("up %s / down %s",
                       ifelse(is.na(m$n_up), "0", as.character(m$n_up)),
                       ifelse(is.na(m$n_down), "0", as.character(m$n_down)))
    m$column <- factor(m$column, levels = levels(df$column)); m
  }
  p_v2 <- ggplot(v2, aes(position, mean, colour = direction, fill = direction)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = panel_counts(v2), aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.08, vjust = 1.4, size = 3, colour = "grey20") +
    facet_wrap(~ column, nrow = 1) +
    scale_colour_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = "Distance to lncRNA TSS (bp)", y = "Mean G4 CUT&Tag signal",
         title = "G4 signal at lncRNA TSS by DEG direction (matched WT vs KO)",
         subtitle = sprintf("lncRNA DEGs (padj<%g & |log2FC|>%g) joined by ENSEMBL gene_id; topologies pooled; auto y-axis",
                            deg_padj, deg_lfc)) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p_v2, "lncrna_g4_deg_pooled", "17_lncrna", fig_root, width = 13, height = 4)
} else {
  message("No lncRNA DEG region sets passed the guard; skipping View 2 figure.")
}

# ===========================================================================
# Counts figure + table (how many lncRNA TSS feed each panel)
# ===========================================================================
c1 <- v1 %>% dplyr::filter(signal_geno == ref_geno) %>%
  dplyr::distinct(ko, category = class, n) %>% dplyr::mutate(view = "gained/lost")
c2 <- NULL
if (nrow(v2) > 0) {
  c2 <- v2 %>% dplyr::filter(signal_geno == ref_geno) %>%
    dplyr::distinct(ko, category = direction, n) %>% dplyr::mutate(view = "DEG up/down")
}
counts_out <- dplyr::bind_rows(c1, c2)
readr::write_csv(counts_out %>% dplyr::select(view, ko, category, n),
                 "results/tables/lncrna_g4_counts.csv")

dz_pal  <- c(gained = "#E41A1C", stable = "grey70", lost = "#377EB8")
p_c1 <- ggplot(c1, aes(ko, n, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = dz_pal, name = "vs WT") +
  labs(x = NULL, y = "lncRNA TSS (n)", title = "Gained/lost view") +
  theme_pub() + theme(legend.position = "bottom")
counts_panels <- list(p_c1)
if (!is.null(c2)) {
  p_c2 <- ggplot(c2, aes(ko, n, fill = category)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
    scale_fill_manual(values = dir_cols, labels = dir_labels, name = "DEG direction") +
    labs(x = NULL, y = "lncRNA DEG TSS (n)", title = "DEG up/down view") +
    theme_pub() + theme(legend.position = "bottom")
  counts_panels <- list(p_c1, p_c2)
}
save_plot(cowplot::plot_grid(plotlist = counts_panels, nrow = 1),
          "lncrna_g4_counts", "17_lncrna", fig_root, width = 11, height = 5)

message("Done. lncRNA G4 metaprofiles (gained/lost + DEG up/down) written.")
