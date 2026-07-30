#!/usr/bin/env Rscript
# ============================================================================
# 25_g4_rloop_promoter_correlation.R   (env: r_g4)
#
# Co-occurrence question: at a given promoter, does high G4 (BG4) CUT&Tag
# signal coincide with high R-loop (S9.6) CUT&Tag signal? Run PER GENOTYPE
# (WT, DHX36KO, FANCJKO, dKO) so the coupling can be compared across the
# helicase knockouts. (The earlier expression-vs-signal test gave a validated
# ~0 correlation, so we pivot to correlating the two CUT&Tag assays.)
#
# Two correlations, computed for each genotype:
#   1. mean_promoter -- one mean G4 + one mean R-loop value per promoter
#      (TSS +/- mean_half_width), correlated across ALL promoters.
#   2. binned        -- each promoter split into bin_size-bp bins over
#      TSS +/- bin_half_width (default 200 bp x 10 = TSS +/- 1 kb); correlate
#      G4 vs R-loop at bin resolution. Heavier -> runs on a seeded promoter
#      subsample. Also reports the distribution of per-promoter r across bins.
#
# Normalisation -- making the two assays directly comparable:
#   Raw bigWig coverage differs vastly between G4 and R-loop (antibody strength,
#   sequencing depth). We depth-normalise EACH replicate to its OWN genomic
#   background (background_median: divide by the median signal of ~20k random
#   std-chrom windows) -> comparable "x background" units; mean across that
#   genotype's reps; log2(+pc). Then we z-score each assay WITHIN a genotype so
#   both axes are 0-centred/unit-var and a y=x diagonal is meaningful. The
#   z-score is monotone-affine, so Pearson r / Spearman rho (reported on the
#   log2 fold-over-background values) are UNCHANGED by it -- it exists purely to
#   put G4 and R-loop on a shared, overlayable axis.
#
# NB on the binned pooled r: it is partly driven by the shared TSS meta-shape
#   (both assays rise toward the TSS), so it captures BOTH between-promoter level
#   co-variation AND within-promoter positional co-variation. The
#   binned_per_promoter summary isolates the within-promoter component.
#
# Inputs:   cfg$paths$promoters_rds (regions_promoter.rds, names = Entrez),
#           per-genotype G4 + R-loop bigWigs (cfg$paths$bigwig_dir).
# Outputs:  results/tables/g4_rloop_promoter_signal.csv               (per promoter x genotype)
#           results/tables/g4_rloop_promoter_correlation_stats.csv    (r/rho/p x genotype)
#           results/figures/18_g4_rloop_correlation/
#             g4_vs_rloop_promoter_mean.{pdf,png}          (facet by genotype)
#             g4_vs_rloop_promoter_binned.{pdf,png}        (facet by genotype)
#             g4_vs_rloop_binned_per_promoter_hist.{pdf,png} (facet by genotype)
#             g4_vs_rloop_correlation_summary.{pdf,png}    (r/rho across genotypes)
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/18_g4_rloop_correlation", "cache")
fig_root <- "results/figures"

gc_cfg     <- cfg$g4_rloop_correlation %||% list()
genotypes  <- as.character(unlist(gc_cfg$genotypes %||% cfg$genotypes))
mean_hw    <- as.integer(gc_cfg$mean_half_width %||% 2000L)
bin_hw     <- as.integer(gc_cfg$bin_half_width  %||% 1000L)
bin_size   <- as.integer(gc_cfg$bin_size        %||% 200L)
pc_cnt     <- gc_cfg$pseudocount     %||% 1
cnt_norm   <- gc_cfg$cnt_norm        %||% "background_median"
bg_n       <- as.integer(gc_cfg$background_n    %||% 20000L)
bg_seed    <- as.integer(gc_cfg$background_seed %||% 7L)
sub_n      <- as.integer(gc_cfg$bin_subsample_n %||% 3000L)
comparab   <- gc_cfg$comparability   %||% "zscore"
use_hexbin <- isTRUE(gc_cfg$hexbin %||% TRUE)
if (use_hexbin && !requireNamespace("hexbin", quietly = TRUE)) {
  message("Package 'hexbin' unavailable; mean-analysis plot falls back to a scatter.")
  use_hexbin <- FALSE
}
n_bins <- as.integer((2L * bin_hw) %/% bin_size)
stopifnot(n_bins >= 2L)

std_chroms <- cfg$std_chroms
assays     <- c(G4_BG4 = "G4_BG4", Rloop_S96 = "Rloop_S96")

# z-score (NA-safe); identity if comparability != "zscore"
zscore <- function(x) {
  if (!identical(comparab, "zscore")) return(x)
  mu <- mean(x, na.rm = TRUE); sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(x - mu)
  (x - mu) / sdv
}

# ---------------------------------------------------------------------------
# Promoters (TSS +/- mean_hw) on std chroms; attach SYMBOL for the table
# ---------------------------------------------------------------------------
prom <- readRDS(cfg$paths$promoters_rds)
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), std_chroms), pruning.mode = "coarse")
# regions_promoter.rds is already TSS +/- 2 kb; only re-centre if a different
# mean_half_width is requested.
if (mean_hw != 2000L) prom <- center_window(prom, width = 2L * mean_hw)

prom_entrez <- names(prom)
if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)
names(prom) <- prom_entrez
message("Promoters on std chroms: ", length(prom))

prom_symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_entrez,
                                     column = "SYMBOL", keytype = "ENTREZID",
                                     multiVals = "first")

# Shared subsample (same promoters for every genotype -> fair cross-genotype comparison)
set.seed(bg_seed)
sub_ix   <- if (length(prom) > sub_n) sort(sample.int(length(prom), sub_n)) else seq_along(prom)
prom_sub <- prom[sub_ix]
win      <- center_window(prom_sub, width = 2L * bin_hw)

# ---------------------------------------------------------------------------
# Per-replicate background level (depth proxy; reused for both analyses)
# ---------------------------------------------------------------------------
random_windows <- function(n, width, chrom_sizes, seed) {
  set.seed(seed)
  sz  <- chrom_sizes[names(chrom_sizes) %in% std_chroms]
  chr <- sample(names(sz), n, replace = TRUE, prob = as.numeric(sz))
  st  <- floor(runif(n, min = 1, max = pmax(2, as.numeric(sz[chr]) - width)))
  gr  <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = st, width = width))
  names(gr) <- sprintf("bg_%07d", seq_len(n))
  gr
}

chrom_sizes <- default_chrom_sizes(std = std_chroms)
bg_gr       <- random_windows(bg_n, 2L * mean_hw, chrom_sizes, bg_seed)

bg_level_for <- function(fp) {
  bg <- compute_region_signal(bg_gr, fp)
  lvl <- stats::median(bg[is.finite(bg) & bg > 0], na.rm = TRUE)
  if (!is.finite(lvl) || lvl <= 0) lvl <- 1
  lvl
}

depth_normalise <- function(sig_mat, fps) {
  if (identical(cnt_norm, "background_median")) {
    for (j in seq_len(ncol(sig_mat))) sig_mat[, j] <- sig_mat[, j] / bg_level_for(fps[j])
  } else if (identical(cnt_norm, "replicate_median")) {
    meds <- apply(sig_mat, 2, function(x) stats::median(x[is.finite(x) & x > 0], na.rm = TRUE))
    meds[!is.finite(meds) | meds <= 0] <- 1
    target <- mean(meds)
    for (j in seq_len(ncol(sig_mat))) sig_mat[, j] <- sig_mat[, j] * target / meds[j]
  } else {
    stop("Unknown cnt_norm: ", cnt_norm)
  }
  sig_mat
}

bw_for <- function(geno, assay) {
  bw <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = assay)
  bw <- bw[as.character(bw$genotype) == geno, , drop = FALSE]
  if (nrow(bw) == 0) stop("No ", geno, " bigWigs for assay ", assay)
  bw
}

# Mean per-promoter depth-normalised signal for one (genotype, assay), cached.
mean_signal <- function(geno, assay) {
  cache_file <- file.path("cache", sprintf("g4rl_sig_mean_%s_%s_%s.rds",
                                           tolower(geno), assay, cnt_norm))
  cache_or_build(cache_file, {
    bw <- bw_for(geno, assay)
    message(geno, " / ", assay, " (mean): ", nrow(bw), " replicate(s)")
    sig_mat <- vapply(bw$filepath, function(fp) compute_region_signal(prom, fp),
                      numeric(length(prom)))
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    sig_mat <- depth_normalise(sig_mat, bw$filepath)
    evict_bigwigs(bw$filepath)
    rowMeans(sig_mat, na.rm = TRUE)
  })
}

# Binned (promoter x bin) mean depth-normalised matrix for one (genotype, assay).
bin_signal <- function(geno, assay) {
  cache_file <- file.path("cache", sprintf("g4rl_sig_bin_%s_%s_%s.rds",
                                           tolower(geno), assay, cnt_norm))
  cache_or_build(cache_file, {
    bw <- bw_for(geno, assay)
    message(geno, " / ", assay, " (binned): ", nrow(bw), " replicate(s), ",
            length(win), " promoters x ", n_bins, " bins")
    mats <- lapply(seq_len(nrow(bw)), function(j) {
      m <- compute_profile_matrix(win, bw$filepath[j], n_bins = n_bins,
                                  half_width = bin_hw)$bin_mat
      m / bg_level_for(bw$filepath[j])
    })
    evict_bigwigs(bw$filepath)
    arr <- simplify2array(mats)
    if (length(dim(arr)) == 3) apply(arr, c(1, 2), mean, na.rm = TRUE) else arr
  })
}

# ---------------------------------------------------------------------------
# Correlation helper
# ---------------------------------------------------------------------------
corr_row <- function(x, y, geno, analysis) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3)
    return(data.frame(genotype = geno, analysis = analysis, n = length(x),
                      pearson_r = NA, pearson_p = NA, spearman_rho = NA, spearman_p = NA))
  pe <- suppressWarnings(stats::cor.test(x, y, method = "pearson"))
  sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  data.frame(genotype = geno, analysis = analysis, n = length(x),
             pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
             spearman_rho = unname(sp$estimate), spearman_p = sp$p.value)
}

# ---------------------------------------------------------------------------
# Per-genotype analysis
# ---------------------------------------------------------------------------
analyse_genotype <- function(geno) {
  # --- mean per promoter ---
  g4_m <- mean_signal(geno, assays[["G4_BG4"]])
  rl_m <- mean_signal(geno, assays[["Rloop_S96"]])
  per_prom <- data.frame(entrez = names(prom), symbol = unname(prom_symbol),
                         genotype = geno, g4_mean = g4_m, rloop_mean = rl_m,
                         stringsAsFactors = FALSE) %>%
    dplyr::mutate(log_g4    = log2(g4_mean + pc_cnt),
                  log_rloop = log2(rloop_mean + pc_cnt),
                  z_g4      = zscore(log_g4),
                  z_rloop   = zscore(log_rloop))

  # --- binned ---
  g4_b <- bin_signal(geno, assays[["G4_BG4"]])
  rl_b <- bin_signal(geno, assays[["Rloop_S96"]])
  bin_long <- data.frame(
    genotype = geno,
    prom_i = rep(seq_len(nrow(g4_b)), times = ncol(g4_b)),
    bin    = rep(seq_len(ncol(g4_b)), each  = nrow(g4_b)),
    g4     = as.numeric(g4_b), rloop = as.numeric(rl_b),
    stringsAsFactors = FALSE) %>%
    dplyr::mutate(log_g4 = log2(g4 + pc_cnt), log_rloop = log2(rloop + pc_cnt)) %>%
    dplyr::filter(is.finite(log_g4), is.finite(log_rloop))
  bin_long$z_g4    <- zscore(bin_long$log_g4)
  bin_long$z_rloop <- zscore(bin_long$log_rloop)

  per_prom_r <- bin_long %>%
    dplyr::group_by(prom_i) %>%
    dplyr::summarise(
      r = if (dplyr::n() >= 3 && stats::sd(log_g4) > 0 && stats::sd(log_rloop) > 0)
            suppressWarnings(stats::cor(log_g4, log_rloop)) else NA_real_,
      .groups = "drop") %>%
    dplyr::filter(is.finite(r)) %>%
    dplyr::mutate(genotype = geno)

  stats <- dplyr::bind_rows(
    corr_row(per_prom$log_g4, per_prom$log_rloop, geno, "mean_promoter"),
    corr_row(bin_long$log_g4, bin_long$log_rloop, geno, "binned_pooled"),
    data.frame(genotype = geno, analysis = "binned_per_promoter", n = nrow(per_prom_r),
               pearson_r = stats::median(per_prom_r$r, na.rm = TRUE), pearson_p = NA_real_,
               spearman_rho = NA_real_, spearman_p = NA_real_)
  )
  list(per_prom = per_prom, bin_long = bin_long, per_prom_r = per_prom_r, stats = stats)
}

res <- lapply(genotypes, analyse_genotype)
names(res) <- genotypes

geno_lvls    <- intersect(main_genotypes, genotypes)   # WT, DHX36KO, FANCJKO, dKO order
mk_factor    <- function(df) { df$genotype <- factor(df$genotype, levels = geno_lvls); df }
per_prom_all <- mk_factor(dplyr::bind_rows(lapply(res, `[[`, "per_prom")))
bin_all      <- mk_factor(dplyr::bind_rows(lapply(res, `[[`, "bin_long")))
per_prom_r_all <- mk_factor(dplyr::bind_rows(lapply(res, `[[`, "per_prom_r")))
stats_tbl    <- mk_factor(dplyr::bind_rows(lapply(res, `[[`, "stats")))

readr::write_csv(per_prom_all, "results/tables/g4_rloop_promoter_signal.csv")
readr::write_csv(stats_tbl,    "results/tables/g4_rloop_promoter_correlation_stats.csv")

# ---------------------------------------------------------------------------
# Plots (faceted by genotype)
# ---------------------------------------------------------------------------
xlab_z <- sprintf("z-scored log2 %s promoter signal", assay_label_map[["G4_BG4"]])
ylab_z <- sprintf("z-scored log2 %s promoter signal", assay_label_map[["Rloop_S96"]])

# per-genotype facet labels carrying r / rho / n
facet_lab <- function(st_sub) {
  setNames(
    sprintf("%s\nr=%.2f  rho=%.2f  n=%s", st_sub$genotype, st_sub$pearson_r,
            st_sub$spearman_rho, format(st_sub$n, big.mark = ",")),
    as.character(st_sub$genotype))
}
lab_mean <- facet_lab(stats_tbl[stats_tbl$analysis == "mean_promoter", ])
lab_bin  <- facet_lab(stats_tbl[stats_tbl$analysis == "binned_pooled", ])

# --- mean-per-promoter (hexbin), facet by genotype ---
p_mean <- ggplot(per_prom_all, aes(x = z_g4, y = z_rloop))
p_mean <- p_mean + (if (use_hexbin)
  list(geom_hex(bins = 60), scale_fill_viridis_c(trans = "log10", name = "promoters"))
  else geom_point(alpha = 0.25, size = 0.4, colour = "grey25"))
p_mean <- p_mean +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE,
              colour = "#E41A1C", fill = "#E41A1C", linewidth = 0.9) +
  facet_wrap(~ genotype, labeller = labeller(genotype = lab_mean)) +
  labs(x = xlab_z, y = ylab_z,
       title = sprintf("G4 vs R-loop promoter signal by genotype (mean over TSS +/- %d bp)", mean_hw),
       subtitle = "Pearson r / Spearman rho per panel; dashed = y=x, red = loess") +
  theme_pub()
save_plot(p_mean, "g4_vs_rloop_promoter_mean", "18_g4_rloop_correlation", fig_root, 9, 8)

# --- binned pooled (scatter), facet by genotype ---
p_bin <- ggplot(bin_all, aes(x = z_g4, y = z_rloop)) +
  geom_point(alpha = 0.25, size = 0.4, colour = "grey25") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_smooth(method = "loess", span = 0.6, se = TRUE,
              colour = "#E41A1C", fill = "#E41A1C", linewidth = 0.9) +
  facet_wrap(~ genotype, labeller = labeller(genotype = lab_bin)) +
  labs(x = xlab_z, y = ylab_z,
       title = sprintf("G4 vs R-loop signal in %d-bp bins by genotype (TSS +/- %d bp, %d promoters)",
                       bin_size, bin_hw, length(prom_sub)),
       subtitle = "pooled r partly reflects the shared TSS profile shape") +
  theme_pub()
save_plot(p_bin, "g4_vs_rloop_promoter_binned", "18_g4_rloop_correlation", fig_root, 9, 8)

# --- per-promoter r distribution, facet by genotype ---
med_r_df <- per_prom_r_all %>% dplyr::group_by(genotype) %>%
  dplyr::summarise(med = stats::median(r, na.rm = TRUE), .groups = "drop")
p_hist <- ggplot(per_prom_r_all, aes(x = r)) +
  geom_histogram(bins = 50, fill = "#377EB8", colour = "white") +
  geom_vline(data = med_r_df, aes(xintercept = med), linetype = 2,
             colour = "#E41A1C", linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey40") +
  facet_wrap(~ genotype) +
  labs(x = "within-promoter Pearson r (G4 vs R-loop across bins)", y = "promoters",
       title = "Within-promoter spatial co-variation of G4 and R-loop by genotype",
       subtitle = sprintf("red dashed = median r per genotype (each promoter has %d bins)", n_bins)) +
  theme_pub()
save_plot(p_hist, "g4_vs_rloop_binned_per_promoter_hist", "18_g4_rloop_correlation", fig_root, 9, 7)

# --- cross-genotype summary of correlation coefficients ---
summ <- stats_tbl %>%
  tidyr::pivot_longer(c(pearson_r, spearman_rho), names_to = "metric", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(metric = dplyr::recode(metric, pearson_r = "Pearson r",
                                       spearman_rho = "Spearman rho"),
                analysis = factor(analysis,
                                  levels = c("mean_promoter", "binned_pooled", "binned_per_promoter")))
p_sum <- ggplot(summ, aes(x = genotype, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_wrap(~ analysis) +
  scale_fill_manual(values = c("Pearson r" = "#E41A1C", "Spearman rho" = "#377EB8"), name = NULL) +
  labs(x = NULL, y = "correlation coefficient",
       title = "G4-R-loop promoter coupling across genotypes",
       subtitle = "binned_per_promoter = median within-promoter Pearson r (no Spearman)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_sum, "g4_vs_rloop_correlation_summary", "18_g4_rloop_correlation", fig_root, 10, 5)

message("Done. G4 vs R-loop promoter-signal correlations (per genotype):")
print(stats_tbl)
