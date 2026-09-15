#!/usr/bin/env Rscript
# ============================================================================
# 26_g4_rloop_cooccurrence.R   (env: r_g4)
#
# Promoter-level CO-OCCURRENCE of G4 (BG4) and R-loop (S9.6) CUT&Tag signal,
# complementing rule 25 (which only correlates mean signal). Four views, all
# per genotype (WT + the 3 helicase KOs), on depth-normalised signal:
#
#   1. Sorted dual heatmap   -- rows = promoters sorted by G4 signal; aligned
#      G4 | R-loop heatmap panels over TSS +/- half_width. Per-locus visual proof.
#   2. Stratified metaprofile -- split promoters into G4 quartiles, overlay the
#      mean R-loop TSS metaprofile per quartile (and the reciprocal). Dose-
#      dependent spatial co-occurrence.
#   2b. TSS signal-vs-signal scatter, per quartile -- the same strata, but
#      plotted as G4-at-TSS vs R-loop-at-TSS with one hex-binned panel per
#      (genotype x quartile). The metaprofile averages promoters together and so
#      destroys the promoter-to-promoter pairing; this keeps it, and answers
#      "inside a stratum, does a promoter with more G4 also carry more R-loop?"
#   3. G4+/- promoter violin -- R-loop signal at G4-peak+ vs G4-peak- promoters
#      (and reciprocal), vs a random-window null. Wilcoxon + Cliff's delta.
#   4. dG4 vs dR-loop KO views -- per-promoter differential (voom, KO vs WT)
#      scatter + an R-loop volcano coloured by G4 differential class. Tests
#      whether helicase loss moves the two marks together.
#
# Normalisation (same as rule 25): depth-normalise each replicate to its own
# genomic background (background_median), mean across a genotype's reps. Scalar
# per-promoter signal is log2(+pc); profile matrices stay in linear "x background"
# units for the metaprofiles and are z-scored per assay for the heatmap fill.
#
# Inputs:   cfg$paths$promoters_rds (names = Entrez), cache/peaks.rds (G4 union),
#           cfg$rloop$union_rds (R-loop union), per-genotype G4 + R-loop bigWigs.
# Outputs:  results/tables/g4_rloop_stratified_profiles.csv
#           results/tables/g4_rloop_violin_stats.csv
#           results/tables/g4_rloop_ko_differential.csv
#           results/tables/g4_rloop_tss_scores.csv
#           results/tables/g4_rloop_tss_quartile_correlation.csv
#           results/figures/19_g4_rloop_cooccurrence/
#             g4_rloop_promoter_heatmap.{pdf,png}
#             g4_rloop_stratified_rloopByG4.{pdf,png}
#             g4_rloop_stratified_g4ByRloop.{pdf,png}
#             g4_rloop_tss_scatter_byG4quartile.{pdf,png}
#             g4_rloop_tss_scatter_byRloopquartile.{pdf,png}
#             g4_rloop_promoter_violin.{pdf,png}
#             g4_rloop_delta_scatter.{pdf,png}
#             g4_rloop_rloop_volcano.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(scales)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
fig_subdir <- "19_g4_rloop_cooccurrence"
ensure_dirs("results/tables", file.path("results/figures", fig_subdir), "cache")
fig_root <- "results/figures"

cc        <- cfg$g4_rloop_cooccurrence %||% list()
genotypes <- as.character(unlist(cc$genotypes %||% cfg$genotypes))
ref_geno  <- cc$ref_genotype %||% cfg$ref_genotype %||% "WT"
half_width <- as.integer(cc$half_width %||% 2000L)
n_bins    <- as.integer(cc$n_bins %||% 100L)
pc_cnt    <- cc$pseudocount %||% 1
cnt_norm  <- cc$cnt_norm %||% "background_median"
bg_n      <- as.integer(cc$background_n %||% 20000L)
bg_seed   <- as.integer(cc$background_seed %||% 7L)
n_strata  <- as.integer(cc$n_strata %||% 4L)
heat_genos <- as.character(unlist(cc$heatmap_genotypes %||% list("WT")))
heat_max  <- cc$heatmap_max_rows %||% 5000L
heat_cap  <- cc$heatmap_cap_z %||% 3
voom_padj <- cc$voom_padj %||% (cfg$differential$voom_padj %||% 0.1)
voom_lfc  <- cc$voom_lfc  %||% (cfg$differential$voom_lfc  %||% 0.585)
prof_sub  <- cc$profile_subsample_n %||% NULL
tss_hw    <- as.integer(cc$tss_scatter_half_width %||% 250L)
hex_bins  <- as.integer(cc$hexbin_bins %||% 60L)
tss_hex_bins <- as.integer(cc$tss_scatter_hexbin_bins %||% 45L)
tss_fill_lim <- as.numeric(unlist(cc$tss_scatter_fill_limits %||% list(1, 100)))
tss_fill_brk <- as.numeric(unlist(cc$tss_scatter_fill_breaks %||% list(1, 10, 100)))
stopifnot(length(tss_fill_lim) == 2, tss_fill_lim[1] > 0, tss_fill_lim[2] > tss_fill_lim[1])
volc_hex_cls <- as.character(unlist(cc$volcano_hex_classes %||% list("stable", "ambiguous")))
use_hexbin <- isTRUE(cc$hexbin %||% TRUE)
if (use_hexbin && !requireNamespace("hexbin", quietly = TRUE)) {
  message("Package 'hexbin' unavailable; the scatter views fall back to points. ",
          "install.packages('hexbin') to get small, Illustrator-editable PDFs.")
  use_hexbin <- FALSE
}

std_chroms <- cfg$std_chroms
assays     <- c(G4_BG4 = "G4_BG4", Rloop_S96 = "Rloop_S96")
kos        <- setdiff(genotypes, ref_geno)
geno_lvls  <- intersect(main_genotypes, genotypes)
positions  <- seq(-half_width, half_width, length.out = n_bins)
g4_lab     <- assay_label_map[["G4_BG4"]]
rl_lab     <- assay_label_map[["Rloop_S96"]]

# ---------------------------------------------------------------------------
# Promoters (TSS +/- half_width), std chroms, names = Entrez
# ---------------------------------------------------------------------------
prom <- readRDS(cfg$paths$promoters_rds)
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), std_chroms), pruning.mode = "coarse")
if (half_width != 2000L) prom <- center_window(prom, width = 2L * half_width)
prom_entrez <- names(prom)
if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)
names(prom) <- prom_entrez
prom_symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_entrez,
                                     column = "SYMBOL", keytype = "ENTREZID",
                                     multiVals = "first")
message("Promoters on std chroms: ", length(prom))

# Promoters used for the profile views (optionally subsampled for speed)
if (!is.null(prof_sub) && length(prom) > prof_sub) {
  set.seed(bg_seed); prof_ix <- sort(sample.int(length(prom), as.integer(prof_sub)))
} else prof_ix <- seq_along(prom)
prom_prof <- prom[prof_ix]

# ---------------------------------------------------------------------------
# Shared depth-normalisation helpers (identical to rule 25)
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
bg_gr       <- random_windows(bg_n, 2L * half_width, chrom_sizes, bg_seed)

bg_level_for <- function(fp) {
  bg <- compute_region_signal(bg_gr, fp)
  lvl <- stats::median(bg[is.finite(bg) & bg > 0], na.rm = TRUE)
  if (!is.finite(lvl) || lvl <= 0) lvl <- 1
  lvl
}
depth_norm_vec <- function(v, fp) v / bg_level_for(fp)

bw_for <- function(geno, assay) {
  bw <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = assay)
  bw <- bw[as.character(bw$genotype) == geno, , drop = FALSE]
  if (nrow(bw) == 0) stop("No ", geno, " bigWigs for assay ", assay)
  bw
}

# Scalar mean per-promoter signal over ALL promoters (shares rule-25 cache).
mean_signal <- function(geno, assay) {
  cache_file <- file.path("cache", sprintf("g4rl_sig_mean_%s_%s_%s.rds",
                                           tolower(geno), assay, cnt_norm))
  cache_or_build(cache_file, {
    bw <- bw_for(geno, assay)
    message(geno, " / ", assay, " (mean): ", nrow(bw), " replicate(s)")
    sig_mat <- vapply(bw$filepath, function(fp) compute_region_signal(prom, fp),
                      numeric(length(prom)))
    if (is.null(dim(sig_mat))) sig_mat <- matrix(sig_mat, ncol = 1)
    if (identical(cnt_norm, "background_median")) {
      for (j in seq_len(ncol(sig_mat))) sig_mat[, j] <- depth_norm_vec(sig_mat[, j], bw$filepath[j])
    } else if (identical(cnt_norm, "replicate_median")) {
      meds <- apply(sig_mat, 2, function(x) stats::median(x[is.finite(x) & x > 0], na.rm = TRUE))
      meds[!is.finite(meds) | meds <= 0] <- 1; target <- mean(meds)
      for (j in seq_len(ncol(sig_mat))) sig_mat[, j] <- sig_mat[, j] * target / meds[j]
    } else stop("Unknown cnt_norm: ", cnt_norm)
    evict_bigwigs(bw$filepath)
    rowMeans(sig_mat, na.rm = TRUE)
  })
}

# Depth-normalised mean profile matrix (prom_prof x n_bins) per (geno, assay).
profile_signal <- function(geno, assay) {
  cache_file <- file.path("cache", sprintf("cooc_prof_%s_%s.rds", tolower(geno), assay))
  cache_or_build(cache_file, {
    bw <- bw_for(geno, assay)
    message(geno, " / ", assay, " (profile): ", nrow(bw), " replicate(s), ",
            length(prom_prof), " promoters x ", n_bins, " bins")
    mats <- lapply(seq_len(nrow(bw)), function(j) {
      m <- compute_profile_matrix(prom_prof, bw$filepath[j], n_bins = n_bins,
                                  half_width = half_width)$bin_mat
      m / bg_level_for(bw$filepath[j])
    })
    evict_bigwigs(bw$filepath)
    arr <- simplify2array(mats)
    if (length(dim(arr)) == 3) apply(arr, c(1, 2), mean, na.rm = TRUE) else arr
  })
}

# ===========================================================================
# View 1 + 2 need profile matrices for all (genotype, assay)
# ===========================================================================
prof <- list()
for (g in genotypes) for (a in assays) prof[[paste(g, a, sep = "|")]] <- profile_signal(g, a)
get_prof <- function(g, a) prof[[paste(g, a, sep = "|")]]

# ---------------------------------------------------------------------------
# View 1: sorted dual heatmap (heat_genos)
# ---------------------------------------------------------------------------
zscore_mat <- function(m) {
  mu <- mean(m, na.rm = TRUE); sdv <- stats::sd(as.numeric(m), na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  (m - mu) / sdv
}

heat_long <- dplyr::bind_rows(lapply(intersect(heat_genos, genotypes), function(g) {
  zg <- zscore_mat(get_prof(g, "G4_BG4"))
  zr <- zscore_mat(get_prof(g, "Rloop_S96"))
  ord <- order(rowMeans(zg, na.rm = TRUE), decreasing = TRUE)   # sort by G4
  if (!is.null(heat_max) && length(ord) > heat_max)
    ord <- ord[round(seq(1, length(ord), length.out = as.integer(heat_max)))]
  zg <- zg[ord, , drop = FALSE]; zr <- zr[ord, , drop = FALSE]
  mk <- function(z, lab) data.frame(
    genotype = g, assay = lab,
    rank = rep(seq_len(nrow(z)), times = ncol(z)),
    position = rep(positions, each = nrow(z)),
    z = as.numeric(z), stringsAsFactors = FALSE)
  dplyr::bind_rows(mk(zg, g4_lab), mk(zr, rl_lab))
}))
heat_long$genotype <- factor(heat_long$genotype, levels = geno_lvls)
heat_long$assay    <- factor(heat_long$assay, levels = c(g4_lab, rl_lab))

p_heat <- ggplot(heat_long, aes(position, rank, fill = z)) +
  geom_raster() +
  facet_grid(genotype ~ assay, switch = "y") +
  scale_fill_viridis_c(limits = c(-heat_cap, heat_cap), oob = scales::squish,
                       name = "signal (z)") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "distance from TSS (bp)", y = "promoters (sorted by G4 signal)",
       title = "G4 and R-loop signal at promoters, rows sorted by G4",
       subtitle = "same row order in both panels -> shared rows = co-occurrence") +
  theme_pub() + theme(panel.grid = element_blank(),
                      axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_plot(p_heat, "g4_rloop_promoter_heatmap", fig_subdir, fig_root,
          width = 4 + 2.2 * length(intersect(heat_genos, genotypes)), height = 7)

# ---------------------------------------------------------------------------
# View 2: stratified metaprofiles (per genotype)
# ---------------------------------------------------------------------------
strata_of <- function(x) {
  br <- stats::quantile(x, probs = seq(0, 1, length.out = n_strata + 1), na.rm = TRUE)
  br <- unique(br); if (length(br) < 3) return(factor(rep(NA, length(x))))
  cut(x, breaks = br, include.lowest = TRUE,
      labels = paste0("Q", seq_len(length(br) - 1)))
}
group_profile <- function(mat, grp) {
  dplyr::bind_rows(lapply(levels(grp), function(q) {
    rows <- which(!is.na(grp) & grp == q)
    if (length(rows) == 0) return(NULL)
    m <- mat[rows, , drop = FALSE]
    data.frame(quartile = q, position = positions,
               mean = colMeans(m, na.rm = TRUE),
               sem  = apply(m, 2, function(z) stats::sd(z, na.rm = TRUE) /
                             sqrt(sum(is.finite(z)))),
               n = length(rows), stringsAsFactors = FALSE)
  }))
}

strat_tbl <- dplyr::bind_rows(lapply(genotypes, function(g) {
  g4m <- get_prof(g, "G4_BG4"); rlm <- get_prof(g, "Rloop_S96")
  s_g4 <- strata_of(rowMeans(g4m, na.rm = TRUE))
  s_rl <- strata_of(rowMeans(rlm, na.rm = TRUE))
  dplyr::bind_rows(
    cbind(genotype = g, stratify_by = "G4_quartile",    signal = rl_lab, group_profile(rlm, s_g4)),
    cbind(genotype = g, stratify_by = "Rloop_quartile", signal = g4_lab, group_profile(g4m, s_rl))
  )
}))
strat_tbl$genotype <- factor(strat_tbl$genotype, levels = geno_lvls)
readr::write_csv(strat_tbl, "results/tables/g4_rloop_stratified_profiles.csv")

plot_strat <- function(df, ylab, ttl, file) {
  p <- ggplot(df, aes(position, mean, colour = quartile, fill = quartile)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ genotype) +
    scale_colour_viridis_d(end = 0.92, name = "quartile") +
    scale_fill_viridis_d(end = 0.92, guide = "none") +
    labs(x = "distance from TSS (bp)", y = ylab, title = ttl) +
    theme_pub()
  save_plot(p, file, fig_subdir, fig_root, width = 9, height = 7)
}
plot_strat(strat_tbl %>% dplyr::filter(stratify_by == "G4_quartile"),
           sprintf("mean %s signal (x background)", rl_lab),
           "R-loop TSS metaprofile by promoter G4 quartile",
           "g4_rloop_stratified_rloopByG4")
plot_strat(strat_tbl %>% dplyr::filter(stratify_by == "Rloop_quartile"),
           sprintf("mean %s signal (x background)", g4_lab),
           "G4 TSS metaprofile by promoter R-loop quartile",
           "g4_rloop_stratified_g4ByRloop")

# ---------------------------------------------------------------------------
# View 2b: G4-at-TSS vs R-loop-at-TSS, one panel per (genotype x quartile)
# ---------------------------------------------------------------------------
# View 2 averages promoters together within a quartile, which is exactly what
# destroys the promoter-to-promoter pairing: two quartiles can have identical
# mean profiles whether or not the individual promoters that are high for G4 are
# the same ones that are high for R-loop. This view keeps the pairing.
#
# "Signal at the TSS" is a deliberately NARROWER window (tss_hw, default 250 bp)
# than the +/- half_width mean that defines the quartiles. If both used the same
# +/- 2 kb mean, x would BE the variable the panels were cut on and every
# panel's r would be attenuated by range restriction on x by construction --
# an artefact shaped exactly like the hypothesis. Using the central bins makes x
# a different measurement. Range restriction is broken, not eliminated (the two
# windows still correlate), so read the r values as a TREND ACROSS PANELS rather
# than as absolute coupling strengths.
#
# Costs no bigWig I/O: these are column subsets of the matrices View 1/2 already
# hold, and the strata come from the same strata_of() call, so panel k here is
# the same promoter set as curve Qk in the metaprofile above.
tss_cols <- which(abs(positions) <= tss_hw)
if (length(tss_cols) < 1)
  stop("tss_scatter_half_width = ", tss_hw, " selects no bin: the profile grid is ",
       signif(positions[2] - positions[1], 3), " bp per bin over +/- ", half_width,
       ". Raise it to at least half a bin width.")
message(sprintf("TSS scatter score: %d of %d bins (|position| <= %d bp, %.1f bp/bin)",
                length(tss_cols), n_bins, tss_hw, positions[2] - positions[1]))

q_lvls <- paste0("Q", seq_len(n_strata))
tss_scores <- dplyr::bind_rows(lapply(genotypes, function(g) {
  g4m <- get_prof(g, "G4_BG4"); rlm <- get_prof(g, "Rloop_S96")
  data.frame(
    genotype       = g,
    entrez         = names(prom_prof),
    g4_tss         = rowMeans(g4m[, tss_cols, drop = FALSE], na.rm = TRUE),
    rloop_tss      = rowMeans(rlm[, tss_cols, drop = FALSE], na.rm = TRUE),
    g4_window      = rowMeans(g4m, na.rm = TRUE),
    rloop_window   = rowMeans(rlm, na.rm = TRUE),
    g4_quartile    = as.character(strata_of(rowMeans(g4m, na.rm = TRUE))),
    rloop_quartile = as.character(strata_of(rowMeans(rlm, na.rm = TRUE))),
    stringsAsFactors = FALSE)
})) %>%
  dplyr::mutate(symbol    = unname(prom_symbol[match(entrez, prom_entrez)]),
                log_g4    = log2(g4_tss    + pc_cnt),
                log_rloop = log2(rloop_tss + pc_cnt)) %>%
  dplyr::filter(is.finite(log_g4), is.finite(log_rloop))
tss_scores$genotype       <- factor(tss_scores$genotype, levels = geno_lvls)
tss_scores$g4_quartile    <- factor(tss_scores$g4_quartile, levels = q_lvls)
tss_scores$rloop_quartile <- factor(tss_scores$rloop_quartile, levels = q_lvls)
readr::write_csv(tss_scores, "results/tables/g4_rloop_tss_scores.csv")

# Correlation within each panel. Spearman as well as Pearson because the TSS
# scores are heavy-tailed and one very high promoter can carry a Pearson r.
tss_corr <- function(qcol) {
  tss_scores %>%
    dplyr::filter(!is.na(.data[[qcol]])) %>%
    dplyr::group_by(genotype, quartile = .data[[qcol]], .drop = FALSE) %>%
    dplyr::summarise(
      n = dplyr::n(),
      pearson_r = if (dplyr::n() >= 3 && stats::sd(log_g4) > 0 && stats::sd(log_rloop) > 0)
                    suppressWarnings(stats::cor(log_g4, log_rloop)) else NA_real_,
      spearman_rho = if (dplyr::n() >= 3)
                    suppressWarnings(stats::cor(log_g4, log_rloop, method = "spearman"))
                    else NA_real_,
      .groups = "drop") %>%
    dplyr::mutate(stratify_by = qcol)
}
tss_stats <- dplyr::bind_rows(tss_corr("g4_quartile"), tss_corr("rloop_quartile"))
readr::write_csv(tss_stats, "results/tables/g4_rloop_tss_quartile_correlation.csv")

# One hex polygon per occupied cell instead of one path per promoter: ~6.1k
# promoters per panel x 16 panels would be ~98k vector objects, which is what
# makes the existing scatter PDFs unopenable in Illustrator.
hex_layer <- function(bins, legend) {
  if (use_hexbin)
    list(geom_hex(bins = bins), scale_fill_viridis_c(trans = "log10", name = legend))
  else geom_point(alpha = 0.2, size = 0.35, colour = "grey25")
}

plot_tss_scatter <- function(qcol, cut_lab, file) {
  df <- tss_scores %>% dplyr::filter(!is.na(.data[[qcol]])) %>%
    dplyr::mutate(quartile = .data[[qcol]])
  st <- tss_stats %>% dplyr::filter(stratify_by == qcol, n > 0)
  st$genotype <- factor(as.character(st$genotype), levels = geno_lvls)
  st$quartile <- factor(as.character(st$quartile), levels = q_lvls)
  p <- ggplot(df, aes(log_g4, log_rloop)) +
  # NOT hex_layer(): that one lets the fill range follow the data, which is right
  # for a pooled scatter and wrong here. Quartering the promoters quarters the
  # count in every hex (measured: median 3-5 per hex, 90th pct 10-19), while the
  # y=0 dropout row -- promoters with no R-loop signal at the TSS -- still piles
  # up to ~2,000 into a handful of hexes. A data-driven ramp is therefore set by
  # that one artefact row and leaves the real cloud in its bottom eighth. Fixed
  # limits put the whole ramp where the structure is and make the 16 panels
  # comparable to each other. oob = squish, so out-of-range hexes are still
  # DRAWN (in the end colours) -- no promoter is dropped from the figure, and no
  # r/rho/n value above is affected.
    (if (use_hexbin)
       list(geom_hex(bins = tss_hex_bins),
            scale_fill_viridis_c(trans = "log10", limits = tss_fill_lim,
                                 breaks = tss_fill_brk, oob = scales::squish,
                                 name = "promoters"))
     else geom_point(alpha = 0.2, size = 0.35, colour = "grey25")) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#E41A1C", linewidth = 0.8) +
    geom_text(data = st, inherit.aes = FALSE, aes(x = -Inf, y = Inf,
              label = sprintf("r = %.2f\nrho = %.2f\nn = %s", pearson_r, spearman_rho,
                              format(n, big.mark = ","))),
              hjust = -0.12, vjust = 1.15, size = 2.4, colour = "grey20", lineheight = 0.95) +
    facet_grid(genotype ~ quartile) +
    labs(x = sprintf("log2 %s signal, TSS +/- %d bp (x background)", g4_lab, tss_hw),
         y = sprintf("log2 %s signal, TSS +/- %d bp (x background)", rl_lab, tss_hw),
         title = sprintf("G4 and R-loop signal at the TSS, within %s quartile", cut_lab),
         subtitle = paste0(
           sprintf("panels cut on mean %s over TSS +/- %d bp; axes are the narrower +/- %d bp score; ",
                   cut_lab, half_width, tss_hw),
           "red = OLS fit\n",
           "r is range-restricted within a panel, so compare r ACROSS panels, ",
           "not against the pooled correlation in rule 25\n",
           sprintf("fill is FIXED at %g-%g promoters per hex and shared by all panels; ",
                   tss_fill_lim[1], tss_fill_lim[2]),
           "hexes below the floor are drawn in the darkest colour, not removed")) +
    theme_pub()
  save_plot(p, file, fig_subdir, fig_root,
            width = 2.9 * n_strata + 2.2, height = 2.6 * length(geno_lvls) + 2)
}
plot_tss_scatter("g4_quartile",    "G4",     "g4_rloop_tss_scatter_byG4quartile")
plot_tss_scatter("rloop_quartile", "R-loop", "g4_rloop_tss_scatter_byRloopquartile")

message("\nWithin-quartile G4 vs R-loop correlation at the TSS (+/- ", tss_hw, " bp):")
print(as.data.frame(tss_stats))

# ===========================================================================
# View 3: G4+/- promoter violin with a random-window null
# ===========================================================================
union_g4 <- readRDS("cache/peaks.rds")$union
GenomeInfoDb::seqlevelsStyle(union_g4) <- "UCSC"
union_rl <- readRDS(file.path(cfg$paths$cache_v2, cfg$rloop$union_rds))
GenomeInfoDb::seqlevelsStyle(union_rl) <- "UCSC"
g4_pos <- IRanges::overlapsAny(prom, union_g4, ignore.strand = TRUE)
rl_pos <- IRanges::overlapsAny(prom, union_rl, ignore.strand = TRUE)
message("Promoters G4-peak+: ", sum(g4_pos), " ; R-loop-peak+: ", sum(rl_pos))

# Depth-normalised mean signal at the random null windows, per (geno, assay).
null_signal <- function(geno, assay) {
  bw <- bw_for(geno, assay)
  m <- vapply(bw$filepath, function(fp) depth_norm_vec(compute_region_signal(bg_gr, fp), fp),
              numeric(length(bg_gr)))
  evict_bigwigs(bw$filepath)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  rowMeans(m, na.rm = TRUE)
}

cliffs_delta <- function(a, b) {  # P(a>b) - P(a<b) via Mann-Whitney U
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 1 || length(b) < 1) return(NA_real_)
  W <- suppressWarnings(stats::wilcox.test(a, b)$statistic)
  unname(2 * W / (length(a) * length(b)) - 1)
}

violin_long <- list(); violin_stats <- list()
for (g in genotypes) {
  g4s <- log2(mean_signal(g, "G4_BG4")   + pc_cnt)
  rls <- log2(mean_signal(g, "Rloop_S96") + pc_cnt)
  g4_null <- log2(null_signal(g, "G4_BG4")   + pc_cnt)
  rl_null <- log2(null_signal(g, "Rloop_S96") + pc_cnt)

  # panel A: R-loop signal stratified by G4 status (+ random null)
  violin_long[[length(violin_long) + 1]] <- dplyr::bind_rows(
    data.frame(genotype = g, panel = "R-loop by G4 status",
               group = ifelse(g4_pos, "G4 peak +", "G4 peak -"), value = rls),
    data.frame(genotype = g, panel = "R-loop by G4 status", group = "random", value = rl_null))
  # panel B: G4 signal stratified by R-loop status (+ random null)
  violin_long[[length(violin_long) + 1]] <- dplyr::bind_rows(
    data.frame(genotype = g, panel = "G4 by R-loop status",
               group = ifelse(rl_pos, "R-loop peak +", "R-loop peak -"), value = g4s),
    data.frame(genotype = g, panel = "G4 by R-loop status", group = "random", value = g4_null))

  violin_stats[[length(violin_stats) + 1]] <- data.frame(
    genotype = g, panel = "R-loop by G4 status",
    median_pos = stats::median(rls[g4_pos], na.rm = TRUE),
    median_neg = stats::median(rls[!g4_pos], na.rm = TRUE),
    median_random = stats::median(rl_null, na.rm = TRUE),
    wilcox_p = suppressWarnings(stats::wilcox.test(rls[g4_pos], rls[!g4_pos])$p.value),
    cliffs_delta = cliffs_delta(rls[g4_pos], rls[!g4_pos]),
    n_pos = sum(g4_pos), n_neg = sum(!g4_pos))
  violin_stats[[length(violin_stats) + 1]] <- data.frame(
    genotype = g, panel = "G4 by R-loop status",
    median_pos = stats::median(g4s[rl_pos], na.rm = TRUE),
    median_neg = stats::median(g4s[!rl_pos], na.rm = TRUE),
    median_random = stats::median(g4_null, na.rm = TRUE),
    wilcox_p = suppressWarnings(stats::wilcox.test(g4s[rl_pos], g4s[!rl_pos])$p.value),
    cliffs_delta = cliffs_delta(g4s[rl_pos], g4s[!rl_pos]),
    n_pos = sum(rl_pos), n_neg = sum(!rl_pos))
}
violin_df    <- dplyr::bind_rows(violin_long)  %>% dplyr::filter(is.finite(value))
violin_stats <- dplyr::bind_rows(violin_stats)
violin_stats$genotype <- factor(violin_stats$genotype, levels = geno_lvls)
readr::write_csv(violin_stats, "results/tables/g4_rloop_violin_stats.csv")

violin_df$genotype <- factor(violin_df$genotype, levels = geno_lvls)
violin_df$group <- factor(violin_df$group, levels = c(
  "G4 peak -", "G4 peak +", "R-loop peak -", "R-loop peak +", "random"))
grp_cols <- c("G4 peak -" = "grey70", "G4 peak +" = "#1B9E77",
              "R-loop peak -" = "grey70", "R-loop peak +" = "#D95F02", "random" = "grey40")
p_violin <- ggplot(violin_df, aes(group, value, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85, colour = "grey30") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.7) +
  facet_grid(panel ~ genotype, scales = "free") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(x = NULL, y = "log2 signal (x background)",
       title = "R-loop signal at G4+/- promoters and vice versa, vs a random-window null",
       subtitle = "Wilcoxon + Cliff's delta in g4_rloop_violin_stats.csv") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_violin, "g4_rloop_promoter_violin", fig_subdir, fig_root, width = 11, height = 7)

# ===========================================================================
# View 4: dG4 vs dR-loop KO views (voom differential at promoters)
# ===========================================================================
assay_differential <- function(assay) {
  cache_file <- file.path("cache", sprintf("cooc_ko_diff_%s.rds", assay))
  cache_or_build(cache_file, {
    bw_all <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = cfg$genotypes, assay = assay)
    out <- dplyr::bind_rows(lapply(kos, function(k) {
      sub <- bw_all[as.character(bw_all$genotype) %in% c(ref_geno, k), , drop = FALSE]
      message("Differential ", assay, " at promoters: ", k, " vs ", ref_geno,
              " (", nrow(sub), " bigWigs)")
      dz  <- peak_voom_lfc(prom, sub, ref_genotype = ref_geno)
      cls <- classify_peak_dz(dz, padj_thresh = voom_padj, lfc_thresh = voom_lfc)
      data.frame(ko = k, entrez = names(prom), logFC = dz$logFC,
                 padj = dz$adj.P.Val, class = as.character(cls), stringsAsFactors = FALSE)
    }))
    evict_bigwigs(bw_all$filepath)
    out
  })
}
g4_dz <- assay_differential("G4_BG4")
rl_dz <- assay_differential("Rloop_S96")

ko_diff <- g4_dz %>%
  dplyr::rename(g4_logFC = logFC, g4_padj = padj, g4_class = class) %>%
  dplyr::inner_join(rl_dz %>% dplyr::rename(rloop_logFC = logFC, rloop_padj = padj,
                                            rloop_class = class),
                    by = c("ko", "entrez")) %>%
  dplyr::mutate(symbol = unname(prom_symbol[match(entrez, prom_entrez)]),
                ko = factor(ko, levels = intersect(geno_lvls, kos)))
readr::write_csv(ko_diff, "results/tables/g4_rloop_ko_differential.csv")

# --- delta-delta scatter (per KO) ---
dd <- ko_diff %>% dplyr::filter(is.finite(g4_logFC), is.finite(rloop_logFC))
dd_stats <- dd %>% dplyr::group_by(ko) %>%
  dplyr::summarise(r = suppressWarnings(stats::cor(g4_logFC, rloop_logFC)),
                   n = dplyr::n(), .groups = "drop")
dd_lab <- setNames(sprintf("%s\nr=%.2f  n=%s", dd_stats$ko, dd_stats$r,
                           format(dd_stats$n, big.mark = ",")), as.character(dd_stats$ko))
p_delta <- ggplot(dd, aes(g4_logFC, rloop_logFC)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  hex_layer(hex_bins, "promoters") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "#E41A1C", linewidth = 0.9) +
  facet_wrap(~ ko, labeller = labeller(ko = dd_lab)) +
  labs(x = sprintf("%s log2FC (KO vs WT)", g4_lab),
       y = sprintf("%s log2FC (KO vs WT)", rl_lab),
       title = "Coordinated change of G4 and R-loop at promoters on helicase loss",
       subtitle = "fill = promoters per hex (log10); r in the strip is over all points, unbinned") +
  theme_pub()
save_plot(p_delta, "g4_rloop_delta_scatter", fig_subdir, fig_root, width = 10, height = 5)

# --- R-loop volcano coloured by G4 differential class (per KO) ---
volc <- ko_diff %>% dplyr::filter(is.finite(rloop_logFC), is.finite(rloop_padj)) %>%
  dplyr::mutate(neglog10p = -log10(pmax(rloop_padj, .Machine$double.xmin)),
                g4_class = factor(ifelse(g4_class %in% c("gained", "lost", "stable"),
                                         g4_class, "ambiguous"),
                                  levels = c("gained", "stable", "lost", "ambiguous")))
# The colour IS the message here, so the whole panel cannot become one hex fill.
# Instead hex only the classes in volcano_hex_classes -- the grey mass, ~70.7k of
# the 73.4k points, which carries no colour information anyway -- and keep
# gained/lost (~2.8k) as real points. Same figure semantics, ~26x fewer vector
# objects. Grey (not viridis) for the hex fill so it stays visually background.
volc_bg <- volc %>% dplyr::filter(as.character(g4_class) %in% volc_hex_cls)
volc_fg <- volc %>% dplyr::filter(!as.character(g4_class) %in% volc_hex_cls)
message("Volcano: ", nrow(volc_bg), " points hex-binned (", paste(volc_hex_cls, collapse = "/"),
        "), ", nrow(volc_fg), " drawn as points")

p_volc <- ggplot(mapping = aes(rloop_logFC, neglog10p)) +
  (if (use_hexbin)
     list(geom_hex(data = volc_bg, bins = hex_bins),
          scale_fill_gradient(low = "grey90", high = "grey45", trans = "log10",
                              name = sprintf("%s
promoters", paste(volc_hex_cls, collapse = "/"))))
   else geom_point(data = volc_bg, aes(colour = g4_class), alpha = 0.5, size = 0.5)) +
  geom_point(data = volc_fg, aes(colour = g4_class), alpha = 0.8, size = 0.7) +
  geom_hline(yintercept = -log10(voom_padj), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = c(-voom_lfc, voom_lfc), linetype = 2, colour = "grey50") +
  facet_wrap(~ ko) +
  scale_colour_manual(values = c(gained = "#E41A1C", stable = "grey75",
                                 lost = "#377EB8", ambiguous = "grey90"),
                      name = "G4 change", drop = FALSE) +
  labs(x = sprintf("%s log2FC (KO vs WT)", rl_lab), y = "-log10 adj.P (R-loop)",
       title = "Differential R-loop at promoters, coloured by G4 differential class",
       subtitle = paste("do R-loop changes concentrate where G4 also changes?",
                        "stable/ambiguous promoters are hex-binned; gained/lost are points")) +
  theme_pub()
save_plot(p_volc, "g4_rloop_rloop_volcano", fig_subdir, fig_root, width = 10, height = 5)

message("Done. Co-occurrence views written. Violin stats:")
print(violin_stats)
message("Delta-delta (G4 vs R-loop log2FC) correlation per KO:")
print(dd_stats)
