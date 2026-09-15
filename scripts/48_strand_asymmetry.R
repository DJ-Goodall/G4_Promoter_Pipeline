#!/usr/bin/env Rscript
# ============================================================================
# 48_strand_asymmetry.R   (system R)   --- Stage L ---
#
# G4 strand-asymmetry analysis. Ported from g4_gloop_strand_analysis.Rmd
# (chunks s1-s8). For each G4 peak in promoter/5'UTR regions: extract flanking
# sequence at each window size, score G4Hunter, assign the peak to the G-rich
# strand, then test genomic-strand bias (binomial), gene-relative bias,
# G4-score distributions (Wilcoxon), detection rate (chi-square), and build
# strand + template-vs-coding metaprofiles.
# gene_name is populated (fixes the blank-label bug used by the Gviz rule 49).
#
# ---------------------------------------------------------------------------
# v1.1 CHANGES (2026-08-07 audit). Three defects, all in this script:
#
# (1) MISSING ignore.strand (the substantive bug). The gene-relative overlap at
#     the old line 135 passed a peak GRanges stamped with `g4_detected_strand`
#     into findOverlaps(), which defaults to ignore.strand = FALSE. A `+` peak
#     could therefore only match a `+` gene, so `gene_strand ==
#     g4_detected_strand` was TRUE by construction: every row of
#     strand_bias_relative.tsv read count_against = 0, pct_with = 100,
#     p ~ 1e-243. Peaks whose nearest gene was on the other strand returned no
#     hit and became NA -- silently DISCARDED rather than counted as "against",
#     which is why n_total there (806) was ~46% of strand_bias_summary (1767).
#     Knock-on: `g4_strand_class` could never take the value "template", so
#     metaprofile_template_vs_coding was structurally a single-curve figure.
#
# (2) WINDOW-MEAN G4Hunter. g4_hunter_score() averaged the run-length score
#     over the whole 501/1001 bp window. A G4 motif occupies ~25 nt, so its
#     signal was diluted ~20-40x: observed scores in v1.0 g4_score_stats.tsv
#     span 0.088-0.153, an order of magnitude below the conventional G4Hunter
#     threshold of 1.2. Setting g4_threshold to 1.2 on the OLD statistic would
#     have returned zero positives in all 32 rows. Fixed by scoring the max
#     over a sliding 25 nt window (g4_hunter_max), which is what the 1.2/1.5
#     literature cutoffs are defined against. The window mean is retained as
#     g4_hunter_mean for continuity with v1.0.
#
# (3) REDUNDANT REVERSE STRAND. Reverse-complementing maps every G-run of
#     length L to a C-run of length L at a mirrored position; rle lengths are
#     preserved, `val` flips sign, and mean() is order-invariant. So
#     g4_hunter_mean(revcomp(s)) == -g4_hunter_mean(s) EXACTLY. `fwd >= rev`
#     was therefore just `fwd >= 0`, and pmax(abs(fwd), abs(rev)) was abs(fwd).
#     The reverse pass is dropped; net-G-richness is now sign(fwd) explicitly.
#
# Also: the header above advertised a Wilcoxon and a chi-square test that were
# never implemented in v1.0. Both are now computed (strand_score_tests.tsv,
# strand_detection_tests.tsv).
# ---------------------------------------------------------------------------
#
# Inputs:   cache/peaks_G4_BG4_{WT,DHX36KO,FANCJKO,dKO}.rds,
#           cache/regions_{promoter,5UTR}.rds, config bigwig_dir
# Outputs:  cache/sequences_scored.rds
#           results/tables/{strand_bias_summary,strand_bias_relative,
#             g4_detection_by_strand,g4_score_stats}.tsv
#           results/tables/{strand_score_tests,strand_detection_tests}.tsv  [v1.1]
#           results/tables/strand_mirrored_profiles.csv                   [v1.2]
#           results/figures/22_strand/{strand_bias_comprehensive,
#             metaprofile_strand_mirrored,metaprofile_strand_mirrored_generel,
#             mean_signal_by_strand}.{pdf,png}
#
# v1.2 (2026-09-15) RETIRES metaprofile_by_strand.* and
# metaprofile_template_vs_coding.*: both were start-anchored with a
# 3-replicate SEM. See the block above the new figures for what replaced them.
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(BSgenome.Mmusculus.UCSC.mm10); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db); library(AnnotationDbi); library(rtracklayer)
  library(IRanges); library(ggplot2); library(cowplot); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures", "results/tables")
fig_root <- "results/figures"
sc <- cfg$strand
window_sizes   <- as.integer(unlist(sc$window_sizes))
g4_threshold   <- sc$g4_threshold   # NULL -> calibrate per window size (v1.1)
genotype_order <- c("WT", "DHX36KO", "FANCJKO", "dKO")
region_types   <- c("promoter", "5'UTR")
n_bins <- sc$n_bins; half_width <- sc$half_width

genome <- BSgenome.Mmusculus.UCSC.mm10
txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene

# --- G4Hunter ---------------------------------------------------------------
# Per-base run-length scoring, shared by both statistics below: a base inside a
# G-run of length L scores +min(L, 4), inside a C-run -min(L, 4), else 0.
g4_hunter_bases <- function(seq_str) {
  if (is.na(seq_str) || nchar(seq_str) == 0) return(NULL)
  chars <- strsplit(toupper(seq_str), "")[[1]]
  r <- rle(chars)
  val <- ifelse(r$values == "G",  pmin(r$lengths, 4L),
         ifelse(r$values == "C", -pmin(r$lengths, 4L), 0L))
  rep(val, r$lengths)
}

# v1.0 statistic: mean over the WHOLE window. Retained for continuity, but it
# is diluted by the flank and is NOT comparable to published G4Hunter cutoffs.
g4_hunter_mean <- function(seq_str) {
  s <- g4_hunter_bases(seq_str)
  if (is.null(s)) return(NA_real_)
  m <- mean(s, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

# v1.1 statistic: max over a sliding `win` nt window (signed, so the sign still
# tells you which strand is G-rich). This is the quantity the conventional
# G4Hunter thresholds (1.2 permissive / 1.5 stringent) are defined against --
# Bedrat et al. 2016 NAR score 25-nt windows, not whole loci. Windows shorter
# than `win` fall back to the plain mean so short inputs still score.
g4_hunter_max <- function(seq_str, win = 25L) {
  s <- g4_hunter_bases(seq_str)
  if (is.null(s)) return(NA_real_)
  n <- length(s)
  if (n < win) {
    m <- mean(s, na.rm = TRUE)
    return(if (is.nan(m)) NA_real_ else m)
  }
  # Rolling mean via cumulative sum: O(n) rather than O(n * win).
  cs   <- c(0, cumsum(s))
  roll <- (cs[(win + 1L):(n + 1L)] - cs[1L:(n - win + 1L)]) / win
  if (!length(roll) || all(is.na(roll))) return(NA_real_)
  # Signed extreme: keep whichever of the G-rich / C-rich peaks is larger in
  # magnitude, so the sign remains an honest "which strand is G-rich" call.
  hi <- max(roll, na.rm = TRUE); lo <- min(roll, na.rm = TRUE)
  if (abs(hi) >= abs(lo)) hi else lo
}

# Back-compat alias: anything that previously called g4_hunter_score() now gets
# the sliding-window statistic.
g4_hunter_score <- g4_hunter_max

# --- Load peaks + assign region_type ----------------------------------------
g4_peaks_list <- lapply(genotype_order, function(g) {
  f <- file.path("cache", sprintf("peaks_G4_BG4_%s.rds", g))
  pk <- readRDS(f); pk$genotype <- g; pk$assay <- "G4_BG4"; pk })
g4_peaks <- Reduce(c, g4_peaks_list)

reg <- list(promoter = readRDS("cache/regions_promoter.rds"),
            `5'UTR`   = readRDS("cache/regions_5UTR.rds"))
g4_peaks$region_type <- NA_character_
for (rt in region_types) {
  ov <- GenomicRanges::findOverlaps(g4_peaks, reg[[rt]])
  idx <- unique(S4Vectors::queryHits(ov))
  idx <- idx[is.na(g4_peaks$region_type[idx])]
  g4_peaks$region_type[idx] <- rt
}
g4_peaks <- g4_peaks[!is.na(g4_peaks$region_type)]

# subsample 2000/genotype
set.seed(sc$seed)
g4_peaks <- do.call(c, lapply(genotype_order, function(g) {
  pg <- g4_peaks[g4_peaks$genotype == g]
  if (length(pg) > sc$subsample) pg[sample(length(pg), sc$subsample)] else pg }))

# gene_name per peak (bug fix for Gviz): overlapping gene -> SYMBOL
genes_gr <- GenomicFeatures::genes(txdb)
ov <- GenomicRanges::findOverlaps(g4_peaks, genes_gr, select = "first")
entrez <- names(genes_gr)[ov]
# NA keys (peaks not overlapping any gene) make mapIds() return a *list*, which
# then explodes into per-element columns inside data.frame() below (-> "differing
# number of rows: 1, 0"). Map only the non-NA keys and keep a plain char vector.
sym <- rep(NA_character_, length(g4_peaks))
valid <- !is.na(entrez)
if (any(valid)) {
  m <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = entrez[valid],
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  sym[valid] <- unname(as.character(m))
}
g4_peaks$gene_name <- sym

# --- Extract sequences (per region_type x window) ---------------------------
scored_cache <- "cache/sequences_scored.rds"
sequences_df <- cache_or_build(scored_cache, {
  extract_seqs <- function(gr, window_bp) {
    ex <- GenomicRanges::trim(GenomicRanges::resize(gr, width = 2 * window_bp + 1, fix = "center"))
    seqs <- BSgenome::getSeq(genome, ex)
    data.frame(seqname = as.character(GenomicRanges::seqnames(gr)),
      peak_start = GenomicRanges::start(gr), peak_end = GenomicRanges::end(gr),
      strand = as.character(GenomicRanges::strand(gr)), genotype = gr$genotype,
      gene_name = gr$gene_name, region_type = gr$region_type, window_size = window_bp,
      sequence = as.character(seqs), stringsAsFactors = FALSE)
  }
  parts <- list()
  for (rt in region_types) for (w in window_sizes) {
    pk <- g4_peaks[g4_peaks$region_type == rt]
    if (length(pk) == 0) next
    parts[[paste(rt, w)]] <- extract_seqs(pk, w)
  }
  sdf <- do.call(rbind, parts); rownames(sdf) <- NULL

  # G4Hunter scoring. v1.0 scored both strands; the reverse pass is provably
  # redundant (g4_hunter(revcomp(s)) == -g4_hunter(s) exactly -- see header note
  # 3), so we score the + strand once and read the - strand off the sign.
  fwd      <- vapply(sdf$sequence, g4_hunter_max,  numeric(1), USE.NAMES = FALSE)
  fwd_mean <- vapply(sdf$sequence, g4_hunter_mean, numeric(1), USE.NAMES = FALSE)

  # Net-G-richness: fwd > 0 -> reference + strand is the G-rich one.
  det <- ifelse(is.na(fwd), "neither", ifelse(fwd >= 0, "+", "-"))

  sdf$g4_score_fwd  <- fwd
  sdf$g4_score_rev  <- -fwd            # exact, by the identity above
  sdf$g4_score_mean <- fwd_mean        # v1.0 window-mean statistic, for comparison
  sdf$g4_detected_strand <- det
  sdf$g4_score <- abs(fwd)
  # NB has_g4_motif is deliberately NOT computed here: it depends on a threshold
  # that is calibrated below, and we do not want a threshold change to
  # invalidate the (expensive) sequence-extraction cache.
  sdf
})

# --- v1.1: calibrate the motif-call threshold against real mm10 -------------
# The audit asked for `g4_threshold` to be set to a "real G4Hunter cutoff
# (1.2 conventional, 1.5 stringent)". Measured against 5,000 random mm10
# windows (verify/calibrate_g4_threshold.R), neither literature value works
# as a per-window classifier of the sliding-25nt max:
#     501 bp windows:  FPR at 1.2 = 0.441,  at 1.5 = 0.230
#    1001 bp windows:  FPR at 1.2 = 0.635,  at 1.5 = 0.396
# because a max over ~477 overlapping 25 nt windows is an extreme-value
# statistic -- the published cutoffs are genome-SCANNING thresholds for
# locating candidate sites, not per-locus decision thresholds. The FPR also
# rises with window width, so a single fixed number cannot serve both sizes.
# We therefore calibrate PER WINDOW SIZE to a stated false-positive rate and
# write the provenance out. Set strand.g4_threshold to a number to override.
calib_thresholds <- local({
  if (!is.null(g4_threshold) && is.numeric(g4_threshold) && g4_threshold > 0) {
    message("Using fixed strand.g4_threshold = ", g4_threshold,
            " (calibration skipped; FPR unknown)")
    return(setNames(rep(as.numeric(g4_threshold), length(window_sizes)),
                    as.character(window_sizes)))
  }
  n_bg   <- if (!is.null(sc$calib_n)) as.integer(sc$calib_n) else 5000L
  fpr    <- if (!is.null(sc$g4_threshold_fpr)) as.numeric(sc$g4_threshold_fpr) else 0.05
  cseed  <- if (!is.null(sc$calib_seed)) as.integer(sc$calib_seed) else 2026L
  sizes  <- GenomeInfoDb::seqlengths(genome)[cfg$std_chroms]
  set.seed(cseed)
  rows <- list(); thr <- numeric(0)
  for (w in window_sizes) {
    width_bp <- 2L * w + 1L
    chr <- sample(names(sizes), n_bg, replace = TRUE, prob = as.numeric(sizes))
    st  <- floor(runif(n_bg) * (as.numeric(sizes[chr]) - width_bp - 1)) + 1
    bg  <- GenomicRanges::GRanges(chr, IRanges::IRanges(as.integer(st), width = width_bp))
    sq  <- BSgenome::getSeq(genome, bg)
    sq  <- sq[as.numeric(Biostrings::letterFrequency(sq, "N")) == 0]
    sc_bg <- abs(vapply(as.character(sq), g4_hunter_max, numeric(1), USE.NAMES = FALSE))
    t_w <- unname(quantile(sc_bg, 1 - fpr, na.rm = TRUE))
    thr[as.character(w)] <- t_w
    rows[[length(rows) + 1]] <- data.frame(
      window_size = w, extracted_width_bp = width_bp, n_background = length(sq),
      target_fpr = fpr, calibrated_threshold = t_w,
      fpr_at_literature_1.2 = mean(sc_bg >= 1.2),
      fpr_at_literature_1.5 = mean(sc_bg >= 1.5))
    message(sprintf("  window %d bp (%d bp extracted): threshold %.3f at FPR %.2f  ",
                    w, width_bp, t_w, fpr),
            sprintf("[literature 1.2 would give FPR %.3f]", mean(sc_bg >= 1.2)))
  }
  readr::write_tsv(dplyr::bind_rows(rows),
                   "results/tables/strand_g4_threshold_calibration.tsv")
  thr
})

sequences_df$g4_threshold_used <- calib_thresholds[as.character(sequences_df$window_size)]
sequences_df$has_g4_motif <- !is.na(sequences_df$g4_score) &
  sequences_df$g4_score >= sequences_df$g4_threshold_used

# --- Gene-relative strand + template/coding ---------------------------------
# v1.1 BUG FIX. gr_peaks carries the G4-detected strand because downstream
# metaprofiles need it for orientation -- but that strand must NOT constrain
# which gene a peak is assigned to. findOverlaps() defaults to
# ignore.strand = FALSE, so in v1.0 a `+` peak could only ever match a `+`
# gene and `gene_strand == g4_detected_strand` was true by construction.
# ignore.strand = TRUE makes the coding-vs-template comparison meaningful.
gr_peaks <- GenomicRanges::GRanges(sequences_df$seqname,
  IRanges::IRanges(sequences_df$peak_start, sequences_df$peak_end),
  strand = ifelse(sequences_df$g4_detected_strand %in% c("+", "-"),
                  sequences_df$g4_detected_strand, "*"))
gs_ov <- GenomicRanges::findOverlaps(gr_peaks, genes_gr, select = "first",
                                     ignore.strand = TRUE)
sequences_df$gene_strand <- as.character(GenomicRanges::strand(genes_gr))[gs_ov]
sequences_df$strand_class <- with(sequences_df, ifelse(
  is.na(gene_strand) | g4_detected_strand == "neither", NA_character_,
  ifelse(g4_detected_strand == gene_strand, "with_strand", "against_strand")))
sequences_df$g4_strand_class <- with(sequences_df, ifelse(
  is.na(gene_strand) | g4_detected_strand == "neither", "unmapped",
  ifelse(g4_detected_strand == gene_strand, "coding", "template")))

# --- Statistics tables ------------------------------------------------------
grid <- expand.grid(genotype = genotype_order, region_type = region_types,
                    window_size = window_sizes, stringsAsFactors = FALSE)
strand_bias_summary <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i], ]
  cp <- sum(s$g4_detected_strand == "+"); cm <- sum(s$g4_detected_strand == "-")
  nt <- cp + cm
  data.frame(grid[i, ], count_plus = cp, count_minus = cm, n_total = nt,
    pct_plus = if (nt > 0) 100 * cp / nt else NA,
    p_value_binom = if (nt > 0) binom.test(cp, nt, 0.5)$p.value else NA,
    fold_plus = if (nt > 0) (cp / nt) / 0.5 else NA)
}))
strand_bias_summary$padj <- p.adjust(strand_bias_summary$p_value_binom, method = "BH")

strand_bias_relative <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i] & !is.na(sequences_df$strand_class), ]
  cw <- sum(s$strand_class == "with_strand"); ca <- sum(s$strand_class == "against_strand")
  nt <- cw + ca
  data.frame(grid[i, ], count_with = cw, count_against = ca, n_total = nt,
    pct_with = if (nt > 0) 100 * cw / nt else NA,
    p_value_binom = if (nt > 0) binom.test(cw, nt, 0.5)$p.value else NA)
}))
strand_bias_relative$padj <- p.adjust(strand_bias_relative$p_value_binom, method = "BH")

g4_score_stats <- sequences_df %>%
  dplyr::filter(g4_detected_strand %in% c("+", "-"), !is.na(g4_score)) %>%
  dplyr::group_by(genotype, region_type, window_size, strand = g4_detected_strand) %>%
  dplyr::summarise(n = dplyr::n(), mean_g4 = mean(g4_score), median_g4 = median(g4_score),
    sd_g4 = sd(g4_score), q25 = quantile(g4_score, .25), q75 = quantile(g4_score, .75),
    .groups = "drop")

g4_detection_by_strand <- sequences_df %>%
  dplyr::filter(g4_detected_strand %in% c("+", "-")) %>%
  dplyr::group_by(genotype, region_type, window_size, strand = g4_detected_strand) %>%
  dplyr::summarise(n_total = dplyr::n(), n_with_motif = sum(has_g4_motif, na.rm = TRUE),
    pct_with_motif = 100 * mean(has_g4_motif, na.rm = TRUE), .groups = "drop")

# --- v1.1: the two tests the v1.0 header advertised but never computed ------
# (a) Wilcoxon: does the G4Hunter score differ between + and - detected strands?
strand_score_tests <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i] &
    sequences_df$g4_detected_strand %in% c("+", "-") & !is.na(sequences_df$g4_score), ]
  xp <- s$g4_score[s$g4_detected_strand == "+"]
  xm <- s$g4_score[s$g4_detected_strand == "-"]
  data.frame(grid[i, ], n_plus = length(xp), n_minus = length(xm),
    median_plus  = if (length(xp)) median(xp) else NA_real_,
    median_minus = if (length(xm)) median(xm) else NA_real_,
    p_value_wilcox = if (length(xp) > 1 && length(xm) > 1)
      suppressWarnings(wilcox.test(xp, xm)$p.value) else NA_real_)
}))
strand_score_tests$padj <- p.adjust(strand_score_tests$p_value_wilcox, method = "BH")

# (b) Chi-square: does the motif-detection RATE depend on detected strand?
#     Only meaningful now that g4_threshold > 0 -- in v1.0 every cell was 100%.
strand_detection_tests <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i] &
    sequences_df$g4_detected_strand %in% c("+", "-"), ]
  tab <- table(factor(s$g4_detected_strand, c("+", "-")),
               factor(s$has_g4_motif, c(TRUE, FALSE)))
  ok <- all(dim(tab) == c(2, 2)) && sum(tab) > 0 && all(rowSums(tab) > 0) &&
        all(colSums(tab) > 0)
  data.frame(grid[i, ],
    n_plus_motif  = tab["+", "TRUE"],  n_plus_nomotif  = tab["+", "FALSE"],
    n_minus_motif = tab["-", "TRUE"],  n_minus_nomotif = tab["-", "FALSE"],
    p_value_chisq = if (ok) suppressWarnings(chisq.test(tab)$p.value) else NA_real_)
}))
strand_detection_tests$padj <- p.adjust(strand_detection_tests$p_value_chisq, method = "BH")

readr::write_tsv(strand_bias_summary,     "results/tables/strand_bias_summary.tsv")
readr::write_tsv(strand_bias_relative,    "results/tables/strand_bias_relative.tsv")
readr::write_tsv(g4_detection_by_strand,  "results/tables/g4_detection_by_strand.tsv")
readr::write_tsv(g4_score_stats,          "results/tables/g4_score_stats.tsv")
readr::write_tsv(strand_score_tests,      "results/tables/strand_score_tests.tsv")
readr::write_tsv(strand_detection_tests,  "results/tables/strand_detection_tests.tsv")

# --- v1.1 sanity guards -----------------------------------------------------
# These are the exact pathologies the audit found in v1.0. Warn loudly rather
# than silently emitting a table that cannot say anything.
n_unmapped <- sum(is.na(sequences_df$gene_strand))
message(sprintf("Peaks with no overlapping gene: %d / %d (%.1f%%) -- excluded from
  the gene-relative test. In v1.0 this also swallowed every against-strand peak.",
  n_unmapped, nrow(sequences_df), 100 * n_unmapped / nrow(sequences_df)))
if (all(strand_bias_relative$count_against == 0, na.rm = TRUE))
  warning("strand_bias_relative: count_against is 0 in EVERY row. The ",
          "ignore.strand fix has not taken effect -- check findOverlaps above.")
if (all(g4_detection_by_strand$pct_with_motif == 100, na.rm = TRUE))
  warning("g4_detection_by_strand: 100% in every row. g4_threshold (", g4_threshold,
          ") is not discriminating -- it must be > 0 and applied to the ",
          "sliding-window score, not the window mean.")

# --- Comprehensive strand-bias figure (multi-panel) -------------------------
gl <- function(df) { df$genotype <- factor(df$genotype, levels = genotype_order); df }
pA <- ggplot(gl(strand_bias_summary), aes(genotype, pct_plus, fill = genotype)) +
  geom_col() + geom_hline(yintercept = 50, linetype = "dashed") +
  geom_text(aes(label = ifelse(padj < 0.05, "*", "")), vjust = -0.2) +
  facet_grid(region_type ~ window_size) + scale_fill_manual(values = condition_colours) +
  labs(title = "A. Genomic-strand G4 bias (% + strand)", y = "% + strand", x = NULL) +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
pB <- ggplot(gl(strand_bias_relative), aes(genotype, pct_with, fill = genotype)) +
  geom_col() + geom_hline(yintercept = 50, linetype = "dashed") +
  facet_grid(region_type ~ window_size) + scale_fill_manual(values = condition_colours) +
  labs(title = "B. Gene-relative bias (% with gene strand)", y = "% with-strand", x = NULL) +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
pC <- ggplot(gl(dplyr::filter(sequences_df, g4_detected_strand %in% c("+", "-"), window_size == 250)),
             aes(genotype, g4_score, fill = g4_detected_strand)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + facet_wrap(~ region_type) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "C. G4Hunter score by strand (250 bp)", y = "G4Hunter score",
       x = NULL, fill = "strand") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pD <- ggplot(gl(g4_detection_by_strand), aes(genotype, pct_with_motif, fill = strand)) +
  geom_col(position = "dodge") + facet_grid(region_type ~ window_size) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "D. G4 detection rate by strand", y = "% with motif", x = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_comp <- cowplot::plot_grid(pA, pB, pC, pD, ncol = 2, labels = NULL)
save_plot(p_comp, "strand_bias_comprehensive", "22_strand", fig_root, width = 16, height = 12)

# --- v1.2: mirrored, region-centred strand metaprofiles ---------------------
# REPLACES metaprofile_by_strand / mean_signal_by_strand /
# metaprofile_template_vs_coding. Three changes, in order of how much they
# mattered:
#
# (1) ANCHOR. The old prof_matrix() built its window from
#     start(region) - half_width, so peaks were aligned by their LEFT EDGE and
#     the centre of the plot was not the G4. Peak widths here span a wide range,
#     so that smears exactly the feature the figure is about. The replacement
#     anchors on the region centre through center_window() + the shared
#     compute_profile_matrix(), which is what rules 24/30/31 use -- so these
#     curves are now on the same grid as every other metaprofile in the repo.
#
# (2) UNIT. The old figure profiled PEAKS. This profiles the promoter / 5'UTR
#     REGION that hosts the peak, anchored at its centre (= the TSS for
#     promoters). A peak-anchored curve necessarily peaks at zero -- that is
#     where the peak is -- so it cannot show WHERE relative to the TSS a
#     strand-assigned G4 sits, which is the question.
#
# (3) SEM. The old prof_ci() took sd across the 3 REPLICATES (n = 3), i.e. a
#     statement about library reproducibility, not about the regions. This uses
#     profile_from_rep_mean(), whose SEM is across regions, matching every other
#     metaprofile here. See the limitation note printed at the end.
#
# LAYOUT. One curve is drawn above y = 0 and the other below it, both as
# POSITIVE intensities -- the y-axis is labelled with absolute values, so a peak
# below the line is a peak, not a depletion. Mirroring costs nothing
# statistically; it removes the overplotting that made the two curves
# indistinguishable where they agreed.
#
# TWO SPLITS, and they are not interchangeable:
#   genomic  -- G4 on the reference + vs - strand. This is the literal +/-
#               question, but note the caveat in the figure subtitle: the window
#               is GENE-ORIENTED (compute_profile_matrix flips minus-strand
#               regions, helpers.R:314), so a "+ strand G4" curve pools
#               coding-strand G4s at + genes with template-strand G4s at -
#               genes. It answers a genome-level question on a gene-level axis.
#   generel  -- coding (non-template) vs template strand, i.e. the same call
#               made relative to the host gene. On a gene-oriented axis this is
#               the coherent split, and it is the one to read for biology.
# Both are row subsets of the same matrix, so the second is free.
mw <- as.integer(sc$mirror_window_size %||% 250L)
if (!mw %in% window_sizes)
  stop("strand.mirror_window_size (", mw, ") is not one of strand.window_sizes (",
       paste(window_sizes, collapse = ", "), ").")
min_n <- as.integer(sc$mirror_min_regions %||% 20L)

bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotype_order, assay = "G4_BG4")
bw_by_genotype <- split(bw_meta$filepath, as.character(bw_meta$genotype))

chrom_sizes <- default_chrom_sizes(cfg$std_chroms)

# Region centres, one row per (genotype x hosted region). A promoter may host
# several scored peaks; the STRONGEST (max |G4Hunter|) supplies the call, which
# is the same "strongest overlapping peak wins" convention assign_region_topology()
# uses for topology -- not a new rule invented here.
build_region_set <- function(g) {
  sdf <- sequences_df[sequences_df$genotype == g & sequences_df$window_size == mw &
                      sequences_df$g4_detected_strand %in% c("+", "-"), ]
  if (!nrow(sdf)) return(NULL)
  pk <- GenomicRanges::GRanges(sdf$seqname,
          IRanges::IRanges(sdf$peak_start, sdf$peak_end))
  out <- list()
  for (rt in region_types) {
    rr <- reg[[rt]]
    ov <- GenomicRanges::findOverlaps(pk, rr, ignore.strand = TRUE)
    if (!length(ov)) next
    q <- S4Vectors::queryHits(ov); sj <- S4Vectors::subjectHits(ov)
    # strongest peak per region
    o <- order(sj, -abs(sdf$g4_score[q]))
    q <- q[o]; sj <- sj[o]
    keep <- !duplicated(sj)
    q <- q[keep]; sj <- sj[keep]
    win_rt <- center_window(rr[sj], width = 2L * half_width)
    S4Vectors::mcols(win_rt) <- NULL
    win_rt$region_type <- rt
    win_rt$genomic     <- sdf$g4_detected_strand[q]
    win_rt$generel     <- sdf$g4_strand_class[q]      # coding / template / unmapped
    out[[rt]] <- win_rt
  }
  if (!length(out)) return(NULL)
  w <- suppressWarnings(do.call(c, unname(out)))
  GenomeInfoDb::seqlevels(w, pruning.mode = "coarse") <-
    intersect(GenomeInfoDb::seqlevels(w), names(chrom_sizes))
  GenomeInfoDb::seqlengths(w) <- chrom_sizes[GenomeInfoDb::seqlevels(w)]
  w <- GenomicRanges::trim(w)
  # rep_mean_matrix() refuses unequal widths (compute_profile_matrix() recycles
  # them, helpers.R:303), so drop anything the chromosome edge clamped.
  n_pre <- length(w)
  w <- w[GenomicRanges::width(w) == 2L * half_width]
  if (length(w) < n_pre)
    message("  ", g, ": dropped ", n_pre - length(w), " edge-clamped window(s).")
  w
}

positions <- seq(-half_width, half_width, length.out = n_bins)
strand_levels <- list(
  genomic = c("+", "-"),
  generel = c("coding", "template"))
strand_cols <- c("+" = "#E41A1C", "-" = "#377EB8",
                 coding = "#D95F02", template = "#1B9E77")

prof_rows <- list()
for (g in genotype_order) {
  w <- build_region_set(g)
  if (is.null(w) || !length(w)) { message("No scored regions for ", g); next }
  message("Profiling ", g, ": ", length(w), " regions over ",
          length(bw_by_genotype[[g]]), " bigWigs ...")
  rm_g <- rep_mean_matrix(w, bw_by_genotype[[g]], n_bins = n_bins,
                          half_width = half_width)
  for (split in names(strand_levels)) {
    for (rt in region_types) {
      for (lv in strand_levels[[split]]) {
        idx <- which(w$region_type == rt & S4Vectors::mcols(w)[[split]] == lv)
        # n is recorded BEFORE the minimum-size guard, so a curve that is too
        # small to draw still reports its true size rather than vanishing.
        n_lv <- length(idx)
        if (n_lv == 0L) {
          prof_rows[[length(prof_rows) + 1L]] <- data.frame(
            bin = NA_integer_, position = NA_real_, mean = NA_real_, sem = NA_real_,
            n_regions = 0L, genotype = g, region_type = rt, split = split,
            level = lv, drawn = FALSE, stringsAsFactors = FALSE)
          next
        }
        p <- profile_from_rep_mean(rm_g, idx, positions = positions)
        p$genotype <- g; p$region_type <- rt; p$split <- split; p$level <- lv
        p$drawn <- n_lv >= min_n
        prof_rows[[length(prof_rows) + 1L]] <- p
      }
    }
  }
  evict_bigwigs(bw_by_genotype[[g]])   # 12 whole-genome RleLists will not fit
}
if (!length(prof_rows)) stop("No strand metaprofiles could be computed.")

strand_prof <- dplyr::bind_rows(prof_rows)
strand_prof$genotype <- factor(strand_prof$genotype, levels = genotype_order)
# The mirror is a PLOTTING transform only: `mean` stays the true intensity and
# `signed` is what is drawn. Anything downstream must read `mean`.
down_level <- c(genomic = "-", generel = "template")
strand_prof$signed <- ifelse(strand_prof$level == down_level[strand_prof$split],
                             -strand_prof$mean, strand_prof$mean)
strand_prof$signed_lo <- strand_prof$signed -
  ifelse(strand_prof$level == down_level[strand_prof$split], -1, 1) * strand_prof$sem
strand_prof$signed_hi <- strand_prof$signed +
  ifelse(strand_prof$level == down_level[strand_prof$split], -1, 1) * strand_prof$sem
readr::write_csv(strand_prof, "results/tables/strand_mirrored_profiles.csv")

sem_note <- paste("ribbon = SEM across regions, NOT a test of the strand",
                  "difference; free y per panel")

mirror_fig <- function(split, fname, title, up_lab, dn_lab, caveat = NULL) {
  df <- strand_prof %>% dplyr::filter(split == !!split, drawn, !is.na(bin))
  if (!nrow(df)) { message("Nothing to draw for split '", split, "'."); return(invisible(NULL)) }
  lv <- strand_levels[[split]]
  df$level <- factor(df$level, levels = lv)
  n_lab <- strand_prof %>%
    dplyr::filter(split == !!split) %>%
    dplyr::distinct(genotype, region_type, level, n_regions, drawn) %>%
    dplyr::group_by(genotype, region_type) %>%
    dplyr::summarise(label = paste(sprintf("%s n=%s%s", level,
                        format(n_regions, big.mark = ",", trim = TRUE),
                        ifelse(drawn, "", " (not drawn)")), collapse = "  |  "),
                     .groups = "drop")
  p <- ggplot(df, aes(position, signed, colour = level, fill = level)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = signed_lo, ymax = signed_hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.03, vjust = 1.4, size = 2.3, colour = "grey25") +
    facet_grid(region_type ~ genotype, scales = "free_y") +
    # Both halves are positive intensities: the sign is layout, not direction.
    scale_y_continuous(labels = function(x) format(abs(x), trim = TRUE)) +
    scale_colour_manual(values = strand_cols, name = NULL, drop = FALSE) +
    scale_fill_manual(values = strand_cols, guide = "none", drop = FALSE) +
    labs(title = title,
         subtitle = paste0(sprintf("%s ABOVE the line, %s BELOW; both axes are positive intensity (%d bp window, %d bins). ",
                                   up_lab, dn_lab, 2L * half_width, n_bins),
                           if (!is.null(caveat)) paste0("\n", caveat, " ") else "",
                           "\n", sem_note),
         x = "Distance to region centre (bp; = TSS for promoters)",
         y = "Mean G4 CUT&Tag signal") +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, "22_strand", fig_root, width = 14, height = 8)
}

mirror_fig("genomic", "metaprofile_strand_mirrored",
  "G4 CUT&Tag signal by the strand carrying the G4 (genomic)",
  "+ strand G4", "- strand G4",
  paste("CAVEAT: the window is gene-oriented (minus-strand regions are flipped),",
        "so the '+' curve pools coding-strand G4s at + genes with template-strand",
        "G4s at - genes. For the gene-relative split see the companion figure."))

mirror_fig("generel", "metaprofile_strand_mirrored_generel",
  "G4 CUT&Tag signal by the strand carrying the G4 (relative to the host gene)",
  "coding (non-template) strand G4", "template strand G4",
  paste("This is the coherent split on a gene-oriented axis, and the one to read",
        "for biology. Peaks with no overlapping gene are 'unmapped' and excluded."))

# Mean signal per strand, recomputed from the SAME profiles as the figures above
# so the bar chart and the curves can no longer disagree.
p_meanbar <- strand_prof %>%
  dplyr::filter(!is.na(bin), n_regions > 0) %>%
  dplyr::group_by(genotype, region_type, split, level) %>%
  dplyr::summarise(mean_signal = mean(mean, na.rm = TRUE),
                   n_regions = dplyr::first(n_regions), .groups = "drop") %>%
  ggplot(aes(level, mean_signal, fill = level)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_regions, big.mark = ",", trim = TRUE)),
            vjust = -0.3, size = 2.4, colour = "grey25") +
  facet_grid(region_type ~ genotype + split, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = strand_cols, guide = "none") +
  labs(title = "Mean G4 signal by the strand carrying the G4",
       subtitle = "window mean of the same profiles plotted above; bar labels are n regions",
       x = NULL, y = "Mean signal") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_meanbar, "mean_signal_by_strand", "22_strand", fig_root,
          width = 15, height = 8)

message("\nRegions per (genotype, region_type, split, level):")
print(as.data.frame(strand_prof %>%
  dplyr::distinct(genotype, region_type, split, level, n_regions, drawn)))
message("\nNOTE: SEM is across regions, so at these n the ribbon is narrow and every ",
        "visible gap looks decisive. It measures region-to-region spread and is not ",
        "a test of the strand difference.")
message("Strand asymmetry analysis done.")
