#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# v1.1: calibrate the G4Hunter motif-call threshold against REAL mm10 sequence.
#
# Why this exists. The audit asked for `g4_threshold` to be set to "a real
# G4Hunter cutoff (1.2 conventional, 1.5 stringent)". Two problems:
#   - Applied to the v1.0 WINDOW-MEAN statistic, 1.2 calls zero of the ten
#     literature control G4s (verified in verify/test_g4hunter.R).
#   - Applied to the v1.1 SLIDING-MAX statistic over a 501 bp window, 1.2 is
#     far too permissive, because a max over ~477 overlapping windows is an
#     extreme-value statistic. The literature 1.2 is a genome-SCANNING cutoff
#     for locating candidate sites, not a per-window classification cutoff.
# So the threshold is calibrated here: draw random mm10 windows of the same
# width, score them the same way, and take the quantile giving the target
# false-positive rate. Reported, not asserted.
#
# Run from the repository root:  Rscript verify/calibrate_g4_threshold.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(BSgenome.Mmusculus.UCSC.mm10); library(GenomicRanges); library(Biostrings)
})

src <- readLines("scripts/48_strand_asymmetry.R")
s <- grep("^g4_hunter_bases <- function", src)
e <- grep("^g4_hunter_score <- g4_hunter_max", src)
eval(parse(text = paste(src[s:e], collapse = "\n")))

genome <- BSgenome.Mmusculus.UCSC.mm10
std <- paste0("chr", c(1:19, "X"))
sizes <- GenomeInfoDb::seqlengths(genome)[std]

n_bg <- 5000
set.seed(2026)

calib <- function(width_bp, target_fpr = 0.05) {
  chr <- sample(std, n_bg, replace = TRUE, prob = as.numeric(sizes))
  st  <- floor(runif(n_bg) * (as.numeric(sizes[chr]) - width_bp - 1)) + 1
  gr  <- GRanges(chr, IRanges(as.integer(st), width = width_bp))
  sq  <- BSgenome::getSeq(genome, gr)
  keep <- as.numeric(Biostrings::letterFrequency(sq, "N")) == 0
  sq <- sq[keep]
  ss <- as.character(sq)
  mx <- vapply(ss, g4_hunter_max,  numeric(1), USE.NAMES = FALSE)
  mn <- vapply(ss, g4_hunter_mean, numeric(1), USE.NAMES = FALSE)
  gc <- as.numeric(Biostrings::letterFrequency(sq, "GC", as.prob = TRUE))
  list(width = width_bp, n = length(ss), gc = mean(gc),
       thr = unname(quantile(abs(mx), 1 - target_fpr)),
       fpr_at_1.2 = mean(abs(mx) >= 1.2),
       fpr_at_1.5 = mean(abs(mx) >= 1.5),
       mean_fpr_at_1.2 = mean(abs(mn) >= 1.2),
       med_max = median(abs(mx)), med_mean = median(abs(mn)))
}

cat("=== Real mm10 background, sliding-25nt-max G4Hunter ===\n")
cat("(window_sizes 250/500 in config -> extracted widths 2*w+1 = 501 / 1001 bp)\n\n")
out <- list()
for (w in c(501L, 1001L)) {
  r <- calib(w)
  out[[as.character(w)]] <- r
  cat(sprintf("width %4d bp  n=%d  mean GC %.3f\n", r$width, r$n, r$gc))
  cat(sprintf("   median |sliding max|          %.3f\n", r$med_max))
  cat(sprintf("   median |window mean| (v1.0)   %.3f\n", r$med_mean))
  cat(sprintf("   FPR at literature 1.2         %.3f   <- too permissive\n", r$fpr_at_1.2))
  cat(sprintf("   FPR at literature 1.5         %.3f\n", r$fpr_at_1.5))
  cat(sprintf("   FPR at 1.2 on v1.0 mean stat  %.3f   <- the audit's fix: nothing passes\n",
              r$mean_fpr_at_1.2))
  cat(sprintf("   CALIBRATED threshold @5%% FPR  %.3f\n\n", r$thr))
}

cat("=== Recommended config ===\n")
cat(sprintf("strand.g4_threshold: %.2f   # calibrated, 5%% FPR vs random mm10 at 501 bp\n",
            out[["501"]]$thr))
cat("strand.g4_threshold_fpr: 0.05        # provenance for the number above\n")

write.csv(do.call(rbind, lapply(out, function(r) data.frame(
  width_bp = r$width, n_background = r$n, mean_gc = r$gc,
  median_abs_sliding_max = r$med_max, median_abs_window_mean = r$med_mean,
  fpr_at_1.2_slidingmax = r$fpr_at_1.2, fpr_at_1.5_slidingmax = r$fpr_at_1.5,
  fpr_at_1.2_windowmean_v10 = r$mean_fpr_at_1.2,
  calibrated_threshold_5pct_fpr = r$thr))),
  "verify/g4_threshold_calibration.csv", row.names = FALSE)
cat("\nWritten: verify/g4_threshold_calibration.csv\n")
