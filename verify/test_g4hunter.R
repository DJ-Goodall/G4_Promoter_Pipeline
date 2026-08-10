#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# v1.1 verification: does the sliding-window G4Hunter statistic actually
# discriminate at the conventional 1.2 cutoff, where the v1.0 window mean did
# not? This is the check the g4_threshold change depends on.
# Run from the repository root:  Rscript verify/test_g4hunter.R
# ---------------------------------------------------------------------------

src <- readLines("scripts/48_strand_asymmetry.R")
s <- grep("^g4_hunter_bases <- function", src)
e <- grep("^g4_hunter_score <- g4_hunter_max", src)
eval(parse(text = paste(src[s:e], collapse = "\n")))

ctl <- read.delim("data/known_topology_controls.tsv", stringsAsFactors = FALSE)
flank <- paste(rep("ATATATATAT", 25), collapse = "")   # 250 bp neutral flank

cat("=== Literature-validated G4 controls, embedded in 500 bp of neutral flank ===\n")
cat(sprintf("%-20s %10s %10s %10s\n", "control", "max(v1.1)", "mean(v1.0)", "calls G4?"))
res <- data.frame()
for (i in seq_len(nrow(ctl))) {
  emb <- paste0(flank, ctl$sequence[i], flank)
  mx <- g4_hunter_max(emb); mn <- g4_hunter_mean(emb)
  cat(sprintf("%-20s %10.3f %10.3f %10s\n", ctl$name[i], mx, mn,
              if (abs(mx) >= 1.2) "yes" else "NO"))
  res <- rbind(res, data.frame(name = ctl$name[i], max_v11 = mx, mean_v10 = mn))
}

set.seed(1)
rnd <- replicate(2000, paste(sample(c("A", "C", "G", "T"), 501, TRUE), collapse = ""))
mx <- vapply(rnd, g4_hunter_max, numeric(1), USE.NAMES = FALSE)
mn <- vapply(rnd, g4_hunter_mean, numeric(1), USE.NAMES = FALSE)

cat("\n=== Random 501 bp sequence, n = 2000 (the background the cutoff must reject) ===\n")
cat(sprintf("v1.1 sliding max : median |score| %.3f   q95 %.3f   frac >= 1.2 = %.4f\n",
            median(abs(mx)), quantile(abs(mx), .95), mean(abs(mx) >= 1.2)))
cat(sprintf("v1.0 window mean : median |score| %.3f   q95 %.3f   frac >= 1.2 = %.4f\n",
            median(abs(mn)), quantile(abs(mn), .95), mean(abs(mn) >= 1.2)))

cat("\n=== Verdict ===\n")
pass_ctl <- mean(abs(res$max_v11) >= 1.2)
pass_v10 <- mean(abs(res$mean_v10) >= 1.2)
fpr      <- mean(abs(mx) >= 1.2)

cat(sprintf("Controls called at 1.2, v1.0 window mean : %.0f%% (%d/%d)\n",
            100 * pass_v10, sum(abs(res$mean_v10) >= 1.2), nrow(res)))
cat(sprintf("Controls called at 1.2, v1.1 sliding max : %.0f%% (%d/%d)\n",
            100 * pass_ctl, sum(abs(res$max_v11) >= 1.2), nrow(res)))
cat(sprintf("Random-sequence FPR at 1.2, v1.1 sliding max: %.1f%%\n", 100 * fpr))

# Claim 1: the audit's prescribed fix (g4_threshold = 1.2 on the v1.0
# statistic) would have called NOTHING. This is the claim REVIEW.md 2.4 makes.
stopifnot(pass_v10 == 0)
cat("\n[PASS] The v1.0 window mean calls ZERO controls at 1.2 -- setting\n")
cat("       g4_threshold to the literature cutoff as the audit prescribed would\n")
cat("       have emptied g4_detection_by_strand.tsv instead of fixing it.\n")

# Claim 2: switching to the sliding max restores sensitivity...
stopifnot(pass_ctl > 0.8)
cat("\n[PASS] The sliding-25nt max recovers >80% of controls at 1.2.\n")

# Claim 3: ...but 1.2 is NOT a usable per-window cutoff, because a max over
# ~477 overlapping windows is an extreme-value statistic. This is why v1.1
# calibrates the threshold rather than hard-coding a literature value.
stopifnot(fpr > 0.10)
cat("\n[PASS] ...but at 1.2 the false-positive rate on random sequence is\n")
cat(sprintf("       %.0f%%, far too permissive for a per-window classifier. The\n",
            100 * fpr))
cat("       published 1.2/1.5 cutoffs are genome-SCANNING thresholds, not\n")
cat("       per-locus decisions. Hence calibration:\n")
cat("       see verify/calibrate_g4_threshold.R for the mm10-derived values\n")
cat("       (~2.28 at 501 bp, ~2.60 at 1001 bp for a 5% FPR).\n")
