#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# v1.1 verification: the width-aware, merged-region expectation.
#
# Checks the analytic expectation against brute-force simulation, then shows
# what it does to the v1.0 fold-enrichment numbers.
# Run from the repository root:  Rscript verify/test_enrichment_null.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges); library(IRanges) })

src <- readLines("scripts/12_peak_region_enrichment.R")
s <- grep("^expected_p_for_width <- function", src)
# closing brace of the function is the first line that is exactly "}" at col 1
e <- s + which(src[(s + 1L):length(src)] == "}")[1]
eval(parse(text = paste(src[s:e], collapse = "\n")))

# --- 1. Analytic expectation vs brute-force simulation ----------------------
set.seed(11)
G <- 1e6L
n_reg <- 200L
st <- sort(sample.int(G - 5000L, n_reg))
rg <- GenomicRanges::reduce(GRanges("chrT", IRanges(st, width = 200L)))

cat("=== Analytic expectation vs simulation (chrT, 1 Mb, 200 regions x 200 bp) ===\n")
cat(sprintf("%8s %14s %14s %8s\n", "peak_w", "analytic", "simulated", "ratio"))
ok <- TRUE
for (w in c(1L, 50L, 200L, 966L, 2000L)) {
  pa <- expected_p_for_width(rg, w, G)
  starts <- sample.int(G - w, 200000L, replace = TRUE)
  pk <- GRanges("chrT", IRanges(starts, width = w))
  ps <- mean(GenomicRanges::countOverlaps(pk, rg) > 0)
  cat(sprintf("%8d %14.6f %14.6f %8.3f\n", w, pa, ps, pa / ps))
  if (abs(pa / ps - 1) > 0.03) ok <- FALSE
}
stopifnot(ok)
cat("PASS: analytic expectation matches simulation within 3% at every width.\n\n")

# --- 2. What this does to the committed v1.0 numbers ------------------------
# Region widths are the mm10 sets the pipeline builds; approximated here from
# the committed exp_mean values so this runs without cache/.
fe <- read.delim("results/tables/region_fold_enrichment.tsv", stringsAsFactors = FALSE)
wt <- fe[fe$assay == "G4_BG4" & fe$genotype == "WT", ]
genome_bp <- 2633776672        # chr1-19 + chrX, mm10
n_peaks   <- 107276            # G4_BG4_WT, peak_summary.csv
peak_w    <- 909               # mean_bp, peak_summary.csv

# v1.0 region footprint implied by exp_mean, and the typical feature width
region_bp <- wt$exp_mean / n_peaks * genome_bp
n_feature <- c(promoter = 24483, enhancer = 104000, `5UTR` = 21000,
               intron1 = 21000, genebody = 21000)[wt$region]
mean_w    <- region_bp / n_feature

cat("=== Effect on the committed v1.0 BG4/WT fold enrichments ===\n")
cat("(feature counts approximate; this shows the DIRECTION and rough size)\n\n")
cat(sprintf("%-10s %8s %10s %12s %10s\n",
            "region", "mean_w", "fold_v1.0", "width_factor", "fold_v1.1"))
for (i in seq_len(nrow(wt))) {
  wf <- (mean_w[i] + peak_w - 1) / mean_w[i]      # width-dilation inflation
  cat(sprintf("%-10s %8.0f %10.2f %12.2f %10.2f\n",
              wt$region[i], mean_w[i], wt$fold[i], wf, wt$fold[i] / wf))
}
cat("\nThe monotone-in-size artefact largely disappears: 5'UTR falls from ~11.9\n")
cat("toward the promoter value rather than dwarfing it. Un-merging promoters\n")
cat("(97.9 Mb claimed footprint vs a smaller true one) pushes the other way and\n")
cat("is corrected simultaneously in the script.\n")
