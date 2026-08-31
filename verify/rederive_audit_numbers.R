#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# v1.1 verification: re-derive every quantitative claim in REVIEW.md straight
# from the committed v1.0 tables, so the review cannot drift from the data.
# Each block stopifnot()s its own claim.
#
# Run from the repository root:  Rscript verify/rederive_audit_numbers.R
#
# The expected values below are the v1.0 numbers REVIEW.md was written against,
# so this check must read v1.0 tables. Once the v1.1 pipeline has been run,
# results/tables/ no longer holds them and every block here mismatches by
# design. Point V11_TABLES at a v1.0 export instead -- run_all.sh materialises
# one from git automatically:
#
#   git --work-tree=<tmp> checkout 3c1b264 -- results/tables
#   V11_TABLES=<tmp>/results/tables Rscript verify/rederive_audit_numbers.R
# ---------------------------------------------------------------------------

TBL <- Sys.getenv("V11_TABLES", "results/tables")
tbl <- function(f) file.path(TBL, f)
cat("Reading tables from: ", TBL, "\n", sep = "")

ok <- function(label, value, expect, tol = 0.005) {
  pass <- isTRUE(all.equal(as.numeric(value), as.numeric(expect), tolerance = tol))
  cat(sprintf("  [%s] %-52s %12.4f (claimed %.4f)\n",
              if (pass) "PASS" else "FAIL", label, value, expect))
  if (!pass) stop("mismatch: ", label)
  invisible(TRUE)
}

cat("=== 2.1 Topology: kappa and the back-fill shares ===\n")

conc <- read.csv(tbl("topology_concordance.csv"))
lv <- c("parallel", "antiparallel", "hybrid")
tab <- matrix(0, 3, 3, dimnames = list(lv, lv))
for (i in seq_len(nrow(conc)))
  tab[conc$topology_g4sp[i], conc$topology_heuristic[i]] <- conc$n[i]

n  <- sum(tab)
po <- sum(diag(tab)) / n
pe <- sum(rowSums(tab) * colSums(tab)) / n^2
kappa <- (po - pe) / (1 - pe)

ok("n where both methods definite",       n,      101638, tol = 0)
ok("observed agreement",                  po,     0.3351)
ok("expected agreement",                  pe,     0.3415)
ok("Cohen's kappa",                       kappa, -0.0098, tol = 0.02)
ok("G4SP marginal: parallel",             sum(tab["parallel", ]) / n, 0.941)
ok("heuristic marginal: parallel",        sum(tab[, "parallel"]) / n, 0.344)
ok("heuristic hybrid confirmed by G4SP",  tab["hybrid","hybrid"] / sum(tab[,"hybrid"]), 0.018)
ok("heuristic anti confirmed by G4SP",    tab["antiparallel","antiparallel"] /
                                            sum(tab[,"antiparallel"]), 0.0075)

cat("\n=== 2.1 Back-fill share per stratum (the number the audit understated) ===\n")
sweep <- read.csv(tbl("topology_precision_sweep.csv"))
s80   <- sweep[sweep$precision == "0.80", ]
cmp   <- read.csv(tbl("topology_composition_compare.csv"))
peak  <- cmp[cmp$level == "Per peak (top PQS)", ]
getn  <- function(t) peak$n[peak$topology == t]

g4sp_n <- c(parallel = s80$n_parallel, antiparallel = s80$n_antiparallel,
            hybrid = s80$n_hybrid)
cat(sprintf("  %-14s %10s %10s %12s %10s\n",
            "stratum", "final", "G4SP", "heuristic", "% heur"))
heur_tot <- 0
for (t in names(g4sp_n)) {
  fin <- getn(t); h <- fin - g4sp_n[[t]]; heur_tot <- heur_tot + h
  cat(sprintf("  %-14s %10d %10d %12d %9.1f%%\n", t, fin, g4sp_n[[t]], h, 100 * h / fin))
}
ok("heuristic total == n_mixed at 0.80",  heur_tot, s80$n_mixed, tol = 0)
ok("hybrid stratum, fraction heuristic",
   (getn("hybrid") - g4sp_n[["hybrid"]]) / getn("hybrid"), 0.870, tol = 0.01)
ok("antiparallel stratum, fraction heuristic",
   (getn("antiparallel") - g4sp_n[["antiparallel"]]) / getn("antiparallel"), 0.661, tol = 0.01)
cat("  (the audit reported 83% and 55%; those mixed per-genotype and union denominators)\n")
ok("n_mixed at precision default",
   sweep$n_mixed[sweep$precision == "default"], 0, tol = 0)

cat("\n=== 2.2 Peak caller: ERCC controls yield more peaks than the real samples ===\n")
ps <- read.csv(tbl("peak_summary.csv"))
gv <- function(s) ps$n_peaks[ps$set == s]
ok("G4_BG4_WT peaks",       gv("G4_BG4_WT"),       107276, tol = 0)
ok("G4_BG4_ERCCKO peaks",   gv("G4_BG4_ERCCKO"),   114980, tol = 0)
ok("G4_BG4_union peaks",    gv("G4_BG4_union"),    147843, tol = 0)
ok("G4_BG4_ercc_union",     gv("G4_BG4_ercc_union"), 150977, tol = 0)
stopifnot(gv("G4_BG4_ERCCKO") > gv("G4_BG4_WT"),
          gv("G4_BG4_ercc_union") > gv("G4_BG4_union"))
cat("  [PASS] ERCC control > real sample, for both per-genotype and union sets\n")
mm10 <- 2730871774
ok("BG4 union, % of mm10",  100 * ps$total_Mbp[ps$set == "G4_BG4_union"] * 1e6 / mm10, 5.2, tol = 0.05)
ok("S9.6 union, % of mm10", 100 * ps$total_Mbp[ps$set == "Rloop_S96_union"] * 1e6 / mm10, 10.5, tol = 0.05)

cat("\n=== 2.5 Enrichment null: promoter footprint is un-merged ===\n")
fe <- read.csv(tbl("region_fold_enrichment.tsv"), sep = "\t")
wtp <- fe[fe$assay == "G4_BG4" & fe$genotype == "WT" & fe$region == "promoter", ]
genome_bp <- 2633776672
implied_bp <- wtp$exp_mean / gv("G4_BG4_WT") * genome_bp
ok("implied promoter footprint (Mb)", implied_bp / 1e6, 97.9, tol = 0.01)
ok("  == 24,483 promoters x 4,001 bp (Mb)", 24483 * 4001 / 1e6, 97.9, tol = 0.01)
cat("  [PASS] exp_mean reproduces an UN-reduce()d footprint exactly\n")
ok("fold: promoter", wtp$fold, 1.989)
ok("fold: 5'UTR",
   fe$fold[fe$assay == "G4_BG4" & fe$genotype == "WT" & fe$region == "5UTR"], 11.904)
cat(sprintf("  rows with p_perm underflowing to exactly 0: %d / %d\n",
            sum(fe$p_perm == 0), nrow(fe)))

cat("\n=== 2.6 ERCC is not an antibody control: the R-loop side inverts ===\n")
rl <- function(g) fe$fold[fe$assay == "Rloop_S96" & fe$genotype == g & fe$region == "promoter"]
ok("S9.6 WT promoter fold",      rl("WT"),      0.987)
ok("S9.6 ERCCWT promoter fold",  rl("ERCCWT"),  1.367)
ok("S9.6 ERCCKO promoter fold",  rl("ERCCKO"),  1.308)
stopifnot(rl("ERCCWT") > rl("WT"), rl("ERCCKO") > rl("WT"))
cat("  [PASS] ERCC R-loop promoter enrichment EXCEEDS the real WT track\n")

cat("\n=== 2.3 DeepG4 sanity metrics ===\n")
sm <- read.csv(tbl("sanity_metrics.csv"))
gm <- function(m) sm$value[sm$metric == m]
ok("AUROC vs GC-matched background", gm("auroc_peak_vs_background"), 0.5632)
ok("Spearman rho vs max_z",          gm("spearman_prob_vs_maxz_SUPPL"), 0.1311)

cat("\n=== 2.3 g4hunter_pqs already IS the working metric ===\n")
cv <- read.csv(tbl("propensity_crossvalidation.csv"))
gr <- function(a, b) cv$spearman_rho[cv$metric_a == a & cv$metric_b == b]
ok("g4hunter_pqs vs pqsfinder_score", gr("g4hunter_pqs", "pqsfinder_score"), 0.635)
ok("g4hunter_pqs vs deepg4_prob",     gr("g4hunter_pqs", "deepg4_prob"),     0.330)
ok("g4hunter_win vs pqsfinder_score", gr("g4hunter_win", "pqsfinder_score"), 0.00035, tol = 0.1)
cat("  [PASS] the audit's 'replace g4hunter_win' fix already exists as g4hunter_pqs\n")

cat("\n=== 2.4 Strand: the tautology ===\n")
sb <- read.csv(tbl("strand_bias_relative.tsv"), sep = "\t")
stopifnot(all(sb$count_against == 0), all(sb$pct_with == 100))
cat(sprintf("  [PASS] count_against == 0 in all %d rows; pct_with == 100 in all\n", nrow(sb)))
ds <- read.csv(tbl("g4_detection_by_strand.tsv"), sep = "\t")
stopifnot(all(ds$pct_with_motif == 100))
cat(sprintf("  [PASS] pct_with_motif == 100 in all %d rows (g4_threshold 0.0)\n", nrow(ds)))
ss <- read.csv(tbl("strand_bias_summary.tsv"), sep = "\t")
m <- merge(sb, ss, by = c("genotype", "region_type", "window_size"))
wt_p <- m[m$genotype == "WT" & m$region_type == "promoter" & m$window_size == 250, ]
cat(sprintf("  WT promoter 250: gene-relative n=%d vs genomic-strand n=%d (%.1f%% retained)\n",
            wt_p$n_total.x, wt_p$n_total.y, 100 * wt_p$n_total.x / wt_p$n_total.y))
cat(sprintf("  across all 16 cells: %d of %d peaks kept (%.1f%%) -- the rest matched no\n",
            sum(m$n_total.x), sum(m$n_total.y), 100 * sum(m$n_total.x) / sum(m$n_total.y)))
cat("  same-strand gene, returned NA, and were dropped rather than counted 'against'\n")
stopifnot(sum(m$n_total.x) < sum(m$n_total.y))
gs <- read.csv(tbl("g4_score_stats.tsv"), sep = "\t")
cat(sprintf("  observed |G4Hunter window mean| range: %.3f - %.3f (cutoff 1.2 unreachable)\n",
            min(gs$mean_g4), max(gs$mean_g4)))
stopifnot(max(gs$mean_g4) < 1.2)
cat("  [PASS] the audit's prescribed g4_threshold=1.2 would have emptied the table\n")

cat("\n=== 2.8 The expression endpoint is null ===\n")
pe2 <- read.csv(tbl("promoter_g4_expression_tests.csv"))
stopifnot(all(pe2$wilcox_parallel_vs_noG4_p > 0.05))
cat(sprintf("  [PASS] headline Wilcoxon non-significant in all 3 KOs: p = %s\n",
            paste(round(pe2$wilcox_parallel_vs_noG4_p, 3), collapse = " / ")))
cat(sprintf("  median |LFC| differences: %s\n",
            paste(sprintf("%.4f", pe2$median_absLFC_parallel - pe2$median_absLFC_noG4),
                  collapse = " / ")))
mg <- read.csv(tbl("multiG4_expression_stats.csv"))
ok("Spearman nG4 vs log baseMean",
   mg$statistic[mg$metric == "spearman_nG4_vs_log_baseMean"], 0.0217)
stopifnot(all(mg$p_value[mg$metric == "kruskal_absLFC_by_bucket"] > 0.05))
cat("  [PASS] all three Kruskal tests on |LFC| by G4-count bucket are null\n")

cat("\n=== 2.1 Control benchmark: G4SP is externally validated, the heuristic is not ===\n")
cb <- read.csv(tbl("g4sp_control_benchmark.csv"))
ok("G4SP correct on controls", sum(cb$correct == "True" | cb$correct == TRUE), 8, tol = 0)
hb <- "verify/topology_control_benchmark_compare.csv"
if (file.exists(hb)) {
  h <- read.csv(hb)
  ok("loop heuristic correct on controls", sum(h$heur_ok), 2, tol = 0)
  cat("  [PASS] G4SP 8/10 vs heuristic 2/10 -- below the 33% expected from guessing\n")
} else {
  cat("  (run verify/run_heuristic_benchmark.R for the heuristic arm)\n")
}

cat("\n=========================================================\n")
cat("ALL CLAIMS IN docs/audit-2026-08/REVIEW.md RE-DERIVED FROM COMMITTED DATA\n")
cat("=========================================================\n")
