#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Standalone runner for the v1.1 heuristic control benchmark, so it can be
# executed without the full pipeline config/cache. Uses the same pqsfinder
# settings and heuristic thresholds as config/config.yaml.
# Run from the repository root:  Rscript verify/run_heuristic_benchmark.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(pqsfinder); library(Biostrings); library(S4Vectors)
})

min_score <- 20L; par_max <- 3L; anti_min <- 2L

heuristic_topology <- function(has_pqs, ll1, ll2, ll3, par_max = 3L, anti_min = 2L) {
  if (isFALSE(has_pqs) || is.na(has_pqs)) return(NA_character_)
  loops <- c(ll1, ll2, ll3)
  if (anyNA(loops)) return(NA_character_)
  long <- sum(loops > par_max)
  if (long == 0) "parallel" else if (long >= anti_min) "antiparallel" else "hybrid"
}

ctrl <- read.delim("data/known_topology_controls.tsv", stringsAsFactors = FALSE)
g4sp <- read.csv("results/tables/g4sp_control_benchmark.csv", stringsAsFactors = FALSE)

out <- data.frame()
for (i in seq_len(nrow(ctrl))) {
  s  <- Biostrings::DNAString(toupper(ctrl$sequence[i]))
  pv <- try(pqsfinder::pqsfinder(s, min_score = min_score, strand = "+"), silent = TRUE)
  if (inherits(pv, "try-error") || length(pv) == 0) {
    h <- NA_character_; l <- c(NA, NA, NA); sc <- NA
  } else {
    md <- S4Vectors::elementMetadata(pv); top <- which.max(md$score)
    l  <- c(md$ll1[top], md$ll2[top], md$ll3[top]); sc <- md$score[top]
    h  <- heuristic_topology(TRUE, l[1], l[2], l[3], par_max, anti_min)
  }
  gp <- g4sp$predicted_topology[match(ctrl$name[i], g4sp$name)]
  out <- rbind(out, data.frame(
    name = ctrl$name[i], reported = ctrl$topology[i],
    g4sp = gp, g4sp_ok = identical(gp, ctrl$topology[i]),
    heuristic = if (is.na(h)) "no_PQS" else h,
    heur_ok = !is.na(h) && identical(h, ctrl$topology[i]),
    loops = paste(l, collapse = "/"), pqs_score = sc,
    stringsAsFactors = FALSE))
}

cat("=== Topology-method accuracy on the 10 literature control G4s ===\n\n")
print(out[, c("name", "reported", "g4sp", "g4sp_ok", "heuristic", "heur_ok", "loops")],
      row.names = FALSE)
cat(sprintf("\nG4ShapePredictor : %d/%d correct\n", sum(out$g4sp_ok), nrow(out)))
cat(sprintf("Loop heuristic   : %d/%d correct\n", sum(out$heur_ok), nrow(out)))
cat(sprintf("\nAgreement between the two methods on controls: %d/%d\n",
            sum(out$g4sp == out$heuristic, na.rm = TRUE), nrow(out)))
write.csv(out, "verify/topology_control_benchmark_compare.csv", row.names = FALSE)
cat("\nWritten: verify/topology_control_benchmark_compare.csv\n")
