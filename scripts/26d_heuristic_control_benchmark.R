#!/usr/bin/env Rscript
# ============================================================================
# 26d_heuristic_control_benchmark.R   (env: r_g4)   --- NEW IN v1.1 ---
#
# Score the pqsfinder LOOP-LENGTH HEURISTIC on exactly the same literature
# control quadruplexes that 26b_control_benchmark.py scores G4ShapePredictor
# on, so the two methods can be compared on a common external standard.
#
# Why this exists
# ---------------
# The 2026-08-07 audit argued that v1.0 should stop back-filling G4SP's
# "mixed" abstentions with the loop-length heuristic, on the grounds that the
# heuristic "would be defensible if it were a validated stand-in. It is not."
# That is an assertion, and the repository already contains the apparatus to
# turn it into a measurement: data/known_topology_controls.tsv plus the
# existing G4SP benchmark (8/10 correct at precision 0.80).
#
# The heuristic had never been run through the same test. This script does it.
# If the heuristic scores comparably to G4SP, back-filling is defensible and
# the audit's headline finding weakens. If it does not, dropping the back-fill
# is an evidence-based decision rather than an assumption.
#
# Method: run pqsfinder with the pipeline's configured min_score on each
# control sequence, take the top-scoring PQS (the same rule rule 18 uses for
# genomic peaks), read ll1/ll2/ll3, and apply heuristic_topology() from
# scripts/_shared/helpers.R with the configured thresholds.
#
# Inputs:   calibration.controls_tsv, pqsfinder.min_score, heuristic.*
#           results/tables/g4sp_control_benchmark.csv  (optional; merged if present)
# Output:   results/tables/topology_control_benchmark_compare.csv
# ============================================================================

suppressPackageStartupMessages({
  library(pqsfinder); library(Biostrings); library(dplyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables")

controls_tsv <- if (!is.null(cfg$calibration$controls_tsv))
  cfg$calibration$controls_tsv else "data/known_topology_controls.tsv"
min_score <- if (!is.null(cfg$pqsfinder$min_score)) as.integer(cfg$pqsfinder$min_score) else 20L
par_max   <- if (!is.null(cfg$heuristic$parallel_max_loop)) cfg$heuristic$parallel_max_loop else 3L
anti_min  <- if (!is.null(cfg$heuristic$antiparallel_min_long_loops))
  cfg$heuristic$antiparallel_min_long_loops else 2L

ctrl <- readr::read_tsv(controls_tsv, show_col_types = FALSE)
message("Controls: ", nrow(ctrl), " sequences; pqsfinder min_score = ", min_score,
        "; heuristic par_max = ", par_max, ", anti_min = ", anti_min)

rows <- lapply(seq_len(nrow(ctrl)), function(i) {
  s <- Biostrings::DNAString(toupper(ctrl$sequence[i]))
  pv <- try(pqsfinder::pqsfinder(s, min_score = min_score, strand = "+"), silent = TRUE)
  base <- data.frame(name = ctrl$name[i], reported_topology = ctrl$topology[i],
                     stringsAsFactors = FALSE)
  if (inherits(pv, "try-error") || length(pv) == 0) {
    return(cbind(base, data.frame(pqs_found = FALSE, pqs_score = NA_integer_,
      ll1 = NA_integer_, ll2 = NA_integer_, ll3 = NA_integer_,
      heuristic_topology = NA_character_, heuristic_correct = NA)))
  }
  md  <- S4Vectors::elementMetadata(pv)
  top <- which.max(md$score)
  h <- heuristic_topology(TRUE, md$ll1[top], md$ll2[top], md$ll3[top],
                          par_max = par_max, anti_min = anti_min)
  cbind(base, data.frame(pqs_found = TRUE, pqs_score = as.integer(md$score[top]),
    ll1 = as.integer(md$ll1[top]), ll2 = as.integer(md$ll2[top]),
    ll3 = as.integer(md$ll3[top]),
    heuristic_topology = h,
    heuristic_correct = !is.na(h) && identical(h, ctrl$topology[i])))
})
heur <- dplyr::bind_rows(rows)

# --- Merge with the G4SP benchmark if rule 26b has already run --------------
g4sp_path <- "results/tables/g4sp_control_benchmark.csv"
if (file.exists(g4sp_path)) {
  g4sp <- readr::read_csv(g4sp_path, show_col_types = FALSE) %>%
    select(name, g4sp_topology = predicted_topology,
           g4sp_confidence = confidence, g4sp_correct = correct)
  cmp <- dplyr::left_join(heur, g4sp, by = "name")
} else {
  message("NOTE: ", g4sp_path, " not found -- writing the heuristic arm only. ",
          "Run rule 26b (control_benchmark) first for the side-by-side.")
  cmp <- heur %>% mutate(g4sp_topology = NA_character_,
                         g4sp_confidence = NA_real_, g4sp_correct = NA)
}

readr::write_csv(cmp, "results/tables/topology_control_benchmark_compare.csv")

# --- Report -----------------------------------------------------------------
n <- nrow(cmp)
h_ok <- sum(cmp$heuristic_correct, na.rm = TRUE)
g_ok <- sum(cmp$g4sp_correct, na.rm = TRUE)
h_called <- sum(!is.na(cmp$heuristic_topology))

message("\n=== Topology method accuracy on literature controls ===")
message(sprintf("  G4ShapePredictor : %d/%d correct", g_ok, n))
message(sprintf("  Loop heuristic   : %d/%d correct (%d of %d had a PQS at min_score %d)",
                h_ok, n, h_called, n, min_score))
print(as.data.frame(cmp[, c("name", "reported_topology", "g4sp_topology",
                            "heuristic_topology", "ll1", "ll2", "ll3")]))

message("\nInterpretation: the loop heuristic classifies on loop length alone, ",
        "which is a real determinant of topology but is measured here on the ",
        "top PQS of an isolated sequence -- no cation, no crowding, no ",
        "chromatin. Compare the two accuracies before deciding whether the ",
        "heuristic may stand in for G4SP where G4SP abstains (v1.1 says no; ",
        "see REVIEW.md 2.1).")
