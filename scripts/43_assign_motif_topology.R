#!/usr/bin/env Rscript
# ============================================================================
# 28_assign_motif_topology.R   (env: r_g4)   --- Phase 13 (all-PQS) ---
#
# Assign one topology per PQS motif from the per-motif G4ShapePredictor calls
# (rule g4sp_topology_all, which runs the unchanged 20 script with motif ids in
# the peak_id column), and report concordance with the per-motif pqsfinder loop
# heuristic (rule 27). Every motif is PQS+ by construction, so there is no
# "no_canonical_PQS" bucket here.
#
# ---------------------------------------------------------------------------
# v1.1 (2026-08-07 audit): this script carried the SAME back-fill as rule 21 --
#
#     topology_final = case_when(
#       g4sp_definite              ~ topology_g4sp,
#       !is.na(topology_heuristic) ~ topology_heuristic,   # <- removed in v1.1
#       TRUE                       ~ "unclassified")
#
# The audit named rule 21, but the motif level has to move with it. These two
# tables are compared directly against each other in
# 44_peak_motif_distributions.R (figure topology_composition_motif_vs_peak,
# table topology_composition_compare.csv), which asks whether counting every
# motif shifts the parallel-dominated peak-level picture. De-chimerising only
# the peak arm would answer that question with one arm of pure G4SP calls and
# one arm still 87%-heuristic -- a worse comparison than leaving both alone.
#
# So: G4SP's abstentions stay "mixed" here too. The heuristic call is still
# computed and kept in its own column for concordance, but never becomes the
# reported topology. Evidence: on the 10 literature control quadruplexes G4SP
# scores 8/10 and the loop heuristic 2/10, below the 33% expected from guessing
# at three classes (26d_heuristic_control_benchmark.R).
# ---------------------------------------------------------------------------
#
# Inputs:   results/tables/motif_catalog.csv          (rule 27; topology_heuristic)
#           results/tables/motif_g4sp_topology.csv     (04 on motif ids)
#           cache/motifs_all_base.rds                   (rule 27)
# Outputs:  results/tables/motif_topology.csv
#           cache/motifs_all.rds                        (base GRanges + $topology)
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(S4Vectors) })

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables")

definite <- c("parallel", "antiparallel", "hybrid")

cat_df <- readr::read_csv("results/tables/motif_catalog.csv", show_col_types = FALSE)
g4sp   <- readr::read_csv("results/tables/motif_g4sp_topology.csv", show_col_types = FALSE) %>%
  dplyr::rename(motif_id = peak_id, topology_g4sp = g4sp_class)   # 04 wrote motif id under peak_id

merged <- cat_df %>%
  select(motif_id, peak_id, topology_heuristic) %>%
  left_join(select(g4sp, motif_id, topology_g4sp, g4sp_confidence), by = "motif_id") %>%
  mutate(
    g4sp_definite = topology_g4sp %in% definite,
    # v1.1: NO heuristic back-fill (see header). "mixed" = G4SP returned a call
    # but declined to commit to a class; "unclassified" = no G4SP row at all.
    topology_final = case_when(
      g4sp_definite            ~ topology_g4sp,
      !is.na(topology_g4sp)    ~ "mixed",
      TRUE                     ~ "unclassified"),
    topology_source = case_when(
      g4sp_definite            ~ "g4sp",
      !is.na(topology_g4sp)    ~ "g4sp_abstained",
      TRUE                     ~ "none"),
    agree = g4sp_definite & !is.na(topology_heuristic) &
            (topology_g4sp == topology_heuristic))

readr::write_csv(
  merged %>% select(motif_id, peak_id, topology_final, topology_g4sp,
                    topology_heuristic, topology_source, g4sp_confidence, agree),
  "results/tables/motif_topology.csv")

# Attach $topology to the motif GRanges for downstream rules.
motifs_gr <- readRDS("cache/motifs_all_base.rds")
S4Vectors::mcols(motifs_gr)$topology <-
  merged$topology_final[match(names(motifs_gr), merged$motif_id)]
saveRDS(motifs_gr, "cache/motifs_all.rds")

message("Done. Motif topology composition:")
print(merged %>% count(topology_final, name = "n") %>% mutate(pct = round(100 * n / sum(n), 1)))
message("Call source:")
print(merged %>% count(topology_source, name = "n"))
