#!/usr/bin/env Rscript
# ============================================================================
# 28_assign_motif_topology.R   (env: r_g4)   --- Phase 13 (all-PQS) ---
#
# Combine the per-motif G4ShapePredictor calls (rule g4sp_topology_all, which runs
# the unchanged 04 script with motif ids in the peak_id column) with the per-motif
# pqsfinder loop heuristic (rule 27) into one final topology per motif. Same rule
# as 05: take G4SP where it gives a definite class, else the loop heuristic. Every
# motif is PQS+ by construction, so there is no "no_canonical_PQS" bucket here.
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
    topology_final = case_when(
      g4sp_definite              ~ topology_g4sp,
      !is.na(topology_heuristic) ~ topology_heuristic,
      TRUE                       ~ "unclassified"),
    topology_source = case_when(
      g4sp_definite              ~ "g4sp",
      !is.na(topology_heuristic) ~ "heuristic",
      TRUE                       ~ "none"),
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
