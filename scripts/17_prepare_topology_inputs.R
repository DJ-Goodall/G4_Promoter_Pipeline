#!/usr/bin/env Rscript
# ============================================================================
# 01_prepare_inputs.R   (env: r_g4)
#
# Load the cached G4 (BG4) CUT&Tag peak sets + promoter regions produced by the
# validated g4_gloop_extended_V3.Rmd pipeline, standardise them, attach a per-peak
# strength value, build 201 bp peak-centred windows for DeepG4, and emit a peak
# catalogue CSV.
#
# Inputs  (paths from config.yaml):
#   cache_V2/peaks_G4_BG4_{WT,DHX36KO,FANCJKO,dKO,union}.rds
#   cache_V2/regions_promoter.rds
# Outputs:
#   cache/peaks.rds                 list(union, per_genotype, union_windows)
#   cache/promoters.rds             strand-aware promoter GRanges (std chroms)
#   results/tables/peak_catalog.csv one row per union peak
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

std_chroms  <- cfg$std_chroms
assay       <- cfg$assay
genotypes   <- cfg$genotypes
deepg4_w    <- cfg$sequences$deepg4_width

ensure_dirs("cache", "results/tables")

# --- Load cached peaks (per-genotype + union) ------------------------------
message("Loading cached G4 peak sets from: ", cfg$paths$cache_v2)
peaks_raw <- load_g4_peaks(cfg$paths$cache_v2, assay = assay,
                           genotypes = genotypes, include_union = TRUE)

union <- normalise_peaks(peaks_raw[["union"]], std_chroms, id_prefix = "g4")
per_genotype <- lapply(genotypes, function(g)
  normalise_peaks(peaks_raw[[g]], std_chroms, id_prefix = tolower(g)))
names(per_genotype) <- genotypes

# Optional smoke-test subsample of union peaks
if (!is.null(cfg$peak_subsample)) {
  n <- as.integer(cfg$peak_subsample)
  if (length(union) > n) {
    set.seed(cfg$peak_subsample_seed %||% 7)
    union <- union[sort(sample(length(union), n))]
    message("Subsampled union peaks to ", n, " (smoke test).")
  }
}

# --- Attach per-peak strength to union peaks -------------------------------
# Union peaks were built via reduce() and lost mean_signal / max_z. Recover a
# representative strength for each union peak = max over per-genotype peaks that
# overlap it (across all genotypes).
union$max_z       <- NA_real_
union$mean_signal <- NA_real_
for (g in genotypes) {
  pg <- per_genotype[[g]]
  if (length(pg) == 0) next
  if (!all(c("max_z", "mean_signal") %in% colnames(mcols(pg)))) next
  ov <- findOverlaps(union, pg, ignore.strand = TRUE)
  if (length(ov) == 0) next
  qh <- queryHits(ov); sh <- subjectHits(ov)
  mz <- tapply(pg$max_z[sh],       qh, max, na.rm = TRUE)
  ms <- tapply(pg$mean_signal[sh], qh, max, na.rm = TRUE)
  qi <- as.integer(names(mz))
  union$max_z[qi]       <- pmax(union$max_z[qi],       as.numeric(mz), na.rm = TRUE)
  union$mean_signal[qi] <- pmax(union$mean_signal[qi], as.numeric(ms), na.rm = TRUE)
}

# --- 201 bp peak-centred windows for DeepG4 --------------------------------
union_windows <- center_window(union, width = deepg4_w)
names(union_windows) <- names(union)

# --- Promoters (keep strand for TSS meta-profiles) -------------------------
message("Loading promoters: ", cfg$paths$promoters_rds)
promoters <- readRDS(cfg$paths$promoters_rds)
GenomeInfoDb::seqlevelsStyle(promoters) <- "UCSC"
promoters <- keepSeqlevels(
  promoters, intersect(seqlevels(promoters), std_chroms), pruning.mode = "coarse")

# --- Save -------------------------------------------------------------------
saveRDS(list(union = union, per_genotype = per_genotype,
             union_windows = union_windows),
        "cache/peaks.rds")
saveRDS(promoters, "cache/promoters.rds")

membership <- as.data.frame(mcols(union)[, intersect(
  c("in_WT", "in_DHX36KO", "in_FANCJKO", "in_dKO", "group"),
  colnames(mcols(union))), drop = FALSE])

catalog <- data.frame(
  peak_id     = names(union),
  chr         = as.character(seqnames(union)),
  start       = start(union),
  end         = end(union),
  width       = width(union),
  max_z       = union$max_z,
  mean_signal = union$mean_signal,
  stringsAsFactors = FALSE
)
catalog <- cbind(catalog, membership)
readr::write_csv(catalog, "results/tables/peak_catalog.csv")

message(sprintf("Done. Union peaks: %d | per-genotype: %s | promoters: %d",
                length(union),
                paste(sprintf("%s=%d", genotypes, vapply(per_genotype, length, 0L)),
                      collapse = ", "),
                length(promoters)))
