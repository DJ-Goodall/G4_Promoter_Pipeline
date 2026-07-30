#!/usr/bin/env Rscript
# ============================================================================
# 07_call_peaks.R   (system R)   --- Stage B ---
#
# Recompute G4 (BG4) and R-loop (S9.6) CUT&Tag peaks from the bigWigs, per assay
# x genotype, via the z-score-on-log1p-binned-signal method. Replaces the
# hand-built cache_V2 peak sets so the whole pipeline reproduces from raw data.
# Ported from 20260522_g4_gloop_extended_V3.Rmd (chunks s2-call-peaks /
# s2-union-peaks / s2-export-bed).
#
# Inputs:   config paths.bigwig_dir  (28 bigWigs: 12 G4 + 8 R-loop + 8 ERCC)
# Outputs (into cache/):
#   peaks_{G4_BG4,Rloop_S96}_{WT,DHX36KO,FANCJKO,dKO,ERCCWT,ERCCKO}.rds   (12)
#   peaks_{G4_BG4,Rloop_S96}_union.rds        (2, main 4-genotype union)
#   peaks_{G4_BG4,Rloop_S96}_ercc_union.rds   (2, ERCC union)
#   bed/<set>.bed                              (one per non-empty set)
#   results/tables/peak_summary.csv
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(rtracklayer); library(dplyr); library(readr)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "cache/bed", "results/tables")

std_chroms  <- cfg$std_chroms
pk          <- cfg$peak_calling
assays      <- as.character(unlist(pk$assays))
genotypes   <- as.character(unlist(pk$genotypes))
chrom_sizes <- default_chrom_sizes(std_chroms)

bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = NULL)
message("BigWig metadata (", nrow(bw_meta), " files):")
print(table(bw_meta$assay, bw_meta$genotype))

# --- Per (assay x genotype) peak calling ------------------------------------
peak_combos <- expand.grid(assay = assays, genotype = genotypes,
                           stringsAsFactors = FALSE)
peak_list <- list()
for (i in seq_len(nrow(peak_combos))) {
  a <- peak_combos$assay[i]; g <- peak_combos$genotype[i]
  key <- paste(a, g, sep = "_")
  cache_file <- file.path("cache", sprintf("peaks_%s.rds", key))

  peak_list[[key]] <- cache_or_build(cache_file, {
    bws <- bw_meta$filepath[bw_meta$assay == a & bw_meta$genotype == g]
    if (length(bws) == 0) {
      message(sprintf("No BigWigs for %s — empty GRanges.", key))
      return(GenomicRanges::GRanges())
    }
    min_rep <- if (length(bws) >= 3) as.integer(pk$min_replicates) else length(bws)
    message(sprintf("[%d/%d] %s (%d reps, min_rep=%d)", i, nrow(peak_combos),
                    key, length(bws), min_rep))
    p <- call_peaks_from_bigwigs(
      bws, chrom_sizes = chrom_sizes,
      bin_size = pk$bin_size, z_thresh = pk$z_thresh,
      min_replicates = min_rep, gap_bp = pk$gap_bp, min_width = pk$min_width)
    S4Vectors::mcols(p)$assay <- a; S4Vectors::mcols(p)$genotype <- g
    p
  })

  # bound RAM across the 12 combos
  evict_bigwigs(bw_meta$filepath[bw_meta$assay == a & bw_meta$genotype == g])
}

# --- Union peaks with per-genotype provenance -------------------------------
build_union <- function(assay_name, geno_set, group_rule) {
  parts <- lapply(geno_set, function(g) {
    p <- peak_list[[paste(assay_name, g, sep = "_")]]
    # keep only assay+genotype so the c() below has matching mcols across all
    # parts (an empty-peak genotype lacks n_bins/mean_signal/max_z)
    S4Vectors::mcols(p) <- S4Vectors::mcols(p)[, c("assay", "genotype"), drop = FALSE]
    p
  })
  names(parts) <- geno_set
  # unname(): do.call(c, <named list>) dispatches to base c() and returns a
  # plain list; unnamed it dispatches to the GRanges c-method (returns GRanges).
  combined <- do.call(c, unname(parts))
  u <- GenomicRanges::reduce(combined, ignore.strand = TRUE)
  for (g in geno_set) {
    ov <- GenomicRanges::findOverlaps(u, parts[[g]], ignore.strand = TRUE)
    S4Vectors::mcols(u)[[paste0("in_", g)]] <-
      seq_along(u) %in% S4Vectors::queryHits(ov)
  }
  S4Vectors::mcols(u)$group <- group_rule(u)
  S4Vectors::mcols(u)$assay <- assay_name
  u
}
# main union: provenance group replicates extended_V3 (shared/WT_only/DHX36KO_only)
main_group <- function(u) dplyr::case_when(
  S4Vectors::mcols(u)$in_WT &  S4Vectors::mcols(u)$in_DHX36KO ~ "shared",
  S4Vectors::mcols(u)$in_WT & !S4Vectors::mcols(u)$in_DHX36KO ~ "WT_only",
  !S4Vectors::mcols(u)$in_WT &  S4Vectors::mcols(u)$in_DHX36KO ~ "DHX36KO_only",
  TRUE ~ "none")
ercc_group <- function(u) dplyr::case_when(
  S4Vectors::mcols(u)$in_ERCCWT &  S4Vectors::mcols(u)$in_ERCCKO ~ "shared",
  S4Vectors::mcols(u)$in_ERCCWT & !S4Vectors::mcols(u)$in_ERCCKO ~ "ERCCWT_only",
  !S4Vectors::mcols(u)$in_ERCCWT &  S4Vectors::mcols(u)$in_ERCCKO ~ "ERCCKO_only",
  TRUE ~ "none")

main_geno <- c("WT", "DHX36KO", "FANCJKO", "dKO")
ercc_geno <- c("ERCCWT", "ERCCKO")
for (a in assays) {
  peak_list[[paste0(a, "_union")]] <- cache_or_build(
    file.path("cache", sprintf("peaks_%s_union.rds", a)),
    build_union(a, main_geno, main_group))
  peak_list[[paste0(a, "_ercc_union")]] <- cache_or_build(
    file.path("cache", sprintf("peaks_%s_ercc_union.rds", a)),
    build_union(a, ercc_geno, ercc_group))
}

# --- Export BED tracks ------------------------------------------------------
for (nm in names(peak_list)) {
  p <- peak_list[[nm]]
  if (length(p) == 0) next
  gr <- p; GenomicRanges::strand(gr) <- "*"
  rtracklayer::export(gr, file.path("cache", "bed", paste0(nm, ".bed")), format = "BED")
}

# --- Summary ----------------------------------------------------------------
peak_summary <- data.frame(
  set       = names(peak_list),
  n_peaks   = vapply(peak_list, length, integer(1)),
  mean_bp   = vapply(peak_list, function(p) if (length(p)) mean(width(p)) else NA_real_, numeric(1)),
  total_Mbp = vapply(peak_list, function(p) sum(as.numeric(width(p))) / 1e6, numeric(1)),
  row.names = NULL)
readr::write_csv(peak_summary, "results/tables/peak_summary.csv")
message("Done. Peak sets: ", nrow(peak_summary))
print(peak_summary)
