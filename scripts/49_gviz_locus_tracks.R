#!/usr/bin/env Rscript
# ============================================================================
# 49_gviz_locus_tracks.R   (system R)   --- Stage L ---
#
# Genome-browser views (Gviz) at 5 selected G4 peak loci: per-genotype BigWig
# histogram tracks (rep 1) with shared y-axis, ideogram, axis, peak, and gene
# model. Ported from g4_gloop_strand_analysis.Rmd chunk s9-gviz-gene-loci.
# gene_name (populated by script 48) drives filenames/titles. The IdeogramTrack
# UCSC fetch is wrapped so the rule still succeeds offline (ideogram skipped).
#
# Inputs:   cache/sequences_scored.rds (from script 48), config bigwig_dir
# Outputs:  results/figures/23_gviz/Peak{1..5}_<gene>.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(Gviz); library(GenomicFeatures); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(GenomicRanges); library(IRanges); library(rtracklayer)
  library(BSgenome.Mmusculus.UCSC.mm10); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures/23_gviz")
out_dir <- "results/figures/23_gviz"
gv <- cfg$gviz
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genotype_order <- c("WT", "DHX36KO", "FANCJKO", "dKO")
genotype_colors <- list(WT = "#4DAF4A", DHX36KO = "#E41A1C", FANCJKO = "#377EB8", dKO = "#984EA3")

sequences_df <- readRDS("cache/sequences_scored.rds")
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotype_order, assay = "G4_BG4")
bw_meta <- bw_meta[order(bw_meta$genotype, bw_meta$replicate), ]
bw_by_genotype <- split(bw_meta$filepath, as.character(bw_meta$genotype))

import_bw_region <- function(bw_path, chrom, s, e) tryCatch({
  bw_abs <- normalizePath(bw_path, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows" && grepl("^[A-Za-z]:", bw_abs)) {
    old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(dirname(bw_abs))
    bw_abs <- basename(bw_abs) }
  gi <- GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10)
  std <- paste0("chr", c(1:19, "X"))
  si <- GenomeInfoDb::Seqinfo(std, as.integer(gi[std]))
  roi <- GenomicRanges::GRanges(chrom, IRanges::IRanges(s, e), seqinfo = si)
  g <- rtracklayer::import(rtracklayer::BigWigFile(bw_abs), which = roi, as = "GRanges")
  if (length(g) == 0) NULL else g
}, error = function(err) NULL)

get_max_signal <- function(paths, chrom, s, e) {
  mx <- 0
  for (p in paths) { g <- import_bw_region(p, chrom, s, e)
    if (!is.null(g)) mx <- max(mx, max(g$score, na.rm = TRUE)) }
  c(0, (if (mx == 0) 1 else mx) * 1.05)
}
make_dtrack <- function(bw_path, name, color, chrom, s, e, ylim) {
  g <- import_bw_region(bw_path, chrom, s, e); if (is.null(g)) return(NULL)
  tr <- Gviz::DataTrack(range = g, name = name, type = "histogram", genome = "mm10",
    chromosome = chrom, fill.histogram = color, col.histogram = color, fill = color,
    col = color, background.title = "lightgray")
  Gviz::displayPars(tr) <- list(ylim = ylim); tr
}

# --- Select 5 peaks (250 bp window, +/- strand, by |score|) -----------------
cand <- sequences_df %>% dplyr::filter(window_size == 250,
  g4_detected_strand %in% c("+", "-")) %>% dplyr::arrange(dplyr::desc(abs(g4_score)))
if (nrow(cand) < 5) { message("Fewer than 5 candidate peaks; nothing to plot."); quit(save = "no") }
n <- nrow(cand)
sel <- cand[unique(c(1, round(n * 0.25), round(n * 0.5), round(n * 0.75), max(5, n - 1))), ]
sel <- dplyr::distinct(sel, peak_start, .keep_all = TRUE)

for (i in seq_len(min(gv$n_loci, nrow(sel)))) {
  pr <- sel[i, ]; chrom <- pr$seqname; mid <- (pr$peak_start + pr$peak_end) / 2
  s <- max(1, mid - gv$padding); e <- mid + gv$padding
  gene_lab <- ifelse(is.na(pr$gene_name) || pr$gene_name == "", "NA", pr$gene_name)
  tracks <- list()
  ideo <- tryCatch(Gviz::IdeogramTrack(genome = "mm10", chromosome = chrom), error = function(x) NULL)
  if (!is.null(ideo)) tracks <- c(tracks, list(ideo))
  tracks <- c(tracks, list(Gviz::GenomeAxisTrack()))

  ylim <- get_max_signal(vapply(genotype_order, function(g) bw_by_genotype[[g]][gv$replicate], character(1)), chrom, s, e)
  dts <- Filter(Negate(is.null), lapply(genotype_order, function(g)
    make_dtrack(bw_by_genotype[[g]][gv$replicate], g, genotype_colors[[g]], chrom, s, e, ylim)))
  tracks <- c(tracks, dts)

  peak_gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pr$peak_start, pr$peak_end),
    strand = pr$g4_detected_strand)
  tracks <- c(tracks, list(Gviz::AnnotationTrack(peak_gr, name = "G4 Peak", chromosome = chrom,
    fill = "#FF6B6B", col = "#FF6B6B", stacking = "dense", shape = "box")))
  gene_track <- Gviz::GeneRegionTrack(txdb, chromosome = chrom, start = s, end = e,
    name = "Genes", transcriptAnnotation = "symbol", collapseTranscripts = "meta",
    fill = "#808080", background.title = "lightgray")
  tracks <- c(tracks, list(gene_track))

  sizes <- c(if (!is.null(ideo)) 0.5, 0.5, rep(1, length(dts)), 0.3, 1.5)
  main <- sprintf("G4 peak at %s (%s, %s-strand, score=%.2f)", gene_lab,
                  pr$region_type, pr$g4_detected_strand, pr$g4_score)
  base <- sprintf("Peak%d_%s", i, gsub("[^A-Za-z0-9_]", "_", gene_lab))

  draw <- function() Gviz::plotTracks(tracks, from = s, to = e, chromosome = chrom,
    sizes = sizes, transcriptAnnotation = "symbol", main = main, cex.main = 1.2)
  tryCatch({
    pdf(file.path(out_dir, paste0(base, ".pdf")), width = 14, height = 10); draw(); dev.off()
    png(file.path(out_dir, paste0(base, ".png")), width = 1400, height = 1000, res = 100); draw(); dev.off()
    message("Wrote ", base)
  }, error = function(err) message("Peak ", i, " failed: ", conditionMessage(err)))
}
message("Gviz locus tracks done.")
