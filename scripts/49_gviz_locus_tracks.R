#!/usr/bin/env Rscript
# ============================================================================
# 49_gviz_locus_tracks.R   (system R)   --- Stage L ---
#
# Genome-browser views (Gviz): per-genotype BigWig histogram tracks (rep 1)
# with a shared y-axis, ideogram, axis, G4 peak and gene model. Ported from
# g4_gloop_strand_analysis.Rmd chunk s9-gviz-gene-loci. The IdeogramTrack UCSC
# fetch is wrapped so the rule still succeeds offline (ideogram skipped).
#
# TWO MODES (gviz.mode):
#   genes -- the curated gviz.genes list, two views each (see below). This is
#            the default: the figure exists to show NAMED loci on a poster, and
#            the peak-ranked sample below cannot be asked for a specific gene.
#   peaks -- v1.0 behaviour: n_loci peaks sampled across the |G4Hunter score|
#            quantiles, named by the gene they fall in.
#
# The two views per gene are not redundant. Loci span 4.9 kb (Hras) to 176 kb
# (Bcl2), so no single window serves both: "tss" fixes the scale and makes peak
# widths comparable BETWEEN genes at the cost of losing the TES on long ones;
# "locus" always frames TSS to TES but changes bp-per-inch 36-fold across the
# list, which makes peak widths comparable WITHIN a gene only.
#
# Symbol resolution is strict. A name that does not map to exactly one Entrez
# id with a single-strand TxDb range is a hard stop, listing what failed and
# what near-matches exist -- because the failure this guards against (a human
# symbol like KRAS or hTERT silently yielding no gene) renders as an
# empty-but-plausible track rather than an error.
#
# Inputs:   cache/sequences_scored.rds (script 48; "peaks" mode + gene labels),
#           cache/peaks_G4_BG4_union.rds, config bigwig_dir
# Outputs:  mode "genes": results/figures/23_gviz/<Symbol>_{tss,locus}.{pdf,png}
#           mode "peaks": results/figures/23_gviz/Peak{1..n}_<gene>.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(Gviz); library(GenomicFeatures); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(GenomicRanges); library(IRanges); library(rtracklayer)
  library(BSgenome.Mmusculus.UCSC.mm10); library(dplyr); library(readr)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures/23_gviz", "results/tables")
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

# --- Union G4 peak track ----------------------------------------------------
# One track from the union peak set rather than four per-genotype tracks: the
# union carries in_<genotype> flags, presence is already readable from the
# signal tracks above it, and four extra annotation rows would crowd a 10 kb
# view. Filled by assigned topology where helpers.R provides a palette.
union_peaks <- tryCatch(readRDS("cache/peaks_G4_BG4_union.rds"),
                        error = function(e) NULL)
if (is.null(union_peaks))
  message("cache/peaks_G4_BG4_union.rds not readable -- peak track omitted.")

make_peak_track <- function(chrom, s, e) {
  if (is.null(union_peaks)) return(NULL)
  roi <- GenomicRanges::GRanges(chrom, IRanges::IRanges(s, e))
  hit <- IRanges::subsetByOverlaps(union_peaks, roi, ignore.strand = TRUE)
  if (length(hit) == 0) return(NULL)
  grp <- as.character(S4Vectors::mcols(hit)$group)
  fill <- if (!is.null(grp) && exists("topology_palette"))
    unname(ifelse(grp %in% names(topology_palette), topology_palette[grp], "#FF6B6B"))
  else "#FF6B6B"
  Gviz::AnnotationTrack(hit, name = "G4 peaks", chromosome = chrom,
    fill = fill, col = NA, stacking = "dense", shape = "box",
    background.title = "lightgray")
}

# --- One locus -> one figure ------------------------------------------------
# `extra` is an optional pre-built peak track (mode "peaks" highlights the single
# peak it selected); otherwise the union track is drawn.
plot_locus <- function(chrom, s, e, base, main, extra = NULL) {
  s <- max(1L, as.integer(s)); e <- as.integer(e)
  tracks <- list()
  ideo <- tryCatch(Gviz::IdeogramTrack(genome = "mm10", chromosome = chrom),
                   error = function(x) NULL)
  if (!is.null(ideo)) tracks <- c(tracks, list(ideo))
  tracks <- c(tracks, list(Gviz::GenomeAxisTrack()))

  rep_paths <- vapply(genotype_order,
    function(g) bw_by_genotype[[g]][gv$replicate], character(1))
  ylim <- get_max_signal(rep_paths, chrom, s, e)
  dts <- Filter(Negate(is.null), lapply(genotype_order, function(g)
    make_dtrack(bw_by_genotype[[g]][gv$replicate], g, genotype_colors[[g]],
                chrom, s, e, ylim)))
  if (length(dts) == 0) {
    message("No signal in any genotype at ", base, " -- skipped.")
    return(NULL)
  }
  tracks <- c(tracks, dts)

  pk <- if (!is.null(extra)) extra else make_peak_track(chrom, s, e)
  if (!is.null(pk)) tracks <- c(tracks, list(pk))
  tracks <- c(tracks, list(Gviz::GeneRegionTrack(txdb, chromosome = chrom,
    start = s, end = e, name = "Genes", transcriptAnnotation = "symbol",
    collapseTranscripts = "meta", fill = "#808080",
    background.title = "lightgray")))

  sizes <- c(if (!is.null(ideo)) 0.5, 0.5, rep(1, length(dts)),
             if (!is.null(pk)) 0.3, 1.5)
  draw <- function() Gviz::plotTracks(tracks, from = s, to = e, chromosome = chrom,
    sizes = sizes, transcriptAnnotation = "symbol", main = main, cex.main = 1.1)
  ok <- tryCatch({
    pdf(file.path(out_dir, paste0(base, ".pdf")), width = 14, height = 10)
    draw(); dev.off()
    png(file.path(out_dir, paste0(base, ".png")), width = 1400, height = 1000, res = 100)
    draw(); dev.off()
    message("Wrote ", base); TRUE
  }, error = function(err) {
    message(base, " failed: ", conditionMessage(err)); FALSE })
  # The y-axis is shared across the four genotype panels, so record its ceiling:
  # without it a reader cannot tell a genuinely low-signal locus from a rescaled
  # one, and these PDFs get pulled into posters away from their log.
  data.frame(figure = base, chrom = chrom, start = s, end = e,
             width_bp = e - s + 1L, ylim_max = ylim[2],
             n_genotype_tracks = length(dts), ok = ok,
             stringsAsFactors = FALSE)
}

mode <- gv$mode %||% "genes"
provenance <- list()

if (identical(mode, "genes")) {
  # --- Resolve the configured SYMBOLS, strictly ------------------------------
  want <- as.character(unlist(gv$genes))
  if (!length(want)) stop("gviz.mode is 'genes' but gviz.genes is empty.")
  ent <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = want,
           column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  genes_gr <- suppressMessages(GenomicFeatures::genes(txdb,
                                 single.strand.genes.only = TRUE))
  bad <- is.na(ent) | !(as.character(ent) %in% names(genes_gr))
  if (any(bad)) {
    # Fail loudly, and suggest the likely correct name. This list is curated
    # precisely because a human symbol (KRAS, hTERT) maps to nothing in mouse
    # and would otherwise render as an empty but entirely plausible track.
    all_sym <- AnnotationDbi::keys(org.Mm.eg.db, "SYMBOL")
    hints <- vapply(want[bad], function(w) {
      cand <- grep(paste0("^", substr(w, 1, 4)), all_sym, ignore.case = TRUE,
                   value = TRUE)
      if (length(cand)) paste(utils::head(cand, 6), collapse = ", ") else "(none)"
    }, character(1))
    stop("gviz.genes: ", sum(bad), " name(s) do not resolve to an mm10 gene:",
         paste0("\n  ", want[bad], " -> no single-strand TxDb gene; near matches: ",
                hints, collapse = ""),
         "\nUse mm10 SYMBOLS (mouse orthologues), not human names.")
  }
  sel <- genes_gr[match(as.character(ent), names(genes_gr))]
  views <- as.character(unlist(gv$views %||% list("tss", "locus")))
  pad   <- as.integer(gv$padding %||% 5000)
  gflk  <- as.integer(gv$gene_flank %||% 2000)
  message("Plotting ", length(sel), " gene(s) x ", length(views), " view(s).")

  for (i in seq_along(sel)) {
    g <- sel[i]; sym_i <- want[i]
    chrom <- as.character(GenomicRanges::seqnames(g))
    stg <- as.character(GenomicRanges::strand(g))
    tss <- if (stg == "-") GenomicRanges::end(g) else GenomicRanges::start(g)
    tes <- if (stg == "-") GenomicRanges::start(g) else GenomicRanges::end(g)
    for (v in views) {
      if (identical(v, "tss")) {
        s <- tss - pad; e <- tss + pad
        main <- sprintf("%s  TSS +/- %s bp   %s:%s-%s (%s strand)", sym_i,
          format(pad, big.mark = ","), chrom, format(s, big.mark = ","),
          format(e, big.mark = ","), stg)
      } else {
        s <- GenomicRanges::start(g) - gflk; e <- GenomicRanges::end(g) + gflk
        main <- sprintf("%s  whole locus (%s bp) +/- %s bp   %s:%s-%s (%s strand)",
          sym_i, format(GenomicRanges::width(g), big.mark = ","),
          format(gflk, big.mark = ","), chrom, format(s, big.mark = ","),
          format(e, big.mark = ","), stg)
      }
      pr <- plot_locus(chrom, s, e, sprintf("%s_%s", sym_i, v), main)
      if (!is.null(pr)) {
        pr$gene <- sym_i; pr$entrez <- as.character(ent)[i]; pr$view <- v
        pr$gene_strand <- stg; pr$gene_length_bp <- GenomicRanges::width(g)
        pr$tss <- tss; pr$tes <- tes
        provenance[[length(provenance) + 1]] <- pr
      }
    }
  }

} else {
  # --- v1.0 behaviour: n_loci peaks across the |score| quantiles -------------
  cand <- sequences_df %>% dplyr::filter(window_size == 250,
    g4_detected_strand %in% c("+", "-")) %>% dplyr::arrange(dplyr::desc(abs(g4_score)))
  if (nrow(cand) < 5) { message("Fewer than 5 candidate peaks; nothing to plot."); quit(save = "no") }
  n <- nrow(cand)
  sel <- cand[unique(c(1, round(n * 0.25), round(n * 0.5), round(n * 0.75), max(5, n - 1))), ]
  sel <- dplyr::distinct(sel, peak_start, .keep_all = TRUE)

  for (i in seq_len(min(gv$n_loci, nrow(sel)))) {
    pr0 <- sel[i, ]; chrom <- pr0$seqname; mid <- (pr0$peak_start + pr0$peak_end) / 2
    s <- max(1, mid - gv$padding); e <- mid + gv$padding
    gene_lab <- ifelse(is.na(pr0$gene_name) || pr0$gene_name == "", "NA", pr0$gene_name)
    peak_gr <- GenomicRanges::GRanges(chrom,
      IRanges::IRanges(pr0$peak_start, pr0$peak_end), strand = pr0$g4_detected_strand)
    extra <- Gviz::AnnotationTrack(peak_gr, name = "G4 Peak", chromosome = chrom,
      fill = "#FF6B6B", col = "#FF6B6B", stacking = "dense", shape = "box")
    main <- sprintf("G4 peak at %s (%s, %s-strand, score=%.2f)", gene_lab,
                    pr0$region_type, pr0$g4_detected_strand, pr0$g4_score)
    base <- sprintf("Peak%d_%s", i, gsub("[^A-Za-z0-9_]", "_", gene_lab))
    pr <- plot_locus(chrom, s, e, base, main, extra = extra)
    if (!is.null(pr)) {
      pr$gene <- gene_lab; pr$entrez <- NA_character_; pr$view <- "peak"
      pr$gene_strand <- NA_character_; pr$gene_length_bp <- NA_integer_
      pr$tss <- NA_integer_; pr$tes <- NA_integer_
      provenance[[length(provenance) + 1]] <- pr
    }
  }
}

# Provenance: exactly what was plotted, where, and at what y ceiling. Without it
# these PDFs are unattributable once they are lifted into a poster.
prov <- if (length(provenance)) dplyr::bind_rows(provenance) else
  data.frame(figure = character(0), stringsAsFactors = FALSE)
readr::write_csv(prov, "results/tables/gviz_loci.csv")
if (nrow(prov)) print(as.data.frame(prov))
if (nrow(prov) && any(!prov$ok))
  stop(sum(!prov$ok), " locus figure(s) failed to render; see the messages above.")
message("Gviz locus tracks done (mode = ", mode, ").")
