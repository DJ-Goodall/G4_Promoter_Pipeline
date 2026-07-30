# ============================================================================
# scripts/_shared/helpers.R
# Shared constants + functions for the DeepG4 / G4ShapePredictor topology pipeline.
#
# The peak-calling / signal / meta-profile / differential helpers are lifted
# (with light edits) from g4_gloop_extended_V3.Rmd so this pipeline reproduces
# the exact behaviour of the validated Sato2025 analysis. New helpers (config
# loading, sequence extraction, topology constants) are added at the end.
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(stringr)
})

# Null-coalescing operator (base R has none)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
condition_colours <- c(
  "WT"      = "#4DAF4A",
  "DHX36KO" = "#E41A1C",
  "FANCJKO" = "#377EB8",
  "dKO"     = "#984EA3",
  "ERCCWT"  = "#FDBF6F",
  "ERCCKO"  = "#FF7F00"
)

main_genotypes <- c("WT", "DHX36KO", "FANCJKO", "dKO")

# bigWig clone code -> genotype name (from filename patterns)
clone_to_genotype <- c(
  "WT" = "WT", "P2D2" = "DHX36KO", "D1D6" = "FANCJKO", "P3D4" = "dKO",
  "ERCCWT" = "ERCCWT", "ERCCKO" = "ERCCKO"
)

assay_label_map <- c("G4_BG4" = "BG4 (G4)", "Rloop_S96" = "S9.6 (R-loop)")

# Topology classes used throughout (G4SP classes + the no-PQS bucket)
topology_levels  <- c("parallel", "antiparallel", "hybrid", "no_canonical_PQS")
topology_labels  <- c("Parallel (4+0)", "Antiparallel (2+2)", "Hybrid (3+1)", "No canonical PQS")
topology_palette <- c(
  "parallel"         = "#1B9E77",
  "antiparallel"     = "#D95F02",
  "hybrid"           = "#7570B3",
  "no_canonical_PQS" = "grey70"
)

# ---------------------------------------------------------------------------
# Config + path helpers
# ---------------------------------------------------------------------------
# Resolve the project root (DeepG4_analysis/). When sourced via Rscript from the
# project root (Snakemake workdir) this is just getwd(); we keep it explicit so
# scripts also work when launched from elsewhere.
load_config <- function(config_path = "config/config.yaml") {
  if (!requireNamespace("yaml", quietly = TRUE))
    stop("Package 'yaml' is required to read config.yaml")
  cfg <- yaml::read_yaml(config_path)
  # yaml sequences come back as R lists; coerce the ones used as character
  # vectors (factor levels, seqname indexing) so downstream code is robust.
  if (!is.null(cfg$std_chroms)) cfg$std_chroms <- as.character(unlist(cfg$std_chroms))
  if (!is.null(cfg$genotypes))  cfg$genotypes  <- as.character(unlist(cfg$genotypes))
  cfg
}

ensure_dirs <- function(...) {
  for (d in c(...)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Plotting helpers (from V3)
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  cowplot::theme_cowplot(font_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold"),
      panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(colour = "grey30", size = base_size - 1),
      legend.position  = "right"
    )
}

# save_plot writes PDF (cairo) + PNG (300 dpi) into <fig_root>/<subdir>/
save_plot <- function(plot_obj, filename, subdir, fig_root, width = 8, height = 6) {
  dir_path <- file.path(fig_root, subdir)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir_path, paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(file.path(dir_path, paste0(filename, ".png")), plot_obj,
         width = width, height = height, dpi = 300, device = "png")
  invisible(plot_obj)
}

save_base_plot <- function(plot_fn, filename, subdir, fig_root, width = 8, height = 6) {
  dir_path <- file.path(fig_root, subdir)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(dir_path, paste0(filename, ".pdf")), width = width, height = height)
  plot_fn(); dev.off()
  png(file.path(dir_path, paste0(filename, ".png")),
      width = width, height = height, units = "in", res = 300)
  plot_fn(); dev.off()
}

cache_or_build <- function(cache_file, expr, rebuild = FALSE) {
  if (!rebuild && file.exists(cache_file)) {
    message("Loading cached object: ", cache_file)
    return(readRDS(cache_file))
  }
  obj <- base::force(expr)
  saveRDS(obj, cache_file)
  message("Cached to: ", cache_file)
  obj
}

# ---------------------------------------------------------------------------
# BigWig signal extraction (from V3) — Windows drive-letter safe, in-memory cached
# ---------------------------------------------------------------------------
.bw_rle_cache <- new.env(parent = emptyenv())

default_chrom_sizes <- function(std = paste0("chr", c(1:19, "X"))) {
  gi <- GenomeInfoDb::seqlengths(
    BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  gi[std]
}

read_bigwig_rle <- function(bw_filepath, chroms = NULL, chrom_sizes = NULL) {
  abs_path <- normalizePath(bw_filepath, winslash = "/", mustWork = TRUE)
  if (exists(abs_path, envir = .bw_rle_cache)) {
    return(get(abs_path, envir = .bw_rle_cache))
  }
  if (is.null(chrom_sizes)) chrom_sizes <- default_chrom_sizes()
  lvl <- names(chrom_sizes)
  if (!is.null(chroms)) lvl <- intersect(lvl, chroms)
  if (length(lvl) == 0) return(S4Vectors::List())

  si <- GenomeInfoDb::Seqinfo(seqnames = lvl,
                              seqlengths = as.integer(chrom_sizes[lvl]))
  which_gr <- GenomicRanges::GRanges(
    seqnames = lvl,
    ranges   = IRanges::IRanges(1L, as.integer(chrom_sizes[lvl])),
    seqinfo  = si)

  if (.Platform$OS.type == "windows" && grepl("^[A-Za-z]:", abs_path)) {
    old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
    setwd(dirname(abs_path)); path_for_import <- basename(abs_path)
  } else {
    path_for_import <- abs_path
  }
  message(sprintf("  Loading BigWig: %s", basename(abs_path)))
  bwf <- rtracklayer::BigWigFile(path_for_import)
  rle <- rtracklayer::import(bwf, which = which_gr, as = "RleList")
  assign(abs_path, rle, envir = .bw_rle_cache)
  rle
}

evict_bigwigs <- function(bw_filepaths) {
  for (fp in normalizePath(bw_filepaths, winslash = "/", mustWork = FALSE)) {
    if (exists(fp, envir = .bw_rle_cache, inherits = FALSE))
      rm(list = fp, envir = .bw_rle_cache)
  }
  gc(verbose = FALSE)
}

compute_region_signal <- function(regions, bw_filepath) {
  cov <- read_bigwig_rle(bw_filepath)
  shared <- intersect(as.character(GenomicRanges::seqnames(regions)), names(cov))
  if (length(shared) == 0) return(rep(NA_real_, length(regions)))
  regions_sub <- regions[as.character(GenomicRanges::seqnames(regions)) %in% shared]
  GenomeInfoDb::seqlevels(regions_sub, pruning.mode = "coarse") <- shared
  sl <- GenomeInfoDb::seqlengths(cov)[shared]
  GenomeInfoDb::seqlengths(regions_sub) <- sl[GenomeInfoDb::seqlevels(regions_sub)]
  regions_sub <- GenomicRanges::trim(regions_sub)
  if (length(regions_sub) == 0) return(rep(NA_real_, length(regions)))
  ba <- GenomicRanges::binnedAverage(regions_sub, cov[shared], "score")
  scores <- rep(NA_real_, length(regions))
  if (!is.null(names(regions)) && !is.null(names(regions_sub))) {
    idx <- match(names(regions_sub), names(regions))
    scores[idx[!is.na(idx)]] <- ba$score
  } else {
    ov <- GenomicRanges::findOverlaps(ba, regions, type = "equal")
    scores[S4Vectors::subjectHits(ov)] <- ba$score[S4Vectors::queryHits(ov)]
  }
  scores
}

# Vectorised, strand-aware meta-profile matrix (from V3)
compute_profile_matrix <- function(regions, bw_filepath, n_bins = 100, half_width = 2000) {
  cov   <- read_bigwig_rle(bw_filepath)
  shared <- intersect(as.character(GenomicRanges::seqnames(regions)), names(cov))
  if (length(shared) == 0)
    return(list(bin_mat = matrix(NA_real_, nrow = length(regions), ncol = n_bins),
                positions = seq(-half_width, half_width, length.out = n_bins)))
  keep <- as.character(GenomicRanges::seqnames(regions)) %in% shared
  regions_sub <- regions[keep]
  GenomeInfoDb::seqlevels(regions_sub, pruning.mode = "coarse") <- shared
  sl <- GenomeInfoDb::seqlengths(cov)[shared]
  GenomeInfoDb::seqlengths(regions_sub) <- sl[GenomeInfoDb::seqlevels(regions_sub)]
  regions_sub <- GenomicRanges::trim(regions_sub)

  bin_mat <- matrix(NA_real_, nrow = length(regions), ncol = n_bins)
  for (chr in shared) {
    chr_idx <- which(as.character(GenomicRanges::seqnames(regions_sub)) == chr)
    if (length(chr_idx) == 0) next
    chr_len  <- length(cov[[chr]])
    orig_idx <- which(keep)[chr_idx]
    s_all <- pmax(1L, GenomicRanges::start(regions_sub[chr_idx]))
    e_all <- pmin(chr_len, GenomicRanges::end(regions_sub[chr_idx]))
    widths <- e_all - s_all + 1L
    ok <- widths >= n_bins
    if (!any(ok)) next
    s_ok <- s_all[ok]; e_ok <- e_all[ok]
    oi_ok <- orig_idx[ok]; ci_ok <- chr_idx[ok]
    vw <- IRanges::Views(cov[[chr]], start = s_ok, end = e_ok)
    sig_matrix <- do.call(rbind, lapply(vw, as.numeric))
    w <- ncol(sig_matrix)
    bin_idx <- cut(seq_len(w), breaks = n_bins, labels = FALSE)
    bin_mat_chr <- matrix(NA_real_, nrow = nrow(sig_matrix), ncol = n_bins)
    for (b in seq_len(n_bins)) {
      cols <- which(bin_idx == b)
      if (length(cols) > 0)
        bin_mat_chr[, b] <- rowMeans(sig_matrix[, cols, drop = FALSE], na.rm = TRUE)
    }
    strands <- as.character(GenomicRanges::strand(regions_sub[ci_ok]))
    minus_strand <- !is.na(strands) & strands == "-"
    if (any(minus_strand))
      bin_mat_chr[minus_strand, ] <- bin_mat_chr[minus_strand, seq(n_bins, 1L), drop = FALSE]
    bin_mat[oi_ok, ] <- bin_mat_chr
  }
  list(bin_mat = bin_mat,
       positions = seq(-half_width, half_width, length.out = n_bins))
}

mean_replicate_signal <- function(regions, bw_filepaths) {
  mat <- vapply(bw_filepaths,
                function(fp) compute_region_signal(regions, fp),
                numeric(length(regions)))
  if (is.null(dim(mat))) return(mat)
  rowMeans(mat, na.rm = TRUE)
}

mean_replicate_profile <- function(regions, bw_filepaths, n_bins = 100, half_width = 2000) {
  profs <- lapply(bw_filepaths, function(fp)
    compute_profile_matrix(regions, fp, n_bins = n_bins, half_width = half_width))
  arr <- simplify2array(lapply(profs, `[[`, "bin_mat"))
  if (length(dim(arr)) == 3) {
    rep_mean <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  } else {
    rep_mean <- arr
  }
  rep_mean[is.nan(rep_mean)] <- NA_real_
  bin_means <- colMeans(rep_mean, na.rm = TRUE)
  bin_sems  <- apply(rep_mean, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  data.frame(bin = seq_len(n_bins), position = profs[[1]]$positions,
             mean = bin_means, sem = bin_sems)
}

# ---------------------------------------------------------------------------
# Peak annotation + bigWig metadata (from V3)
# ---------------------------------------------------------------------------
annotate_peak_regions <- function(peaks, region_list, intergenic_label = "intergenic") {
  out <- rep(intergenic_label, length(peaks))
  for (nm in names(region_list)) {
    gr <- region_list[[nm]]
    if (length(gr) == 0) next
    hits <- GenomicRanges::findOverlaps(peaks, gr, ignore.strand = TRUE)
    idx  <- S4Vectors::queryHits(hits)
    idx  <- idx[out[idx] == intergenic_label]
    out[idx] <- nm
  }
  S4Vectors::mcols(peaks)$region <-
    factor(out, levels = c(names(region_list), intergenic_label))
  peaks
}

build_bw_meta <- function(bigwig_dir, genotypes = main_genotypes, assay = NULL) {
  bw_files <- list.files(bigwig_dir, pattern = "\\.bw$", full.names = TRUE)
  meta <- data.frame(filepath = bw_files, filename = basename(bw_files),
                     stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      gsm = stringr::str_extract(filename, "GSM\\d+"),
      assay = dplyr::case_when(
        grepl("_G4_",    filename) ~ "G4_BG4",
        grepl("_Rloop_", filename) ~ "Rloop_S96",
        TRUE ~ "other"),
      clone_code = dplyr::case_when(
        grepl("_WT_",     filename) ~ "WT",
        grepl("_P2D2_",   filename) ~ "P2D2",
        grepl("_D1D6_",   filename) ~ "D1D6",
        grepl("_P3D4_",   filename) ~ "P3D4",
        grepl("_ERCCWT_", filename) ~ "ERCCWT",
        grepl("_ERCCKO_", filename) ~ "ERCCKO",
        TRUE ~ "unknown"),
      genotype = unname(clone_to_genotype[clone_code]),
      replicate = as.integer(stringr::str_extract(filename, "(?:R|rep)(\\d+)", group = 1))
    ) %>%
    dplyr::filter(genotype %in% genotypes, assay != "other")
  if (!is.null(assay)) meta <- dplyr::filter(meta, assay == !!assay)
  meta$genotype <- factor(meta$genotype, levels = genotypes)
  meta
}

# ---------------------------------------------------------------------------
# Differential peak signal (gained / stable / lost) — from V3
# ---------------------------------------------------------------------------
peak_signal_matrix <- function(peaks, bw_filepath, chrom_sizes = NULL) {
  n <- length(peaks)
  mat <- matrix(0, nrow = n, ncol = length(bw_filepath))
  colnames(mat) <- basename(bw_filepath)
  for (j in seq_along(bw_filepath)) {
    rle <- read_bigwig_rle(bw_filepath[j], chrom_sizes = chrom_sizes)
    shared <- intersect(as.character(GenomicRanges::seqnames(peaks)), names(rle))
    if (length(shared) == 0) next
    for (chr in shared) {
      idx <- which(as.character(GenomicRanges::seqnames(peaks)) == chr)
      if (length(idx) == 0) next
      chr_len <- length(rle[[chr]])
      ir_start <- pmax(1L, GenomicRanges::start(peaks[idx]))
      ir_end   <- pmin(chr_len, GenomicRanges::end(peaks[idx]))
      # regions entirely off the chromosome clamp to end < start; skip them
      # (they keep their initialised 0 signal). In-bounds peaks all pass, so
      # this is behaviour-preserving for the peak-based callers (rules 07/12).
      ok <- ir_end >= ir_start
      if (!any(ok)) next
      vw <- IRanges::Views(rle[[chr]], IRanges::IRanges(start = ir_start[ok], end = ir_end[ok]))
      mat[idx[ok], j] <- IRanges::viewMeans(vw)
    }
  }
  mat[is.na(mat)] <- 0
  mat
}

peak_voom_lfc <- function(peaks, bw_meta, ref_genotype = "WT") {
  stopifnot(c("filepath", "genotype") %in% colnames(bw_meta))
  geno_levels <- c(ref_genotype,
                   setdiff(unique(as.character(bw_meta$genotype)), ref_genotype))
  bw_meta <- bw_meta[order(factor(bw_meta$genotype, levels = geno_levels)), , drop = FALSE]
  bw_meta$genotype <- factor(as.character(bw_meta$genotype), levels = geno_levels)

  message("Computing per-peak signal matrix (", nrow(bw_meta),
          " bigwigs x ", length(peaks), " peaks)...")
  sig_mat <- peak_signal_matrix(peaks, bw_meta$filepath)
  width_kb <- as.numeric(GenomicRanges::width(peaks)) / 1000
  count_mat <- round(sig_mat * width_kb)
  storage.mode(count_mat) <- "integer"
  keep <- rowSums(count_mat >= 1, na.rm = TRUE) >= 2
  message("Retained ", sum(keep), " / ", length(peaks),
          " peaks with >= 1 count in >= 2 samples.")

  na_df <- data.frame(logFC = rep(NA_real_, length(peaks)), AveExpr = NA_real_,
                      t = NA_real_, P.Value = NA_real_, adj.P.Val = NA_real_, B = NA_real_)
  if (sum(keep) < 10) {
    warning("Too few peaks pass the count filter; returning all-NA result.")
    return(na_df)
  }
  dge <- edgeR::DGEList(counts = count_mat[keep, , drop = FALSE])
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  design <- stats::model.matrix(~ bw_meta$genotype)
  colnames(design) <- c("(Intercept)", paste0("genotype", levels(bw_meta$genotype)[-1]))
  v   <- limma::voom(dge, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  coef_name <- paste0("genotype", setdiff(levels(bw_meta$genotype), ref_genotype))
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  out <- na_df
  out[keep, ] <- tt[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  out
}

classify_peak_dz <- function(dz_df, padj_thresh = 0.1, lfc_thresh = log2(1.5)) {
  cls <- rep("stable", nrow(dz_df))
  cls[is.na(dz_df$adj.P.Val)] <- "ambiguous"
  ix_g <- which(!is.na(dz_df$adj.P.Val) & dz_df$adj.P.Val < padj_thresh & dz_df$logFC >  lfc_thresh)
  ix_l <- which(!is.na(dz_df$adj.P.Val) & dz_df$adj.P.Val < padj_thresh & dz_df$logFC < -lfc_thresh)
  cls[ix_g] <- "gained"; cls[ix_l] <- "lost"
  factor(cls, levels = c("gained", "stable", "lost", "ambiguous"))
}

# ---------------------------------------------------------------------------
# New helpers for this pipeline
# ---------------------------------------------------------------------------
# Load the cached per-genotype + union G4 peak GRanges into a named list.
load_g4_peaks <- function(cache_v2_dir, assay = "G4_BG4",
                          genotypes = main_genotypes, include_union = TRUE) {
  keys <- genotypes
  if (include_union) keys <- c(keys, "union")
  out <- lapply(keys, function(g) {
    f <- file.path(cache_v2_dir, sprintf("peaks_%s_%s.rds", assay, g))
    if (!file.exists(f)) stop("Missing cached peak file: ", f)
    readRDS(f)
  })
  names(out) <- keys
  out
}

# Standardise a peak GRanges to UCSC std chroms and give stable peak IDs.
normalise_peaks <- function(gr, std_chroms, id_prefix = "pk") {
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepSeqlevels(
    gr, intersect(GenomeInfoDb::seqlevels(gr), std_chroms), pruning.mode = "coarse")
  GenomicRanges::strand(gr) <- "*"
  if (is.null(names(gr)) || any(duplicated(names(gr)))) {
    names(gr) <- sprintf("%s_%06d", id_prefix, seq_along(gr))
  }
  gr
}

# Resize a GRanges to fixed-width windows centred on each range's midpoint.
center_window <- function(gr, width) {
  mid <- GenomicRanges::resize(gr, width = 1, fix = "center")
  GenomicRanges::resize(mid, width = width, fix = "center")
}

# Assign each region the topology of its STRONGEST overlapping union G4 peak
# (by max_z); regions with no overlapping peak get "none". This is the single
# source of truth for the "topology of strongest peak" rule used by the
# meta-profile modules (rules 08, 14, 15). `union` must carry $topology (and
# ideally $max_z); the caller decides which groups to keep (e.g. definite + none).
# Returns `regions` with a character mcols column `group`.
assign_region_topology <- function(regions, union) {
  topo <- as.character(S4Vectors::mcols(union)$topology)
  mz <- S4Vectors::mcols(union)$max_z
  if (is.null(mz)) mz <- rep(0, length(union))
  mz[is.na(mz)] <- -Inf
  grp <- rep("none", length(regions))
  ov <- GenomicRanges::findOverlaps(regions, union, ignore.strand = TRUE)
  if (length(ov) > 0) {
    qh <- S4Vectors::queryHits(ov); sh <- S4Vectors::subjectHits(ov)
    # strongest peak per region: sort by query asc, max_z desc, keep first per query
    o  <- order(qh, -mz[sh])
    qo <- qh[o]; so <- sh[o]
    first <- !duplicated(qo)
    grp[qo[first]] <- topo[so[first]]
  }
  S4Vectors::mcols(regions)$group <- grp
  regions
}

# Collapse a multi-fragment-per-gene GRanges (mcols$gene_id, e.g. regions_5UTR.rds
# which holds several 5'UTR exon pieces per gene) to ONE TSS-proximal fragment per
# gene: + strand = lowest start, - strand = highest end (the piece closest to the
# transcript 5' end). Sets names() = gene_id for downstream joins; keeps mcols.
collapse_to_tss_proximal <- function(gr, id_col = "gene_id") {
  gid <- as.character(S4Vectors::mcols(gr)[[id_col]])
  st  <- as.character(GenomicRanges::strand(gr))
  # rank key within gene: + strand by start asc, - strand by end desc
  key <- ifelse(st == "-", -GenomicRanges::end(gr), GenomicRanges::start(gr))
  sel <- gr[order(gid, key)]
  sel <- sel[!duplicated(as.character(S4Vectors::mcols(sel)[[id_col]]))]
  names(sel) <- as.character(S4Vectors::mcols(sel)[[id_col]])
  sel
}

# Loop-length topology heuristic from pqsfinder loop lengths (ll1, ll2, ll3).
# Shared by rule 05 (top PQS per peak) and rule 27 (every PQS per peak): count how
# many of the three loops exceed `par_max`. 0 long loops -> parallel (4+0);
# >= anti_min long loops -> antiparallel (2+2); exactly one -> hybrid (3+1).
# Returns NA when there is no PQS or any loop length is missing. Defaults match
# config heuristic.{parallel_max_loop, antiparallel_min_long_loops}. NB: rule 05
# keeps its own identical local copy (closure over cfg) so it stays byte-untouched;
# this version takes the thresholds as explicit arguments for rule 27.
heuristic_topology <- function(has_pqs, ll1, ll2, ll3, par_max = 3L, anti_min = 2L) {
  if (isFALSE(has_pqs) || is.na(has_pqs)) return(NA_character_)
  loops <- c(ll1, ll2, ll3)
  if (anyNA(loops)) return(NA_character_)
  long <- sum(loops > par_max)
  if (long == 0) "parallel" else if (long >= anti_min) "antiparallel" else "hybrid"
}

# Extract DNA sequences for a GRanges from BSgenome mm10 as a DNAStringSet.
get_sequences <- function(gr, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) {
  seqs <- BSgenome::getSeq(genome, gr)
  if (!is.null(names(gr))) names(seqs) <- names(gr)
  seqs
}

write_fasta <- function(dna_string_set, path) {
  Biostrings::writeXStringSet(dna_string_set, filepath = path, format = "fasta")
  invisible(path)
}

# ---------------------------------------------------------------------------
# G4Hunter score (third sequence-propensity metric for cross-validation)
# ---------------------------------------------------------------------------
# Faithful R port of G4ShapePredictor's G4Hunter() (G4ShapePredictor.py:53):
# each base in a G-run scores +min(run,4), each base in a C-run -min(run,4),
# everything else 0; the score is the mean over the sequence length. Accepts a
# character vector or a DNAStringSet; returns a numeric vector (NA for empties).
g4hunter_one <- function(s) {
  s <- toupper(as.character(s))
  n <- nchar(s)
  if (is.na(n) || n == 0) return(NA_real_)
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  r <- rle(chars)
  vals <- ifelse(r$values == "G",  pmin(r$lengths, 4L),
          ifelse(r$values == "C", -pmin(r$lengths, 4L), 0L))
  sum(rep(vals, r$lengths)) / n
}

g4hunter_score <- function(seqs) {
  if (methods::is(seqs, "XStringSet")) seqs <- as.character(seqs)
  vapply(seqs, g4hunter_one, numeric(1), USE.NAMES = FALSE)
}

# ===========================================================================
# ADDED FOR THE COMPLETE PIPELINE — peak calling, RNA-seq loading
# (ported from 20260522_g4_gloop_extended_V3.Rmd and 20260521_..._analysis_V2.Rmd)
# ===========================================================================

# ---------------------------------------------------------------------------
# Peak calling from bigWig replicates (z-score on log1p binned signal). Ported
# verbatim from g4_gloop_extended_V3.Rmd `s1-setup-helpers`. Deterministic.
#   bin genome -> replicate-mean signal per bin -> global z on log1p(mean) ->
#   seed bins (z >= z_thresh AND signal>0 in >= min_replicates) -> merge within
#   gap_bp -> keep peaks >= min_width. Attaches n_bins / mean_signal / max_z.
# ---------------------------------------------------------------------------
call_peaks_from_bigwigs <- function(bw_filepaths, chrom_sizes, bin_size = 500L,
                                    z_thresh = 2, min_replicates = 2,
                                    gap_bp = 500, min_width = 500) {
  std_chroms <- names(chrom_sizes)
  n_reps     <- length(bw_filepaths)
  if (n_reps == 0) return(GenomicRanges::GRanges())

  cov_list <- lapply(bw_filepaths, function(fp)
    read_bigwig_rle(fp, chroms = std_chroms, chrom_sizes = chrom_sizes))

  shared <- Reduce(intersect, c(list(std_chroms), lapply(cov_list, names)))
  if (length(shared) == 0) return(GenomicRanges::GRanges())

  per_chr <- lapply(shared, function(chr) {
    chr_len <- as.integer(chrom_sizes[chr])
    starts  <- seq(1L, chr_len, by = bin_size)
    ends    <- pmin(starts + bin_size - 1L, chr_len)
    sig_mat <- vapply(cov_list, function(cov) {
      rle_chr <- cov[[chr]]
      rle_len <- length(rle_chr)
      ends_ok <- pmin(ends, rle_len)
      valid   <- ends_ok >= starts
      sig     <- numeric(length(starts))
      if (any(valid)) {
        vw         <- IRanges::Views(rle_chr, starts[valid], ends_ok[valid])
        sig[valid] <- as.numeric(IRanges::viewMeans(vw, na.rm = TRUE))
      }
      sig
    }, numeric(length(starts)))
    if (!is.matrix(sig_mat)) sig_mat <- matrix(sig_mat, ncol = n_reps)
    list(chr = chr, starts = starts, ends = ends,
         mean_sig = rowMeans(sig_mat), support = rowSums(sig_mat > 0))
  })

  all_mean    <- unlist(lapply(per_chr, `[[`, "mean_sig"), use.names = FALSE)
  all_support <- unlist(lapply(per_chr, `[[`, "support"),  use.names = FALSE)
  log_sig <- log1p(all_mean)
  mu  <- mean(log_sig, na.rm = TRUE)
  sdv <- sd(log_sig, na.rm = TRUE)
  z   <- if (sdv > 0) (log_sig - mu) / sdv else rep(0, length(log_sig))

  seed_mask <- z >= z_thresh & all_support >= min_replicates & all_mean > 0
  if (!any(seed_mask)) return(GenomicRanges::GRanges())

  chr_n      <- vapply(per_chr, function(x) length(x$starts), integer(1))
  chr_ends   <- cumsum(chr_n)
  chr_starts <- c(1L, chr_ends[-length(chr_ends)] + 1L)

  seed_parts <- lapply(seq_along(per_chr), function(k) {
    ii <- seq.int(chr_starts[k], chr_ends[k])
    mk <- seed_mask[ii]
    if (!any(mk)) return(NULL)
    GenomicRanges::GRanges(
      seqnames    = per_chr[[k]]$chr,
      ranges      = IRanges::IRanges(per_chr[[k]]$starts[mk], per_chr[[k]]$ends[mk]),
      mean_signal = per_chr[[k]]$mean_sig[mk],
      z           = z[ii][mk])
  })

  seed_gr <- do.call(c, Filter(Negate(is.null), seed_parts))
  GenomeInfoDb::seqlevels(seed_gr)  <- shared
  GenomeInfoDb::seqlengths(seed_gr) <- chrom_sizes[shared]

  merged  <- GenomicRanges::reduce(seed_gr, min.gapwidth = gap_bp, with.revmap = TRUE)
  revmap  <- S4Vectors::mcols(merged)$revmap
  S4Vectors::mcols(merged)$n_bins      <- lengths(revmap)
  S4Vectors::mcols(merged)$mean_signal <- vapply(revmap,
    function(i) mean(S4Vectors::mcols(seed_gr)$mean_signal[i]), numeric(1))
  S4Vectors::mcols(merged)$max_z       <- vapply(revmap,
    function(i) max(S4Vectors::mcols(seed_gr)$z[i]), numeric(1))
  S4Vectors::mcols(merged)$revmap <- NULL
  merged[GenomicRanges::width(merged) >= min_width]
}

# ---------------------------------------------------------------------------
# RNA-seq count loading + genotype assignment. Consolidates the identical
# count-read / metadata-regex / gene-filter logic from analysis_V2.Rmd (§1.3-1.5)
# and DeepG4 13a_deseq2.R. Returns count matrix, colData (genotype factor with
# ref first), and gene_meta (gene_id, gene_name). Columns whose prefix (before
# "-") is not in `sample_map`, or whose mapped genotype is not in `genotypes`,
# are dropped (this excludes ERCC RNA-seq columns, matching the main analyses).
# ---------------------------------------------------------------------------
load_rnaseq_counts <- function(count_table, sample_map, genotypes,
                               ref_genotype = "WT", min_count = 10) {
  if (!requireNamespace("readr", quietly = TRUE))
    stop("Package 'readr' is required for load_rnaseq_counts()")
  raw <- readr::read_tsv(count_table, show_col_types = FALSE)
  stopifnot(all(c("gene_name", "gene_id") %in% colnames(raw)))
  sample_map  <- unlist(sample_map)
  sample_cols <- setdiff(colnames(raw), c("gene_name", "gene_id"))
  prefix <- sub("-.*$", "", sample_cols)
  geno   <- unname(sample_map[prefix])
  keep   <- !is.na(geno) & geno %in% genotypes
  sample_cols <- sample_cols[keep]; geno <- geno[keep]

  cnt <- as.matrix(raw[, sample_cols]); storage.mode(cnt) <- "integer"
  # disambiguate duplicate gene names by appending gene_id (analysis_V2 §1.5)
  gname <- raw$gene_name
  dup   <- gname %in% gname[duplicated(gname)]
  gname[dup] <- paste(gname[dup], raw$gene_id[dup], sep = "_")
  rownames(cnt) <- raw$gene_id
  gene_meta <- data.frame(gene_id = raw$gene_id, gene_name = raw$gene_name,
                          label = gname, stringsAsFactors = FALSE)

  keep_g <- rowSums(cnt, na.rm = TRUE) >= min_count
  cnt <- cnt[keep_g, , drop = FALSE]; cnt[is.na(cnt)] <- 0L
  gene_meta <- gene_meta[keep_g, , drop = FALSE]

  lv <- c(ref_genotype, setdiff(genotypes, ref_genotype))
  col_data <- data.frame(genotype = factor(geno, levels = lv),
                         row.names = sample_cols)
  list(counts = cnt, coldata = col_data, gene_meta = gene_meta,
       sample_cols = sample_cols)
}

# ---------------------------------------------------------------------------
# Region constants (from g4_gloop_extended_V3.Rmd s1-setup-constants) — used by
# the regional-enrichment scripts (12-14) and the strand module (48).
# ---------------------------------------------------------------------------
extended_genotypes <- c("WT", "DHX36KO", "FANCJKO", "dKO", "ERCCWT", "ERCCKO")
region_levels  <- c("promoter", "enhancer", "5UTR", "intron1", "genebody", "intergenic")
region_labels  <- c("Promoter", "Enhancer (dELS)", "5'UTR", "1st intron",
                    "Gene body", "Intergenic")
region_palette <- setNames(
  RColorBrewer::brewer.pal(length(region_levels), "Set2"), region_levels)
