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

# Topology classes used throughout (G4SP classes + the no-PQS bucket).
#
# v1.1: "mixed" added. It is G4ShapePredictor's own abstention class -- no
# class probability cleared its precision-calibrated threshold. v1.0 never
# reached the figures with this label because rule 21 back-filled every mixed
# call with the pqsfinder loop-length heuristic; at precision 0.80 that was
# 29,225 of 130,863 PQS-bearing union peaks, and it supplied 87% of the hybrid
# and 66% of the antiparallel stratum. Adding it here means every downstream
# figure that factors on topology_levels picks it up automatically.
topology_levels  <- c("parallel", "antiparallel", "hybrid", "mixed",
                      "no_canonical_PQS")
topology_labels  <- c("Parallel (4+0)", "Antiparallel (2+2)", "Hybrid (3+1)",
                      "Mixed / unresolved", "No canonical PQS")
topology_palette <- c(
  "parallel"         = "#1B9E77",
  "antiparallel"     = "#D95F02",
  "hybrid"           = "#7570B3",
  "mixed"            = "#E7B416",
  "no_canonical_PQS" = "grey70"
)
# The three classes G4SP will commit to. Kept here so scripts stop redefining
# it locally (rule 21 and rule 43 each had their own copy).
topology_definite <- c("parallel", "antiparallel", "hybrid")

# v1.1: which topology facets a metaprofile should display, given the labels it
# actually holds.
#
# Rules 24, 30 and 33 each hard-coded `keep_groups <- c(definite, "none")` and
# then subset to it. That was harmless while rule 21 back-filled every "mixed"
# call, but under v1.1 it would DELETE every region whose strongest overlapping
# peak is mixed -- and delete it silently, since those regions do not fall
# through to "none", they leave the denominator entirely. Include the facet when
# it is populated, so the figure says "unresolved" rather than "absent".
topology_keep_groups <- function(group_vec, incl_none = TRUE) {
  c(topology_definite,
    if (any(group_vec %in% "mixed", na.rm = TRUE)) "mixed",
    if (isTRUE(incl_none)) "none")
}

# Facet labels for the above. "none" is not a topology class (it means "no G4
# peak overlaps this region"), so it is not in topology_levels and is labelled
# explicitly.
topology_group_labels <- function(keep_groups) {
  lab <- setNames(topology_labels[match(keep_groups, topology_levels)], keep_groups)
  if ("none" %in% keep_groups) lab[["none"]] <- "No G4 peak"
  lab
}

# v1.1: grid figures now annotate every panel with the n behind it. Removing the
# back-fill left some topology strata with a handful of regions, and a curve drawn
# from n = 10 is otherwise indistinguishable from one drawn from n = 216.
#
# Returns one row per facet combination with an "up N / down M" label, ready for
# geom_text(aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE). Rules 31
# and 32 each carried a near-identical local panel_counts(); this replaces both.
panel_n_labels <- function(df, facet_cols, n_col, dir_col = "direction") {
  keep <- unique(as.data.frame(df)[, c(facet_cols, dir_col, n_col), drop = FALSE])
  # Each (facet, direction) must carry ONE n. Otherwise the merge below becomes a
  # many-to-many join and silently invents labels -- e.g. passing a multi-topology
  # frame but faceting on column alone yields every up-n crossed with every down-n.
  if (anyDuplicated(keep[, c(facet_cols, dir_col), drop = FALSE]))
    stop("panel_n_labels(): ", n_col, " is not unique per (",
         paste(c(facet_cols, dir_col), collapse = ", "),
         "). Aggregate before calling.")
  side <- function(d) {
    s <- keep[as.character(keep[[dir_col]]) == d, c(facet_cols, n_col), drop = FALSE]
    names(s)[names(s) == n_col] <- paste0("n_", d)
    s
  }
  m <- merge(side("up"), side("down"), by = facet_cols, all = TRUE)
  m$label <- sprintf("up %s / down %s",
                     ifelse(is.na(m$n_up),   "0", as.character(m$n_up)),
                     ifelse(is.na(m$n_down), "0", as.character(m$n_down)))
  # merge() drops factor levels; restore them or the labels land on the wrong facet.
  for (fc in facet_cols)
    if (is.factor(df[[fc]]))
      m[[fc]] <- factor(as.character(m[[fc]]), levels = levels(df[[fc]]))
  m
}

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
  if (is.null(chrom_sizes)) chrom_sizes <- default_chrom_sizes()
  lvl <- names(chrom_sizes)
  if (!is.null(chroms)) lvl <- intersect(lvl, chroms)
  if (length(lvl) == 0) return(S4Vectors::List())
  # v1.1: key the cache on path AND the requested seqlevels. v1.0 keyed on
  # abs_path alone, so the FIRST call decided the chromosome set for every
  # later call on the same file -- a call restricted to a subset silently
  # poisoned subsequent whole-genome calls. Latent in v1.0 (all callers happen
  # to pass the same std_chroms) but a trap for anyone adding a rule.
  cache_key <- paste0(abs_path, "|", paste(lvl, collapse = ","))
  if (exists(cache_key, envir = .bw_rle_cache)) {
    return(get(cache_key, envir = .bw_rle_cache))
  }

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
  assign(cache_key, rle, envir = .bw_rle_cache)
  rle
}

evict_bigwigs <- function(bw_filepaths) {
  # v1.1: cache keys are now "<abs_path>|<seqlevels>", so evict every entry
  # whose key starts with this path. Matching on the bare path (as v1.0 did)
  # would silently free nothing and the eviction calls after each genotype
  # loop are load-bearing for memory, not hygiene.
  keys <- ls(envir = .bw_rle_cache, all.names = TRUE)
  for (fp in normalizePath(bw_filepaths, winslash = "/", mustWork = FALSE)) {
    hit <- keys[keys == fp | startsWith(keys, paste0(fp, "|"))]
    if (length(hit)) rm(list = hit, envir = .bw_rle_cache)
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
# Replicate-averaged profile MATRIX, and profiles derived from it by row subset.
#
# mean_replicate_profile() above re-reads every bigWig for every region set it is
# handed. Rules 24/30/31/32 get away with that because their strata partition a
# few thousand regions. A rule that stratifies ALL ~24k promoters five ways per
# genotype (rules 50/51) would re-read the same 12 files five times over for no
# new I/O, so it splits the work instead: rep_mean_matrix() does the read ONCE
# per genotype over ALL regions, and profile_from_rep_mean() then derives any
# stratum by subsetting rows.
#
# The two together are mean_replicate_profile() cut at its own line 339, so a
# subset reproduces it EXACTLY (bit-identical, not all.equal): colMeans() and
# sd/sqrt(n) are per-column reductions over rows, so subsetting before or after
# is the same arithmetic in the same order.
#
# That holds only if every region has the SAME WIDTH. compute_profile_matrix()
# rbinds the per-chromosome Views (line 303) and RECYCLES unequal ones, making
# `w` -- and therefore the binning grid and the short rows' contents -- depend on
# which OTHER regions share that chromosome. rep_mean_matrix() refuses unequal
# widths rather than letting a cached matrix silently diverge; callers must drop
# edge-clamped windows first.
#
# NB rule 40's local group_profile() (40_g4_rloop_cooccurrence.R:224) is the same
# idea with sqrt(sum(is.finite(z))); this uses sqrt(sum(!is.na(x))) to match
# mean_replicate_profile() at line 341 exactly (they differ only on +/-Inf).
# Folding rule 40 onto this helper is a separate change.
# ---------------------------------------------------------------------------
rep_mean_matrix <- function(regions, bw_filepaths, n_bins = 100, half_width = 2000) {
  stopifnot(length(bw_filepaths) > 0, length(regions) > 0)
  w <- unique(GenomicRanges::width(regions))
  if (length(w) != 1L)
    stop("rep_mean_matrix(): all regions must have the same width (got ",
         paste(range(w), collapse = "-"), "). compute_profile_matrix() recycles ",
         "unequal Views (helpers.R:303); drop edge-clamped windows first.")
  profs <- lapply(bw_filepaths, function(fp)
    compute_profile_matrix(regions, fp, n_bins = n_bins, half_width = half_width))
  arr <- simplify2array(lapply(profs, `[[`, "bin_mat"))
  # Single-replicate genotypes come back 2-D; same branch as line 334.
  rep_mean <- if (length(dim(arr)) == 3) apply(arr, c(1, 2), mean, na.rm = TRUE) else arr
  rep_mean[is.nan(rep_mean)] <- NA_real_          # line 339, kept verbatim
  dimnames(rep_mean) <- list(names(regions), NULL)
  attr(rep_mean, "positions") <- profs[[1]]$positions
  rep_mean
}

# (rep_mean, rows) -> data.frame(bin, position, mean, sem, n_regions).
# `rows` may be integer, logical or character (matched against rownames).
# Returns NULL for an empty selection so callers can drop the stratum rather than
# emitting a row of NaN that plots as an invisible gap.
profile_from_rep_mean <- function(rep_mean, rows = NULL, positions = NULL,
                                  half_width = 2000) {
  if (is.null(positions)) positions <- attr(rep_mean, "positions")
  if (!is.null(rows)) {
    if (is.character(rows)) {
      ix <- match(rows, rownames(rep_mean))
      if (anyNA(ix))
        stop("profile_from_rep_mean(): ", sum(is.na(ix)),
             " row name(s) not present in the matrix.")
      rows <- ix
    }
    rep_mean <- rep_mean[rows, , drop = FALSE]
  }
  if (nrow(rep_mean) == 0L) return(NULL)
  n_bins <- ncol(rep_mean)
  if (is.null(positions)) positions <- seq(-half_width, half_width, length.out = n_bins)
  stopifnot(length(positions) == n_bins)
  data.frame(bin       = seq_len(n_bins),
             position  = positions,
             mean      = colMeans(rep_mean, na.rm = TRUE),
             sem       = apply(rep_mean, 2, function(x)
                                 stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
             n_regions = nrow(rep_mean),
             row.names = NULL)
}

# Per-curve shape summary: WHERE the signal sits, as numbers rather than a
# picture. At hw2000/b80 the bin step is 50.6 bp, so peak_position differences
# below that are not resolvable -- centre_of_mass is, which is why both are here.
profile_shape_stats <- function(df) {
  m <- pmax(df$mean, 0); m[!is.finite(m)] <- 0
  tot <- sum(m)
  data.frame(
    peak_position   = if (tot > 0) df$position[which.max(m)] else NA_real_,
    centre_of_mass  = if (tot > 0) sum(df$position * m) / tot else NA_real_,
    frac_downstream = if (tot > 0) sum(m[df$position > 0]) / tot else NA_real_,
    window_mean     = mean(df$mean, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# SCALED-REGION profiles: 1 kb flank | gene body compressed to fixed bins | 1 kb
# flank. The deepTools `computeMatrix scale-regions` layout, for rules 52/53.
#
# Why this cannot reuse compute_profile_matrix(). That function bins a
# FIXED-WIDTH window, and its whole binning grid is derived from one shared
# `w <- ncol(sig_matrix)` (line 303-304) -- which is exactly the recycling trap
# rep_mean_matrix() guards against. Gene bodies here span 1 kb to 35.8 Mb
# (measured on the 24,483 promoter genes: median 15.8 kb, 5th pct 193 bp, max
# Mir3084-2 at 35,832,227 bp, an annotation span over scattered miRNA copies
# rather than a real transcription unit). Equal widths are impossible by
# construction, so the body is binned PER GENE instead.
#
# It is also a different algorithm, not just a different grid. compute_profile_
# matrix() materialises every window with as.numeric() and rbinds them -- for
# rule 50 that is 24,464 x 4000 doubles per chromosome pass. Here the span is
# ~1.3 Gb per bigWig (13x rule 50), so materialising is not an option: bins are
# built as explicit ranges and reduced with viewMeans(), which works on the Rle
# run structure and never expands it. Cost therefore scales with the number of
# RUNS touched, not with base pairs.
#
# COORDINATES. `positions` is in gene-fraction units, not bp, because bp is not
# a shared axis once bodies are scaled:
#       -1 = TSS - flank      0 = TSS      1 = TES      2 = TES + flank
# so the flanks are in units of `flank` and the body in units of gene length.
# A body bin is width/n_body_bins bp, which is 17 bp at the 1 kb minimum and
# 3.4 kb at the 95th percentile -- the x-axis is NOT linear in bp and no
# distance may be read off it. That is the price of putting TSS and TES at
# fixed positions for every gene.
# ---------------------------------------------------------------------------
scaled_positions <- function(flank = 1000L, n_flank_bins = 20L, n_body_bins = 60L) {
  c(-1 + (seq_len(n_flank_bins) - 0.5) / n_flank_bins,
         (seq_len(n_body_bins)  - 0.5) / n_body_bins,
     1 + (seq_len(n_flank_bins) - 0.5) / n_flank_bins)
}

# Bin edges for every region, in GENOMIC left-to-right order (the strand flip
# happens after the signal is read, exactly as compute_profile_matrix() does it
# at line 314). Returns integer start/end vectors of length length(regions) *
# n_bins, laid out bin-major so matrix(x, nrow = length(regions)) reassembles it.
scaled_bin_ranges <- function(regions, flank = 1000L, n_flank_bins = 20L,
                              n_body_bins = 60L) {
  s <- as.numeric(GenomicRanges::start(regions))
  e <- as.numeric(GenomicRanges::end(regions))
  w <- e - s + 1
  if (any(w < n_body_bins))
    stop("scaled_bin_ranges(): ", sum(w < n_body_bins), " region(s) are shorter ",
         "than n_body_bins (", n_body_bins, " bp), so a body bin would be empty. ",
         "Filter on width first (metaprofile.expr_gene_min_length).")
  fb <- flank / n_flank_bins
  # Left edge of every bin, plus the one final right edge: n x (n_bins + 1).
  edges <- cbind(
    outer(s - flank, (seq_len(n_flank_bins) - 1) * fb, "+"),
    s + outer(w / n_body_bins, seq_len(n_body_bins) - 1),
    outer(e + 1,     (seq_len(n_flank_bins) - 1) * fb, "+"),
    e + 1 + flank)
  nb <- 2L * n_flank_bins + n_body_bins
  st <- floor(edges[, seq_len(nb),      drop = FALSE])
  en <- floor(edges[, seq_len(nb) + 1L, drop = FALSE]) - 1
  # A scaled body bin can round to zero width on a short gene; never let a bin
  # invert, or Views() silently returns a zero-length view and the mean is NaN.
  en <- pmax(en, st)
  list(start = as.integer(st), end = as.integer(en), n_bins = nb)
}

# One bigWig -> (region x bin) matrix on the scaled grid. The read loop, the
# out-of-range policy and the strand flip now live in bin_matrix_from_ranges()
# below, shared with the composite layout so the two cannot diverge.
scaled_profile_matrix <- function(regions, bw_filepath, flank = 1000L,
                                  n_flank_bins = 20L, n_body_bins = 60L)
  bin_matrix_from_ranges(regions, bw_filepath,
    scaled_bin_ranges(regions, flank, n_flank_bins, n_body_bins))

# Replicate average of the above -- the scaled counterpart of rep_mean_matrix(),
# and consumed by the same profile_from_rep_mean(). No equal-width guard here:
# unequal widths are the POINT, and the per-gene binning makes each row
# independent of which other genes share its chromosome (which is what the
# equal-width guard buys on the fixed-width path).
scaled_rep_mean_matrix <- function(regions, bw_filepaths, flank = 1000L,
                                   n_flank_bins = 20L, n_body_bins = 60L) {
  stopifnot(length(bw_filepaths) > 0, length(regions) > 0)
  mats <- lapply(bw_filepaths, function(fp)
    scaled_profile_matrix(regions, fp, flank, n_flank_bins, n_body_bins))
  arr <- simplify2array(mats)
  rep_mean <- if (length(dim(arr)) == 3) apply(arr, c(1, 2), mean, na.rm = TRUE) else arr
  rep_mean[is.nan(rep_mean)] <- NA_real_
  dimnames(rep_mean) <- list(names(regions), NULL)
  attr(rep_mean, "positions") <- scaled_positions(flank, n_flank_bins, n_body_bins)
  rep_mean
}

# ---------------------------------------------------------------------------
# COMPOSITE (TSS-anchored | scaled interior | TES-anchored) gene profiles
#
# WHY THIS REPLACED THE PURE SCALED LAYOUT (scaled_* above), measured 2026-09-15.
# The scaled layout put 50 bp flank bins next to body bins of MEDIAN 323 bp --
# a 6.5x resolution step at exactly the TSS. That step, not biology, produced
# the vertical cliff at the TSS in every Stage N figure:
#
#   rule 50, same genes, same bigWigs, 50 bp bins, no scaling:
#     TSS -25 bp = 4.59   +25 bp = 4.53   +275 bp = 3.58   +530 bp = 2.95
#   i.e. G4 signal carries ~500 bp PAST the TSS and decays smoothly.
#
#   the scaled layout reported that entire decay as two numbers, 3.67 and 2.96,
#   and drew it as a vertical drop.
#
# Worse, a "body bin" is the gene's own length / n_body_bins, so it is a
# DIFFERENT genomic distance for every gene: 32 bp at the 5th percentile,
# 323 bp at the median, 3,352 bp at the 95th, 597 kb at the maximum. 37% of
# genes had a first body bin wider than the 500 bp decay and 21% wider than
# 1 kb, so for those the promoter-proximal shoulder was averaged into intronic
# background BEFORE the group mean was taken. The median mm10 5'UTR is 161 bp;
# 68.9% of genes had a first body bin wider than that. The 5'UTR had no
# dedicated bin anywhere on the axis and could not be read off the figure.
#
# THE FIX: only the deep interior is scaled. Both anchors get a fixed-bp window
# at a fixed bin size, so the TSS shoulder and the TES region are linear in bp
# and directly comparable to rule 50's figure (same +/-2 kb, same 50 bp bins).
#
#   [ TSS - W .. TSS + W ]  [ interior, scaled ]  [ TES - W .. TES + W ]
#     2W/anchor_bin bins       n_body_bins           2W/anchor_bin bins
#         REAL bp                 NOT bp                  REAL bp
#
# Cost: a gene needs width >= 2W + min_interior for the blocks not to overlap,
# so the short genes leave the set. At W = 2000 that is 18,219 of 21,946 (83%),
# and the loss is NOT uniform across expression groups (silent -24.2%, Q3
# -8.3%), because silent genes are shorter -- script 52 reports this per group
# rather than letting it quietly reweight the curves.
#
# PLOT COORDINATES. `x` is a single continuous axis so one ggplot geom_line can
# cross all three blocks, but it is only bp WITHIN a block:
#       -W = TSS - W      0 = TSS      +W = TSS + W
#                      [ interior occupies W .. 3W of plot space ]
#      3W = TES - W     4W = TES      5W = TES + W
# The interior is given 2W of plot width so it reads as one block of the same
# size as each anchor, not because 2W bp of anything is in there.
# ---------------------------------------------------------------------------

# Bin centres + the block each bin belongs to. `bp` is the real offset from the
# block's anchor and is NA across the interior, which is the machine-readable
# form of "no distance may be read off this stretch".
composite_positions <- function(anchor_window = 2000L, anchor_bin = 50L,
                                n_body_bins = 60L) {
  W  <- as.numeric(anchor_window)
  na <- as.integer(2 * W / anchor_bin)
  if (!isTRUE(all.equal(2 * W / anchor_bin, as.numeric(na))))
    stop("composite_positions(): 2 * anchor_window (", 2 * W, ") is not a whole ",
         "number of anchor_bin (", anchor_bin, ") bins.")
  bp_a <- -W + (seq_len(na) - 0.5) * anchor_bin
  frac <- (seq_len(n_body_bins) - 0.5) / n_body_bins
  data.frame(
    bin   = seq_len(2L * na + n_body_bins),
    x     = c(bp_a, W + frac * 2 * W, 4 * W + bp_a),
    block = factor(rep(c("tss", "interior", "tes"), c(na, n_body_bins, na)),
                   levels = c("tss", "interior", "tes")),
    bp    = c(bp_a, rep(NA_real_, n_body_bins), bp_a),
    frac  = c(rep(NA_real_, na), frac, rep(NA_real_, na)),
    stringsAsFactors = FALSE)
}

# Axis furniture for the composite x. Returned rather than hard-coded in the
# scripts so rules 52 and 53 cannot drift apart on where "TES" is drawn.
composite_axis <- function(anchor_window = 2000L) {
  W <- as.numeric(anchor_window)
  # Signed, because each anchor block has a -W and a +W edge and an unsigned
  # "2 kb" at both ends of the figure reads as the same place.
  kb <- function(v) sprintf("%s%s kb", if (v > 0) "+" else "-",
                            format(abs(v) / 1000, trim = TRUE))
  list(breaks    = c(-W, 0, W, 2 * W, 3 * W, 4 * W, 5 * W),
       labels    = c(kb(-W), "TSS", kb(W), "interior\n(scaled)", kb(-W), "TES", kb(W)),
       anchors   = c(0, 4 * W),          # dashed: the TSS and the TES
       dividers  = c(W, 3 * W),          # solid:  where bp stops/resumes meaning
       interior  = c(W, 3 * W),          # shaded span
       limits    = c(-W, 5 * W))
}

# Bin edges for every region in GENOMIC left-to-right order; the strand flip
# happens after the signal is read, exactly as compute_profile_matrix():314 and
# scaled_bin_ranges() do it. Reversing a row maps the left anchor block onto the
# right one, which is only correct because both blocks have the same bin count
# -- that is why there is a single anchor_window rather than one per end.
composite_bin_ranges <- function(regions, anchor_window = 2000L, anchor_bin = 50L,
                                 n_body_bins = 60L, min_interior = 1000L) {
  W  <- as.numeric(anchor_window)
  ab <- as.numeric(anchor_bin)
  na <- as.integer(2 * W / ab)
  s  <- as.numeric(GenomicRanges::start(regions))
  e  <- as.numeric(GenomicRanges::end(regions))
  interior <- (e - W) - (s + W) + 1
  if (any(interior < min_interior))
    stop("composite_bin_ranges(): ", sum(interior < min_interior), " region(s) leave ",
         "less than min_interior (", min_interior, " bp) between the two ",
         anchor_window, " bp anchor windows. Filter on width first ",
         "(metaprofile.expr_gene_min_length must be >= ", 2 * W + min_interior, ").")
  if (any(interior < n_body_bins))
    stop("composite_bin_ranges(): ", sum(interior < n_body_bins), " region(s) have an ",
         "interior shorter than n_body_bins (", n_body_bins, " bp), so an interior ",
         "bin would be empty.")
  # Left edge of every bin, plus the single final right edge: n x (n_bins + 1).
  # Blocks tile [s-W, s+W-1], [s+W, e-W], [e-W+1, e+W] with no overlap and no gap.
  edges <- cbind(
    outer(s - W,     (seq_len(na) - 1) * ab, "+"),
    s + W + outer(interior / n_body_bins, seq_len(n_body_bins) - 1),
    outer(e - W + 1, (seq_len(na) - 1) * ab, "+"),
    e + W + 1)
  nb <- 2L * na + n_body_bins
  st <- floor(edges[, seq_len(nb),      drop = FALSE])
  en <- floor(edges[, seq_len(nb) + 1L, drop = FALSE]) - 1
  en <- pmax(en, st)   # never let a rounded interior bin invert into a NaN view
  list(start = as.integer(st), end = as.integer(en), n_bins = nb, n_anchor_bins = na)
}

# Shared reader for every per-region binning scheme (scaled_* and composite_*).
# Split out so the out-of-range policy and the strand flip exist once: both
# layouts must treat a bigWig that ends before the chromosome does the same way,
# and a divergence there would be invisible in the figures.
bin_matrix_from_ranges <- function(regions, bw_filepath, gr) {
  cov <- read_bigwig_rle(bw_filepath)
  nb  <- gr$n_bins
  n   <- length(regions)
  bin_mat <- matrix(NA_real_, nrow = n, ncol = nb)
  n_oob <- 0L
  chrs <- as.character(GenomicRanges::seqnames(regions))
  for (chr in intersect(unique(chrs), names(cov))) {
    ci <- which(chrs == chr)
    if (!length(ci)) next
    chr_len <- length(cov[[chr]])
    # `idx` walks the n x nb layout column-major, which is how the *_bin_ranges
    # functions laid it out.
    idx <- as.vector(outer(ci, (seq_len(nb) - 1L) * n, "+"))
    st  <- gr$start[idx]; en <- gr$end[idx]
    # A bigWig's Rle can be SHORTER than the BSgenome chromosome -- measured on
    # GSM8305989 chr16, coverage ends at 98,008,968 against a BSgenome length of
    # 98,207,768, a 199 kb shortfall. Such bins stay NA rather than being
    # clamped: clamping would shrink that gene's bin and silently rescale it
    # against every other gene, the exact distortion per-gene binning avoids.
    vm  <- rep(NA_real_, length(idx))
    inb <- st >= 1L & en <= chr_len
    if (any(inb))
      vm[inb] <- IRanges::viewMeans(
        IRanges::Views(cov[[chr]], start = st[inb], end = en[inb]), na.rm = TRUE)
    n_oob <- n_oob + sum(!inb)
    bin_mat[ci, ] <- matrix(vm, nrow = length(ci), ncol = nb)
  }
  if (n_oob > 0L)
    message(sprintf("  %s: %d of %d bins (%.3f%%) lie past the end of the bigWig's ",
                    basename(bw_filepath), n_oob, n * nb, 100 * n_oob / (n * nb)),
            "coverage and are NA (not clamped).")
  # Minus-strand genes are read right-to-left, so reversing the row puts the TSS
  # on the left for every gene. Same operation as compute_profile_matrix():314.
  minus <- as.character(GenomicRanges::strand(regions)) == "-"
  if (any(minus)) bin_mat[minus, ] <- bin_mat[minus, seq(nb, 1L), drop = FALSE]
  bin_mat[is.nan(bin_mat)] <- NA_real_
  dimnames(bin_mat) <- list(names(regions), NULL)
  bin_mat
}

composite_profile_matrix <- function(regions, bw_filepath, anchor_window = 2000L,
                                     anchor_bin = 50L, n_body_bins = 60L,
                                     min_interior = 1000L)
  bin_matrix_from_ranges(regions, bw_filepath,
    composite_bin_ranges(regions, anchor_window, anchor_bin, n_body_bins, min_interior))

# Replicate average, consumed by the same profile_from_rep_mean() as every other
# layout in the repo.
composite_rep_mean_matrix <- function(regions, bw_filepaths, anchor_window = 2000L,
                                      anchor_bin = 50L, n_body_bins = 60L,
                                      min_interior = 1000L) {
  stopifnot(length(bw_filepaths) > 0, length(regions) > 0)
  mats <- lapply(bw_filepaths, function(fp)
    composite_profile_matrix(regions, fp, anchor_window, anchor_bin, n_body_bins,
                             min_interior))
  arr <- simplify2array(mats)
  rep_mean <- if (length(dim(arr)) == 3) apply(arr, c(1, 2), mean, na.rm = TRUE) else arr
  rep_mean[is.nan(rep_mean)] <- NA_real_
  dimnames(rep_mean) <- list(names(regions), NULL)
  attr(rep_mean, "positions") <-
    composite_positions(anchor_window, anchor_bin, n_body_bins)$x
  rep_mean
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
  # Fail fast on a bigWig directory that is missing or holds no .bw at all.
  #
  # Without this the pipeline does not stop -- it proceeds on an empty file set
  # and dies much later somewhere unrelated. Observed on 2026-08-10: rule 11
  # sampled 500,000 genome windows, CACHED an all-empty signal matrix, and only
  # then failed inside tidyr::pivot_longer with "`cols` must select at least one
  # column"; rule 07 reported "BigWig metadata (0 files)" and failed with "no
  # function to return from". Neither message names the real problem, and rule
  # 11 left a corrupt cache behind.
  #
  # Note this checks the DIRECTORY, not the filtered result: a legitimate call
  # can return zero rows after filtering on genotype/assay, and that is fine.
  if (length(bw_files) == 0) {
    stop("No .bw files found in bigwig_dir: '", bigwig_dir, "'\n",
         if (!dir.exists(bigwig_dir))
           "  The directory does not exist.\n"
         else
           "  The directory exists but contains no .bw files.\n",
         "  Fix paths.bigwig_dir in the config, or fetch the data with\n",
         "  scripts/download_data.ps1 (Windows) / download_data.sh.\n",
         "  To run against bigWigs you already have elsewhere, copy\n",
         "  config/config.local.example.yaml to config/config.local.yaml, set the\n",
         "  paths, and pass --configfile config/config.local.yaml to snakemake.",
         call. = FALSE)
  }
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
#
# ---------------------------------------------------------------------------
# v1.1 (2026-08-07 audit): `z_per_replicate`.
#
# In v1.0 the replicate criterion was `rowSums(sig_mat > 0) >= min_replicates`
# -- signal strictly greater than ZERO in >= 2 replicates. At CUT&Tag depth
# almost every mappable 500 bp bin satisfies that, so it is close to a no-op,
# and the z-test itself was applied once to the replicate MEAN, letting a
# single strong replicate carry a bin. The README described this as
# "replicate-aware peak calling" and "peaks present in >= 2 replicates";
# neither was accurate.
#
# The evidence that this matters is already in the v1.0 output: the ERCC
# control tracks produce MORE peaks than the real samples (G4_BG4_ERCCKO
# 114,980 vs G4_BG4_WT 107,276; ercc_union 150,977 vs union 147,843). Because
# the threshold is a relative global z on each sample's own signal, there is no
# null in which a peak can fail to be called -- the caller returns ~110k peaks
# whatever it is given.
#
# z_per_replicate = TRUE z-scores each replicate column independently and
# requires z >= z_thresh in >= min_replicates of them. The v1.0 path is kept
# (z_per_replicate = FALSE, the default) so v1.0 output stays exactly
# reproducible; set it in config as peak_calling.z_per_replicate.
#
# This does NOT add an input/IgG null, which is the deeper fix -- there is no
# IgG track in the GEO submission. Use the ERCC tracks as the empirical control
# arm instead and report the overlap (see docs/audit-2026-08/REVIEW.md section 2.2).
# ---------------------------------------------------------------------------
call_peaks_from_bigwigs <- function(bw_filepaths, chrom_sizes, bin_size = 500L,
                                    z_thresh = 2, min_replicates = 2,
                                    gap_bp = 500, min_width = 500,
                                    z_per_replicate = FALSE) {
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
    list(chr = chr, starts = starts, ends = ends, sig_mat = sig_mat,
         mean_sig = rowMeans(sig_mat), support = rowSums(sig_mat > 0))
  })

  all_mean    <- unlist(lapply(per_chr, `[[`, "mean_sig"), use.names = FALSE)
  all_support <- unlist(lapply(per_chr, `[[`, "support"),  use.names = FALSE)
  log_sig <- log1p(all_mean)
  mu  <- mean(log_sig, na.rm = TRUE)
  sdv <- sd(log_sig, na.rm = TRUE)
  z   <- if (sdv > 0) (log_sig - mu) / sdv else rep(0, length(log_sig))

  if (isTRUE(z_per_replicate)) {
    # v1.1: z-score each replicate independently, then require the threshold to
    # be met in >= min_replicates of them. This is what "peaks present in >= 2
    # replicates" actually means.
    all_sig <- do.call(rbind, lapply(per_chr, `[[`, "sig_mat"))
    rep_pass <- vapply(seq_len(ncol(all_sig)), function(j) {
      lj <- log1p(all_sig[, j])
      mj <- mean(lj, na.rm = TRUE); sj <- sd(lj, na.rm = TRUE)
      if (!is.finite(sj) || sj <= 0) return(rep(FALSE, length(lj)))
      ((lj - mj) / sj) >= z_thresh
    }, logical(nrow(all_sig)))
    if (!is.matrix(rep_pass)) rep_pass <- matrix(rep_pass, ncol = ncol(all_sig))
    n_rep_pass <- rowSums(rep_pass)
    seed_mask  <- n_rep_pass >= min_replicates & all_mean > 0
    message(sprintf(
      "  z_per_replicate: %d bins pass in >= %d/%d replicates (v1.0 rule would give %d)",
      sum(seed_mask), min_replicates, ncol(all_sig),
      sum(z >= z_thresh & all_support >= min_replicates & all_mean > 0)))
  } else {
    seed_mask <- z >= z_thresh & all_support >= min_replicates & all_mean > 0
  }
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
