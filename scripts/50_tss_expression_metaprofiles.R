#!/usr/bin/env Rscript
# ============================================================================
# 50_tss_expression_metaprofiles.R   (env: r_g4)   --- Stage M, point 1 ---
#
# Do genes with LOW vs HIGH expression carry their promoter G4 signal at
# different POSITIONS? TSS-centred G4 CUT&Tag meta-profiles over ALL promoters,
# stratified by RNA-seq expression quartile plus a silent group, one panel per
# genotype (WT, DHX36KO, FANCJKO, dKO). This is a different question from rules
# 31/32 (up- vs down-regulated DEGs); here every annotated gene participates.
#
# THIS IS THE ONLY BIGWIG PASS OF THE 50/51 PAIR. One replicate-averaged
# (promoter x bin) matrix per genotype is computed once over ALL promoters and
# cached; every stratum is then a row subset of it (rep_mean_matrix() /
# profile_from_rep_mean(), helpers.R). Calling mean_replicate_profile() per
# (genotype x stratum) instead would re-read the same 12 bigWigs five times over
# for no new I/O. Rule 51 reads the cache and never opens a bigWig.
#
# Expression groups (see the metaprofile.expr_* keys in config.yaml):
#   Q1..Q4  quartiles of mean rlog, cut per genotype from its own 3 replicates
#   silent  gene IS in the raw count table but was dropped by rnaseq.min_count
# Promoters whose symbol is absent from the RNA-seq vocabulary have UNKNOWN
# expression, not zero: they are dropped and counted, never folded into silent.
#
# Inputs:   cfg$paths$promoters_rds, cache/rld.rds, cache/sample_info.rds,
#           cache/gene_meta.rds, cfg$rnaseq$count_table, G4 bigWigs
# Outputs:  cache/tss_expr_repmean_<assay>_hw<HW>_b<N>.rds   (rule 51 reads this)
#           results/tables/tss_expression_gene_groups.csv
#           results/tables/tss_expression_counts.csv
#           results/tables/tss_expression_quartile_migration.csv
#           results/tables/tss_expression_profiles.csv
#           results/tables/tss_expression_shape_stats.csv
#           results/figures/24_expression_metaprofiles/
#             tss_signal_by_expression.{pdf,png}
#             tss_signal_by_expression_shape.{pdf,png}
#             tss_signal_by_expression_fixedcuts.{pdf,png}
#             tss_expression_group_sizes.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(SummarizedExperiment)
  library(dplyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

fig_subdir <- "24_expression_metaprofiles"
fig_root   <- "results/figures"
ensure_dirs("cache", "results/tables", file.path(fig_root, fig_subdir))

mp         <- cfg$metaprofile
rc         <- cfg$rnaseq
genotypes  <- cfg$genotypes
ref_geno   <- cfg$ref_genotype
half_width <- mp$half_width %||% 2000
n_bins     <- mp$n_bins %||% 80
n_q        <- as.integer(mp$expr_n_quantiles %||% 4)
per_geno   <- isTRUE(mp$expr_per_genotype %||% TRUE)
fixed_arm  <- isTRUE(mp$expr_fixed_cuts_arm %||% TRUE) && per_geno
basis_mode <- mp$expr_quantile_basis %||% "matched"
incl_silent <- isTRUE(mp$expr_include_silent %||% TRUE)
q_levels   <- paste0("Q", seq_len(n_q))

cache_file <- sprintf("cache/tss_expr_repmean_%s_hw%d_b%d.rds",
                      cfg$assay, half_width, n_bins)

# --- 1. Promoter TSS windows -------------------------------------------------
prom <- readRDS(cfg$paths$promoters_rds)
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), cfg$std_chroms), pruning.mode = "coarse")
GenomeInfoDb::seqlengths(prom) <- default_chrom_sizes(cfg$std_chroms)[GenomeInfoDb::seqlevels(prom)]
win   <- GenomicRanges::trim(center_window(prom, width = 2 * half_width))
n_pre <- length(win)
# Every window must be exactly 2*half_width. compute_profile_matrix() rbinds the
# per-chromosome Views (helpers.R:303) and RECYCLES unequal ones, which would make
# a row's bins depend on which other rows share its chromosome -- and so make the
# cached matrix disagree with mean_replicate_profile() on a subset. Measured
# 2026-09-01: 0 of 24483 windows are clamped, and this keeps it that way.
win <- win[GenomicRanges::width(win) == 2L * half_width]
n_edge_clamped <- n_pre - length(win)
if (n_edge_clamped > 0)
  message("Dropped ", n_edge_clamped, " promoter window(s) clamped at a chromosome ",
          "edge (unequal widths break the cached-subset equivalence).")

# --- 2. Entrez -> SYMBOL (the join key to RNA-seq) ---------------------------
prom_entrez <- names(win)
if (is.null(prom_entrez)) {
  prom_entrez <- as.character(S4Vectors::mcols(win)$gene_id)
  names(win)  <- prom_entrez     # rownames of the cached matrix are these ids
}
sym <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_entrez, column = "SYMBOL",
                             keytype = "ENTREZID", multiVals = "first")
n_no_symbol <- sum(is.na(sym))
keep_sym <- !is.na(sym)
win <- win[keep_sym]; sym <- unname(sym[keep_sym])
# Zero collisions measured; assert rather than assume, because a symbol claimed by
# two Entrez promoters would double-count that gene in every stratum.
if (anyDuplicated(sym))
  stop("Entrez -> SYMBOL is not 1:1 (", sum(duplicated(sym)),
       " duplicate symbols). The expression join would double-count those genes.")
message(sprintf("Promoter windows: %d -> %d with a SYMBOL (%d unmapped, %d edge-clamped).",
                n_pre, length(win), n_no_symbol, n_edge_clamped))

# --- 3. Expression universe: rlog matrix + the silent set --------------------
expr <- local({
  rld <- readRDS("cache/rld.rds")
  m   <- SummarizedExperiment::assay(rld)
  rm(rld); gc(verbose = FALSE)
  m
})
si <- readRDS("cache/sample_info.rds")
gm <- readRDS("cache/gene_meta.rds")
stopifnot(identical(rownames(expr), gm$gene_id), !anyDuplicated(gm$gene_name))
rownames(expr) <- gm$gene_name

# "silent" = in the raw count table but below rnaseq.min_count summed over all 12
# libraries. Read through load_rnaseq_counts() with min_count = 0 rather than a
# bespoke read_tsv(), so the sample-map and duplicate-name handling cannot drift
# from rule 01's. cache/dds.rds cannot substitute: helpers.R applies the filter
# BEFORE caching, so dds holds only the survivors.
silent_syms <- character(0)
if (incl_silent && isTRUE(mp$expr_silent_from_counts %||% TRUE)) {
  if (!file.exists(rc$count_table)) {
    warning("rnaseq.count_table not readable ('", rc$count_table, "'); the silent ",
            "group is DROPPED. It is deliberately NOT widened to all promoters ",
            "missing from the DESeq2 set: ~40% of those are symbol-vocabulary ",
            "mismatches whose expression is unknown, not zero.", call. = FALSE)
  } else {
    rl0 <- load_rnaseq_counts(rc$count_table, rc$sample_map, genotypes,
                              ref_genotype = ref_geno, min_count = 0)
    below <- rowSums(rl0$counts, na.rm = TRUE) < (rc$min_count %||% 10)
    silent_syms <- setdiff(unique(rl0$gene_meta$gene_name[below]), gm$gene_name)
    rm(rl0); gc(verbose = FALSE)
  }
}
n_quantiled <- sum(sym %in% rownames(expr))
n_silent    <- sum(sym %in% silent_syms)
n_unknown   <- length(sym) - n_quantiled - n_silent
message(sprintf(paste0("Expression partition of %d promoters: %d quantiled, ",
                       "%d silent (<%d reads over all libraries), %d unknown (symbol ",
                       "absent from the RNA-seq vocabulary -> dropped, not silent)."),
                length(sym), n_quantiled, n_silent, rc$min_count %||% 10, n_unknown))

# --- 4. Expression group per promoter, per genotype --------------------------
# Cuts come from `basis`: "matched" uses only genes that HAVE a promoter, so the
# four quartiles carry equal n and their SEM ribbons are comparable; "all" cuts on
# the full DESeq2 set and yields unequal panels.
group_for <- function(g) {
  v <- rowMeans(expr[, si$condition == g, drop = FALSE])
  basis <- if (identical(basis_mode, "matched")) v[names(v) %in% sym] else v
  br <- stats::quantile(basis, seq(0, 1, length.out = n_q + 1L), na.rm = TRUE)
  br[1] <- -Inf; br[length(br)] <- Inf          # nothing may fall outside the cuts
  q <- cut(v, breaks = br, labels = q_levels, include.lowest = TRUE)
  out <- rep(NA_character_, length(sym))
  hit <- match(sym, names(v))
  out[!is.na(hit)] <- as.character(q[hit[!is.na(hit)]])
  if (incl_silent) out[is.na(out) & sym %in% silent_syms] <- "silent"
  factor(out, levels = c(if (incl_silent) "silent", q_levels))
}
rlog_for <- function(g) unname(rowMeans(expr[, si$condition == g, drop = FALSE])[match(sym, rownames(expr))])

groups_geno  <- setNames(lapply(genotypes, group_for), genotypes)
groups_fixed <- group_for(ref_geno)
arms <- list(per_genotype = if (per_geno) groups_geno
                            else setNames(rep(list(groups_fixed), length(genotypes)), genotypes))
if (fixed_arm)
  arms$fixed <- setNames(rep(list(groups_fixed), length(genotypes)), genotypes)
arm_labels <- c(per_genotype = if (per_geno) "cuts per genotype" else sprintf("cuts from %s", ref_geno),
                fixed = sprintf("cuts from %s, applied to all", ref_geno))

# Per-gene table: the join key rule 51 uses, plus every group assignment.
gene_groups <- data.frame(entrez = names(win), symbol = sym,
                          in_deseq2 = sym %in% rownames(expr),
                          in_silent = sym %in% silent_syms,
                          stringsAsFactors = FALSE)
for (g in genotypes) {
  gene_groups[[paste0("rlog_", g)]]  <- rlog_for(g)
  gene_groups[[paste0("group_", g)]] <- as.character(groups_geno[[g]])
}
gene_groups$group_fixed <- as.character(groups_fixed)
readr::write_csv(gene_groups, "results/tables/tss_expression_gene_groups.csv")

# Counts, including the two dropped classes, so the table sums to the promoter
# universe rather than quietly losing 3773 genes (the lesson from rules 31/32).
counts <- dplyr::bind_rows(lapply(names(arms), function(a)
  dplyr::bind_rows(lapply(genotypes, function(g) {
    tb <- table(arms[[a]][[g]], useNA = "no")
    data.frame(arm = a, genotype = g, expr_group = names(tb),
               n_promoters = as.integer(tb), stringsAsFactors = FALSE)
  }))))
counts <- dplyr::bind_rows(counts, data.frame(
  arm = "dropped", genotype = NA_character_,
  expr_group = c("no_symbol", "not_in_rnaseq", "edge_clamped"),
  n_promoters = c(n_no_symbol, n_unknown, n_edge_clamped), stringsAsFactors = FALSE))
readr::write_csv(counts, "results/tables/tss_expression_counts.csv")

# How much does quartile membership actually move between genotypes? This is the
# confound of per-genotype cuts, quantified rather than assumed.
migration <- dplyr::bind_rows(lapply(setdiff(genotypes, ref_geno), function(g) {
  tb <- as.data.frame(table(ref = groups_geno[[ref_geno]], ko = groups_geno[[g]]),
                      stringsAsFactors = FALSE)
  tb$ko_genotype <- g; tb
})) %>% dplyr::rename(n = Freq)
readr::write_csv(migration, "results/tables/tss_expression_quartile_migration.csv")

# --- 5. THE bigWig pass: one (promoter x bin) matrix per genotype ------------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)
bw_by_geno <- split(bw_meta$filepath, as.character(bw_meta$genotype))
missing_bw <- setdiff(genotypes, names(bw_by_geno))
if (length(missing_bw) > 0) stop("No ", cfg$assay, " bigWigs for: ", paste(missing_bw, collapse = ", "))

cache_ok <- function(x)
  is.list(x) && identical(x$half_width, half_width) && identical(x$n_bins, n_bins) &&
  identical(x$assay, cfg$assay) && identical(x$regions, win) &&
  setequal(names(x$rep_mean), genotypes)

cached <- NULL
if (file.exists(cache_file) && !isTRUE(mp$expr_rebuild_cache %||% FALSE)) {
  cached <- readRDS(cache_file)
  if (!cache_ok(cached)) {
    message("Cache ", cache_file, " does not match the current promoter set / window / ",
            "assay -- rebuilding. (Delete it manually if this repeats.)")
    cached <- NULL
  } else {
    message("Reusing cached profile matrices from ", cache_file, " (no bigWig I/O).")
  }
}

if (is.null(cached)) {
  positions <- seq(-half_width, half_width, length.out = n_bins)
  rep_means <- list()
  for (g in genotypes) {
    bws <- bw_by_geno[[g]]
    message("Profiling ", g, " (", length(bws), " bigWigs) over ", length(win), " promoters ...")
    rm_g <- rep_mean_matrix(win, bws, n_bins = n_bins, half_width = half_width)
    # Prove the cached path equals mean_replicate_profile() while the RleLists are
    # still hot: this is the whole correctness claim of the 50/51 split.
    if (isTRUE(mp$expr_verify_equivalence %||% TRUE)) {
      chk <- sample(seq_len(nrow(rm_g)), min(400L, nrow(rm_g)))
      ref <- mean_replicate_profile(win[chk], bws, n_bins = n_bins, half_width = half_width)
      new <- profile_from_rep_mean(rm_g, chk, positions = positions)
      stopifnot(isTRUE(all.equal(ref$mean, new$mean, tolerance = 0)),
                isTRUE(all.equal(ref$sem,  new$sem,  tolerance = 0)),
                identical(ref$position, new$position))
      message("  equivalence check passed (n = ", length(chk), ", tolerance = 0)")
    }
    rep_means[[g]] <- rm_g
    evict_bigwigs(bws)          # load-bearing: 12 whole-genome RleLists will not fit
  }
  cached <- list(rep_mean = rep_means, regions = win, positions = positions,
                 assay = cfg$assay, half_width = half_width, n_bins = n_bins,
                 bigwigs = lapply(bw_by_geno, basename), built = Sys.time(),
                 built_by = "scripts/50_tss_expression_metaprofiles.R")
  # Atomic write: an interrupted saveRDS would otherwise leave a truncated cache
  # that could still pass the metadata check on the next run.
  tmp <- paste0(cache_file, ".tmp")
  saveRDS(cached, tmp); file.rename(tmp, cache_file)
  message("Wrote ", cache_file)
}
positions <- cached$positions

# --- 6. Profiles: every stratum is a row subset ------------------------------
profiles <- dplyr::bind_rows(lapply(names(arms), function(a)
  dplyr::bind_rows(lapply(genotypes, function(g) {
    grp <- arms[[a]][[g]]
    dplyr::bind_rows(lapply(levels(grp), function(lv) {
      idx <- which(!is.na(grp) & grp == lv)
      p <- profile_from_rep_mean(cached$rep_mean[[g]], idx, positions = positions)
      if (is.null(p)) return(NULL)
      p$arm <- a; p$genotype <- g; p$expr_group <- lv
      p
    }))
  }))))
if (nrow(profiles) == 0) stop("No expression-stratified profiles computed.")

grp_levels <- c(if (incl_silent) "silent", q_levels)
profiles$expr_group <- factor(profiles$expr_group, levels = grp_levels)
profiles$genotype   <- factor(profiles$genotype, levels = genotypes)
profiles <- profiles %>% dplyr::rename(n_promoters = n_regions)
# Shape = the curve divided by its own window mean: isolates WHERE the signal sits
# from HOW MUCH there is. Expression correlates with open chromatin, so the raw
# amplitude ordering is expected and is not the finding.
profiles <- profiles %>%
  dplyr::group_by(arm, genotype, expr_group) %>%
  dplyr::mutate(window_mean = mean(mean, na.rm = TRUE),
                shape     = mean / window_mean,
                shape_sem = sem  / window_mean) %>%
  dplyr::ungroup()
readr::write_csv(profiles, "results/tables/tss_expression_profiles.csv")

# --- 7. Shape statistics (the quantitative answer) --------------------------
shape_stats <- profiles %>%
  dplyr::group_by(arm, genotype, expr_group) %>%
  dplyr::group_modify(~ cbind(data.frame(n_promoters = .x$n_promoters[1]),
                              profile_shape_stats(.x))) %>%
  dplyr::ungroup()
readr::write_csv(shape_stats, "results/tables/tss_expression_shape_stats.csv")

# --- 8. Figures --------------------------------------------------------------
sub_base <- sprintf("+/- %d bp, %d bins (bin = %.1f bp); %s",
                    half_width, n_bins, (2 * half_width) / (n_bins - 1),
                    if (identical(basis_mode, "matched"))
                      "quartile cuts from the promoter-matched genes (equal n per quartile)"
                    else "quartile cuts from all DESeq2 genes (unequal n per quartile)")
sem_note <- paste("ribbon = SEM across promoters (promoter-to-promoter spread),",
                  "NOT a test of the quartile or genotype difference")
lty_vals <- c(setNames(rep("solid", n_q), q_levels),
              silent = as.character(mp$expr_silent_linetype %||% "22"))

build_expr_fig <- function(df, ycol, semcol, fname, title, subtitle, height = 7) {
  # panel_n_labels() is hard-wired to a two-level up/down direction column and
  # stops on duplicate keys (helpers.R:98-106), so the n annotation is built here.
  n_lab <- df %>%
    dplyr::distinct(genotype, expr_group, n_promoters) %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarise(label = paste(sprintf("%s %s", expr_group,
                                           format(n_promoters, big.mark = ",", trim = TRUE)),
                                   collapse = " | "), .groups = "drop")
  p <- ggplot(df, aes(position, .data[[ycol]], colour = expr_group, fill = expr_group)) +
    geom_ribbon(aes(ymin = .data[[ycol]] - .data[[semcol]],
                    ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = expr_group), linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.03, vjust = 1.4, size = 2.4, colour = "grey25") +
    facet_wrap(~ genotype, scales = "free_y", nrow = 2) +
    scale_colour_viridis_d(end = 0.92, name = "Expression", drop = FALSE) +
    scale_fill_viridis_d(end = 0.92, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = lty_vals, guide = "none", drop = FALSE) +
    labs(x = "Distance to TSS (bp)", y = NULL, title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, fig_subdir, fig_root, width = 11, height = height)
}

prim <- profiles %>% dplyr::filter(arm == "per_genotype")
build_expr_fig(prim, "mean", "sem", "tss_signal_by_expression",
  "Promoter G4 signal at the TSS, by RNA-seq expression quartile",
  sprintf("%s; %s\nfree y per genotype panel: the 5 curves WITHIN a panel are directly comparable, heights BETWEEN panels are not; %s",
          sub_base, arm_labels[["per_genotype"]], sem_note))

if (isTRUE(mp$expr_shape_normalise %||% TRUE))
  build_expr_fig(prim, "shape", "shape_sem", "tss_signal_by_expression_shape",
    "Promoter G4 signal SHAPE at the TSS, by RNA-seq expression quartile",
    sprintf("%s; each curve divided by its own window mean, so only POSITION differs\nread differences off tss_expression_shape_stats.csv: below one bin (%.1f bp) is not resolvable",
            sub_base, (2 * half_width) / (n_bins - 1)))

if (fixed_arm)
  build_expr_fig(profiles %>% dplyr::filter(arm == "fixed"), "mean", "sem",
    "tss_signal_by_expression_fixedcuts",
    "Promoter G4 signal at the TSS, by expression quartile (fixed cuts)",
    sprintf("%s; %s\nsensitivity arm: gene membership is IDENTICAL across panels, so any panel difference here is signal, not relabelling",
            sub_base, arm_labels[["fixed"]]))

p_n <- ggplot(counts %>% dplyr::filter(arm == "per_genotype"),
              aes(expr_group, n_promoters, fill = expr_group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = format(n_promoters, big.mark = ",", trim = TRUE)),
            vjust = -0.3, size = 2.6) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_viridis_d(end = 0.92, guide = "none", drop = FALSE) +
  labs(x = NULL, y = "Promoters (n)",
       title = "Promoters per expression group",
       subtitle = sprintf(paste("Excluded and reported separately: %d without a SYMBOL,",
                                "%d whose symbol is absent from the RNA-seq vocabulary",
                                "(expression unknown, NOT silent)"), n_no_symbol, n_unknown)) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_n, "tss_expression_group_sizes", fig_subdir, fig_root, width = 11, height = 4.5)

# --- 9. Report ---------------------------------------------------------------
message("\nPromoters per (arm, genotype, expression group):")
print(as.data.frame(counts))
message("\nCurve shape statistics (this is where a positional difference is read off):")
print(as.data.frame(shape_stats %>% dplyr::filter(arm == "per_genotype")))
if (per_geno) {
  stable <- migration %>% dplyr::filter(ref == ko) %>% dplyr::group_by(ko_genotype) %>%
    dplyr::summarise(n_stable = sum(n), .groups = "drop")
  tot <- sum(!is.na(groups_geno[[ref_geno]]))
  message("\nQuartile membership stability vs ", ref_geno,
          " (per-genotype cuts move genes between panels):")
  print(as.data.frame(stable %>% dplyr::mutate(pct = round(100 * n_stable / tot, 2))))
}
message("\nDone. Rule 51 reads ", cache_file, " and opens no bigWigs.")
