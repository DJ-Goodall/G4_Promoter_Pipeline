#!/usr/bin/env Rscript
# ============================================================================
# 52_gene_body_expression_metaprofiles.R   (system R)   --- Stage N, point 1 ---
#
# Rule 50's question over the whole transcription unit instead of the promoter:
# do genes with LOW vs HIGH expression carry their G4 and R-loop CUT&Tag signal
# at different positions across the GENE? One panel per genotype, stratified by
# the same RNA-seq expression groups rule 50 used.
#
# ---------------------------------------------------------------------------
# v2 (2026-09-15): THE AXIS CHANGED, because v1's was misleading.
#
# v1 used a pure scale-regions layout: 50 bp flank bins either side of a body
# compressed to 60 bins. That put a 6.5x resolution STEP at exactly the TSS
# (50 bp vs a median 323 bp body bin), and the step -- not biology -- drew the
# vertical cliff that appeared at the TSS in every v1 figure. Measured on the
# SAME genes and the SAME bigWigs at a flat 50 bp (rule 50):
#
#     bp from TSS   -25    +25   +125   +275   +425   +530
#     G4 (WT Q4)   4.59   4.53   4.17   3.58   3.14   2.95
#
# G4 signal carries ~500 bp PAST the TSS and decays smoothly. v1 reported that
# entire decay as two numbers, 3.67 and 2.96, and drew it as a cliff.
#
# Worse, a v1 body bin was the gene's own length / 60, a DIFFERENT genomic
# distance for every gene: 32 bp at the 5th percentile, 323 bp median, 3,352 bp
# at the 95th, 597 kb at the maximum. 37% of genes had a first body bin wider
# than the 500 bp decay and 21% wider than 1 kb, so for those the
# promoter-proximal shoulder was averaged into intronic background BEFORE the
# group mean was taken. Median mm10 5'UTR is 161 bp and 68.9% of genes had a
# first body bin wider than that: the 5'UTR had no dedicated bin anywhere.
#
# v2 scales ONLY the deep interior:
#
#   [ TSS-2kb .. TSS+2kb ]   [ interior, scaled ]   [ TES-2kb .. TES+2kb ]
#      80 bins of 50 bp           60 bins             80 bins of 50 bp
#          REAL bp                 NOT bp                 REAL bp
#
# so the 5'UTR shoulder and the TES region are linear in bp and the TSS block
# sits on the same grid as rule 50's figure. Only the stretch between TSS+2kb
# and TES-2kb is scaled, and nothing there is positionally stereotyped anyway.
#
# PRICE, measured and NOT uniform: a gene needs >= 2*2000 + 1000 = 5000 bp for
# the two anchor windows not to overlap, so 21,946 -> 18,219 genes (83.0%).
# Dropped per group: silent 1,273 (24.2%), Q1 603 (16.1%), Q2 421 (11.2%),
# Q3 314 (8.3%), Q4 420 (11.2%). Silent genes are shorter, so that curve is
# reweighted hardest -- the attrition is written into the counts table and
# named in every figure subtitle rather than left to shift the curves quietly.
# Short genes are not lost to the project: rule 50 profiles them at the TSS.
#
# WHY THE INTERIOR IS FLAT, and why that is not a bug (measured 2026-09-15):
# 52.2% of the 62,913 G4 peaks and 40.9% of the 42,960 R-loop peaks fall inside
# gene bodies, but gene bodies are 43% of the genome -- so body peaks are only
# ~1.2x enriched against ~2.5x at the TSS, and they sit at idiosyncratic
# positions per gene. Averaged over ~18,000 genes they do not align into a bump;
# they ARE the raised interior floor (G4 ~2.5 vs a silent-gene baseline ~1.9).
# "Flat interior" means "not positionally stereotyped", NOT "no G4s".
# ---------------------------------------------------------------------------
#
# TWO ASSAYS, ONE GENE SET. metaprofile.expr_gene_assays lists them (G4_BG4 +
# Rloop_S96). Both are profiled over the IDENTICAL genes, bins and expression
# groups, so a G4 panel and an R-loop panel for the same genotype and quartile
# are built from exactly the same genes -- only the bigWigs differ. One cache
# per assay; rule 53 reads only the G4_BG4 one, because topology is a property
# of a G4 peak and an R-loop topology figure would be meaningless.
#
# EXPRESSION GROUPS ARE NOT RECOMPUTED. They are read from rule 50's
# tss_expression_gene_groups.csv, so a gene is in the same quartile in both
# figures by construction and the TSS-window and whole-gene views are directly
# comparable.
#
# GENE BODIES COME FROM TxDb, NOT FROM THE PROMOTER GRanges. The Entrez join is
# a clean 1:1 here, unlike rule 50's SYMBOL bridge. Note that genes() returns
# the OUTERMOST span over all isoforms, so the TES is the most distal poly(A)
# site of any transcript -- one reason a TES feature would be blurred.
#
# THE LENGTH-TERTILE FIGURE (gene_signal_by_length*) exists because the interior
# is still scaled: an interior bin is 50 bp for a 8 kb gene and 8 kb for a
# 500 kb one. Splitting on length makes that confound visible rather than
# assumed, and separates it from expression, which correlates with length
# (median: silent 20.6 kb, Q1 25.9, Q2 28.7, Q3 31.2, Q4 27.4 kb).
#
# THE GENOTYPE OVERLAY (gene_signal_by_genotype*) is the one figure here that
# compares AMPLITUDE between genotypes, and nothing in the metaprofile rules
# depth-normalises a library -- so it ships as a raw figure carrying that caveat
# plus a shape companion in which each genotype is divided by its own window
# mean and depth cancels. Read the shape one for genotype effects.
#
# A bigWig's coverage can end before the BSgenome chromosome does (measured:
# chr16 stops 199 kb short). Those bins are NA, not clamped -- clamping would
# rescale that gene's interior against the others.
#
# Inputs:   results/tables/tss_expression_gene_groups.csv  (rule 50),
#           G4 bigWigs, R-loop bigWigs
# Outputs:  cache/gene_composite_repmean_<assay>_w<W>_ab<AB>_bb<BB>.rds
#             (one per assay; the G4_BG4 one is rule 53's input)
#           results/tables/gene_body_expression_genes.csv
#           results/tables/gene_body_expression_counts.csv
#           results/tables/gene_body_expression_profiles.csv
#           results/tables/gene_body_expression_shape_stats.csv
#           results/tables/gene_body_genotype_profiles.csv
#           results/tables/gene_body_genotype_shape_stats.csv
#           results/tables/gene_body_length_profiles.csv
#           results/figures/25_tss_tes_metaprofiles/
#             gene_signal_by_expression{,_rloop}{,_shape,_fixedcuts}.{pdf,png}
#             gene_signal_by_genotype{,_rloop}{,_shape}.{pdf,png}
#             gene_signal_by_length{,_rloop}.{pdf,png}
#             gene_expression_group_sizes.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(GenomicFeatures)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

fig_subdir <- "25_tss_tes_metaprofiles"
fig_root   <- "results/figures"
ensure_dirs("cache", "results/tables", file.path(fig_root, fig_subdir))

mp        <- cfg$metaprofile
genotypes <- cfg$genotypes
ref_geno  <- cfg$ref_genotype
W         <- as.integer(mp$expr_gene_anchor_window %||% 2000L)
abin      <- as.integer(mp$expr_gene_anchor_bin    %||% 50L)
n_bbins   <- as.integer(mp$expr_gene_body_bins     %||% 60L)
min_int   <- as.integer(mp$expr_gene_min_interior  %||% 1000L)
min_len   <- as.integer(mp$expr_gene_min_length    %||% 5000L)
max_len   <- if (is.null(mp$expr_gene_max_length)) Inf else as.numeric(mp$expr_gene_max_length)
n_tert    <- as.integer(mp$expr_gene_length_strata %||% 3L)
n_q       <- as.integer(mp$expr_n_quantiles %||% 4)
per_geno  <- isTRUE(mp$expr_per_genotype %||% TRUE)
fixed_arm <- isTRUE(mp$expr_fixed_cuts_arm %||% TRUE) && per_geno
incl_sil  <- isTRUE(mp$expr_include_silent %||% TRUE)
q_levels  <- paste0("Q", seq_len(n_q))
grp_levels <- c(if (incl_sil) "silent", q_levels)
n_abins   <- as.integer(2L * W / abin)
n_bins    <- 2L * n_abins + n_bbins

# The layout floor is derived, never assumed: if someone lowers min_length
# without lowering the anchor window, composite_bin_ranges() would throw deep
# inside the bigWig pass after minutes of I/O. Fail here instead.
if (min_len < 2L * W + min_int)
  stop("metaprofile.expr_gene_min_length (", min_len, ") is below ",
       "2 * expr_gene_anchor_window + expr_gene_min_interior (",
       2L * W + min_int, "). The two anchor windows would overlap.")

cache_file <- function(a)
  sprintf("cache/gene_composite_repmean_%s_w%d_ab%d_bb%d.rds", a, W, abin, n_bbins)

# --- 1. Expression groups from rule 50 --------------------------------------
# entrez forced to character: readr would otherwise coerce the ids to numeric
# and every match() against TxDb names would half-succeed (helpers contract).
gg_path <- "results/tables/tss_expression_gene_groups.csv"
if (!file.exists(gg_path))
  stop("Missing ", gg_path, ". Run rule tss_expression_metaprofiles (50) first -- ",
       "this rule deliberately reuses its quartile assignment rather than ",
       "recomputing it, so the TSS and whole-gene figures stay comparable.")
gg <- readr::read_csv(gg_path,
        col_types = readr::cols(entrez = readr::col_character(),
                                .default = readr::col_guess()))
n_universe <- nrow(gg)
message("Rule 50 gene groups: ", n_universe, " promoters.")

# --- 2. Gene bodies ----------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_all <- suppressMessages(
  GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE))
GenomeInfoDb::seqlevelsStyle(genes_all) <- "UCSC"
genes_all <- GenomeInfoDb::keepSeqlevels(
  genes_all, intersect(GenomeInfoDb::seqlevels(genes_all), cfg$std_chroms),
  pruning.mode = "coarse")

ix <- match(gg$entrez, names(genes_all))
n_no_gene <- sum(is.na(ix))
gg <- gg[!is.na(ix), ]; genes <- genes_all[ix[!is.na(ix)]]
names(genes) <- gg$entrez
stopifnot(identical(names(genes), gg$entrez))

gl <- GenomicRanges::width(genes)
keep_len <- gl >= min_len & gl <= max_len
n_short <- sum(gl < min_len); n_long <- sum(gl > max_len)

# Chromosome-edge feasibility. The anchor windows reach W bp beyond each end,
# and bin_matrix_from_ranges() refuses a bin that falls off a chromosome rather
# than clamping it (clamping would silently rescale that gene against others).
chrom_sizes <- default_chrom_sizes(cfg$std_chroms)
GenomeInfoDb::seqlengths(genes) <- chrom_sizes[GenomeInfoDb::seqlevels(genes)]
keep_edge <- (GenomicRanges::start(genes) - W) >= 1L &
             (GenomicRanges::end(genes) + W) <=
               chrom_sizes[as.character(GenomicRanges::seqnames(genes))]
n_edge <- sum(!keep_edge)

# The per-group attrition of the length floor, recorded BEFORE the filter. This
# is the one real cost of the composite axis and it is not uniform (silent genes
# are shorter, so that curve is reweighted hardest); leaving it implicit would
# let a composition change read as a signal change.
drop_by_group <- as.data.frame(table(
  group = factor(gg$group_fixed, levels = grp_levels),
  dropped = !(keep_len & keep_edge)))
attrition <- drop_by_group %>%
  dplyr::filter(dropped == TRUE) %>%
  dplyr::transmute(arm = "dropped_by_length_floor", genotype = NA_character_,
                   expr_group = as.character(group), n_genes = as.integer(Freq))

keep <- keep_len & keep_edge
genes <- genes[keep]; gg <- gg[keep, ]
message(sprintf(paste0("Gene bodies: %d promoters -> %d profiled ",
                       "(%d with no single-strand TxDb gene, %d shorter than %d bp, ",
                       "%d longer than the cap, %d within %d bp of a chromosome end)."),
                nrow(gg) + n_no_gene + sum(!keep), length(genes), n_no_gene,
                n_short, min_len, n_long, n_edge, W))
if (!length(genes)) stop("No genes survived the length / edge filters.")

bl  <- GenomicRanges::width(genes)
ilen <- bl - 2L * W                       # the stretch that is actually scaled
message(sprintf("Interior bin width (bp): min %d, median %d, 95th pct %d, max %d",
                round(min(ilen) / n_bbins), round(stats::median(ilen) / n_bbins),
                round(stats::quantile(ilen, 0.95) / n_bbins), round(max(ilen) / n_bbins)))
message(sprintf("Anchor bins are a fixed %d bp; %d bins per anchor block.", abin, n_abins))

# Length tertiles, cut on the profiled set so every stratum is non-empty.
tert_br <- stats::quantile(bl, probs = seq(0, 1, length.out = n_tert + 1L))
tert_br[1] <- -Inf; tert_br[length(tert_br)] <- Inf
tert_lab <- vapply(seq_len(n_tert), function(i) {
  lo <- stats::quantile(bl, (i - 1) / n_tert); hi <- stats::quantile(bl, i / n_tert)
  sprintf("%s (%.0f-%.0f kb)", c("short", "medium", "long", paste0("L", seq_len(n_tert)))[i],
          lo / 1000, hi / 1000)
}, character(1))
length_tertile <- cut(bl, breaks = tert_br, labels = tert_lab, include.lowest = TRUE)

gene_tbl <- data.frame(entrez = names(genes), symbol = gg$symbol,
  chrom = as.character(GenomicRanges::seqnames(genes)),
  start = GenomicRanges::start(genes), end = GenomicRanges::end(genes),
  strand = as.character(GenomicRanges::strand(genes)),
  length_bp = bl, interior_bp = ilen, interior_bin_bp = round(ilen / n_bbins),
  length_tertile = as.character(length_tertile), stringsAsFactors = FALSE)
for (g in genotypes) {
  cn <- paste0("group_", g)
  gene_tbl[[cn]] <- if (cn %in% names(gg)) as.character(gg[[cn]]) else NA_character_
}
gene_tbl$group_fixed <- as.character(gg$group_fixed)
readr::write_csv(gene_tbl, "results/tables/gene_body_expression_genes.csv")

# --- 3. Expression arms (same convention as rule 50) ------------------------
group_col <- function(g) {
  cn <- paste0("group_", g)
  as.character(if (per_geno && cn %in% names(gg)) gg[[cn]] else gg$group_fixed)
}
arms <- list(per_genotype = setNames(lapply(genotypes, group_col), genotypes))
if (fixed_arm)
  arms$fixed <- setNames(rep(list(as.character(gg$group_fixed)), length(genotypes)),
                         genotypes)
arm_labels <- c(per_genotype = if (per_geno) "cuts per genotype"
                               else sprintf("cuts from %s", ref_geno),
                fixed = sprintf("cuts from %s, applied to all", ref_geno))

# "unknown_expression" is an explicit ROW, not an implicit NA: rule 50's third
# class is a promoter whose symbol is absent from the RNA-seq vocabulary, so its
# expression is UNKNOWN rather than zero, and folding it into "silent" would make
# a no-expression claim about genes that were never measured.
counts <- dplyr::bind_rows(lapply(names(arms), function(a)
  dplyr::bind_rows(lapply(genotypes, function(g) {
    v  <- arms[[a]][[g]]
    tb <- table(factor(v, levels = grp_levels))
    data.frame(arm = a, genotype = g,
               expr_group = c(names(tb), "unknown_expression"),
               n_genes = c(as.integer(tb), sum(is.na(v))), stringsAsFactors = FALSE)
  }))))
counts <- dplyr::bind_rows(counts, data.frame(
  arm = "dropped", genotype = NA_character_,
  expr_group = c("no_txdb_gene", "shorter_than_min", "longer_than_max", "chrom_edge"),
  n_genes = c(n_no_gene, n_short, n_long, n_edge), stringsAsFactors = FALSE), attrition)
# The table must account for every promoter rule 50 handed us. An assertion, not
# a comment, because a silently unbalanced denominator is the failure class the
# v1.1 audit was about. `attrition` is a BREAKDOWN of rows already counted in
# `dropped`, so it is excluded from the closure sum.
per_geno_total <- sum(counts$n_genes[counts$arm == "per_genotype"]) / length(genotypes)
dropped_total  <- sum(counts$n_genes[counts$arm == "dropped"])
stopifnot(per_geno_total == length(genes))
if (per_geno_total + dropped_total != n_universe)
  stop("Counts do not close: ", per_geno_total, " profiled + ", dropped_total,
       " dropped = ", per_geno_total + dropped_total, ", but rule 50 supplied ",
       n_universe, " promoters.")
message(sprintf("Counts close: %d profiled + %d dropped = %d (rule 50 universe).",
                per_geno_total, dropped_total, n_universe))
readr::write_csv(counts, "results/tables/gene_body_expression_counts.csv")

attr_txt <- paste(sprintf("%s -%s", attrition$expr_group,
                          format(attrition$n_genes, big.mark = ",", trim = TRUE)),
                  collapse = ", ")
message("Dropped by the ", min_len, " bp length floor, per group: ", attr_txt)

# --- 4. THE bigWig pass: one (gene x bin) matrix per genotype, per assay ----
assay_keys <- as.character(unlist(mp$expr_gene_assays %||% cfg$assay))
if (!cfg$assay %in% assay_keys)
  stop("metaprofile.expr_gene_assays (", paste(assay_keys, collapse = ", "),
       ") does not include the primary assay ", cfg$assay,
       ", but rule 53 reads that assay's cache. Add it back.")

# Short, filename-safe token per assay. The G4 figures keep the names they had,
# so nothing downstream or on an existing poster has to be renamed.
assay_tag <- function(a) if (identical(a, cfg$assay)) "" else
  paste0("_", tolower(gsub("[^A-Za-z0-9]+", "", sub("_.*$", "", a))))
# %||% cannot help here: assay_label_map[<unknown>] returns a named NA, not
# NULL, so an unmapped assay would be printed as "NA" in every title.
assay_name <- function(a)
  if (a %in% names(assay_label_map)) unname(assay_label_map[a]) else a

ramp_for <- function(a) {
  anchors <- if (identical(a, cfg$assay))
    mp$expr_palette %||% c("#2166AC", "#4393C3", "#D6604D", "#B2182B")
  else
    mp$expr_rloop_palette %||% c("#542788", "#998EC3", "#F1A340", "#B35806")
  c(setNames(grDevices::colorRampPalette(as.character(unlist(anchors)))(n_q), q_levels),
    silent = as.character(mp$expr_silent_colour %||% "grey55"))
}

bw_for <- function(a) {
  bm <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = a)
  by_g <- split(bm$filepath, as.character(bm$genotype))
  missing_bw <- setdiff(genotypes, names(by_g))
  if (length(missing_bw) > 0)
    stop("No ", a, " bigWigs for: ", paste(missing_bw, collapse = ", "),
         ". Check paths.bigwig_dir and metaprofile.expr_gene_assays.")
  by_g
}

# Cache key encodes EVERY parameter that changes the result -- assay, anchor
# window, anchor bin, interior bins -- plus the region set itself. The v1 caches
# (cache/gene_expr_repmean_*) are a DIFFERENT layout and a different name, so
# they are orphaned rather than silently reused; delete them.
cache_ok <- function(x, a)
  is.list(x) && identical(x$anchor_window, W) && identical(x$anchor_bin, abin) &&
  identical(x$n_body_bins, n_bbins) && identical(x$min_interior, min_int) &&
  identical(x$assay, a) && identical(x$regions, genes) &&
  setequal(names(x$rep_mean), genotypes)

load_or_build <- function(a) {
  cf <- cache_file(a)
  cached <- NULL
  if (file.exists(cf) && !isTRUE(mp$expr_gene_rebuild_cache %||% FALSE)) {
    cached <- readRDS(cf)
    if (!cache_ok(cached, a)) {
      message("Cache ", cf, " does not match the current gene set / layout / ",
              "assay -- rebuilding. (Delete it manually if this repeats.)")
      cached <- NULL
    } else {
      message("Reusing cached ", a, " profile matrices from ", cf, " (no bigWig I/O).")
    }
  }
  if (!is.null(cached)) return(cached)

  by_g <- bw_for(a)
  rep_means <- list()
  for (g in genotypes) {
    bws <- by_g[[g]]
    message("Profiling ", a, " / ", g, " (", length(bws), " bigWigs) over ",
            length(genes), " genes, ", n_bins, " bins each ...")
    rep_means[[g]] <- composite_rep_mean_matrix(genes, bws, anchor_window = W,
                        anchor_bin = abin, n_body_bins = n_bbins,
                        min_interior = min_int)
    evict_bigwigs(bws)      # load-bearing: 12 whole-genome RleLists will not fit
  }
  cached <- list(rep_mean = rep_means, regions = genes,
                 positions = composite_positions(W, abin, n_bbins),
                 assay = a, anchor_window = W, anchor_bin = abin,
                 n_body_bins = n_bbins, min_interior = min_int,
                 bigwigs = lapply(by_g, basename), built = Sys.time(),
                 built_by = "scripts/52_gene_body_expression_metaprofiles.R")
  # Atomic write: an interrupted saveRDS would otherwise leave a truncated cache
  # that could still pass the metadata check on the next run.
  tmp <- paste0(cf, ".tmp")
  saveRDS(cached, tmp); file.rename(tmp, cf)
  message("Wrote ", cf)
  cached
}

caches <- setNames(lapply(assay_keys, load_or_build), assay_keys)
layout <- composite_positions(W, abin, n_bbins)
positions <- layout$x
axis_spec <- composite_axis(W)

# --- 5. Profiles: every stratum is a row subset ------------------------------
add_shape <- function(df, keys)
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::mutate(window_mean = mean(mean, na.rm = TRUE),
                  shape     = mean / window_mean,
                  shape_sem = sem  / window_mean) %>%
    dplyr::ungroup()
# Every profile carries the block it sits in, so a consumer can tell a real-bp
# bin from a scaled one without re-deriving the layout.
# profile_from_rep_mean() calls its axis column `position`, but on a composite
# axis that value is PLOT SPACE, not a genomic position -- the real distance is
# `bp`, and it is NA across the interior. Renaming here keeps the two apart for
# anyone reading the CSV.
with_layout <- function(df)
  dplyr::left_join(dplyr::rename(df, x = position),
                   layout[, c("bin", "block", "bp", "frac")], by = "bin")

strata_profiles <- function(split_vec_by_geno, label_levels, key) {
  dplyr::bind_rows(lapply(assay_keys, function(a)
    dplyr::bind_rows(lapply(genotypes, function(g) {
      v <- split_vec_by_geno[[g]]
      dplyr::bind_rows(lapply(label_levels, function(lv) {
        idx <- which(!is.na(v) & v == lv)
        p <- profile_from_rep_mean(caches[[a]]$rep_mean[[g]], idx, positions = positions)
        if (is.null(p)) return(NULL)
        p$assay <- a; p$genotype <- g; p[[key]] <- lv
        p
      }))
    }))))
}

profiles <- dplyr::bind_rows(lapply(names(arms), function(arm) {
  d <- strata_profiles(arms[[arm]], grp_levels, "expr_group"); d$arm <- arm; d
}))
if (nrow(profiles) == 0) stop("No expression-stratified gene profiles computed.")
profiles$expr_group <- factor(profiles$expr_group, levels = grp_levels)
profiles$genotype   <- factor(profiles$genotype, levels = genotypes)
profiles$assay      <- factor(profiles$assay, levels = assay_keys)
profiles <- profiles %>% dplyr::rename(n_genes = n_regions)
# Shape = curve / its own window mean: isolates WHERE the signal sits from HOW
# MUCH there is. Expression correlates with open chromatin, so the raw amplitude
# ordering is expected and is not the finding.
profiles <- add_shape(profiles, c("assay", "arm", "genotype", "expr_group")) %>%
  with_layout() %>% dplyr::relocate(assay)
readr::write_csv(profiles, "results/tables/gene_body_expression_profiles.csv")

# Length tertiles: the same matrix cut on gene length instead of expression.
tert_by_geno <- setNames(rep(list(as.character(length_tertile)), length(genotypes)),
                         genotypes)
lenprof <- strata_profiles(tert_by_geno, tert_lab, "length_tertile")
lenprof$length_tertile <- factor(lenprof$length_tertile, levels = tert_lab)
lenprof$genotype <- factor(lenprof$genotype, levels = genotypes)
lenprof$assay    <- factor(lenprof$assay, levels = assay_keys)
lenprof <- lenprof %>% dplyr::rename(n_genes = n_regions) %>%
  add_shape(c("assay", "genotype", "length_tertile")) %>%
  with_layout() %>% dplyr::relocate(assay)
readr::write_csv(lenprof, "results/tables/gene_body_length_profiles.csv")

# The genotype overlay: NO expression split at all, every profiled gene in one
# curve per genotype. rows = NULL takes the whole matrix, which is the honest
# "overall" -- it includes the unknown-expression genes the quartile figures
# cannot place.
overlay <- dplyr::bind_rows(lapply(assay_keys, function(a)
  dplyr::bind_rows(lapply(genotypes, function(g) {
    p <- profile_from_rep_mean(caches[[a]]$rep_mean[[g]], rows = NULL,
                               positions = positions)
    p$assay <- a; p$genotype <- g; p
  }))))
overlay$genotype <- factor(overlay$genotype, levels = genotypes)
overlay$assay    <- factor(overlay$assay, levels = assay_keys)
overlay <- overlay %>% dplyr::rename(n_genes = n_regions) %>%
  add_shape(c("assay", "genotype")) %>% with_layout() %>% dplyr::relocate(assay)
readr::write_csv(overlay, "results/tables/gene_body_genotype_profiles.csv")

# --- 6. Shape statistics ----------------------------------------------------
# These are the numbers the v1 axis could not produce. Because both anchor
# blocks are real bp, tss_half_decay_bp answers "how far into the gene does the
# promoter G4 shoulder actually carry?" in bp -- the question the cliff hid --
# and tes_peak_bp says whether anything is happening at the terminator.
# The frac_* columns are bin-weighted, not bp-weighted (interior bins cover more
# genome than anchor bins), so read them as "share of the drawn curve".
composite_stats <- function(d) {
  d <- d[order(d$bin), ]
  tss_b <- d[d$block == "tss", ]; tes_b <- d[d$block == "tes", ]
  int_b <- d[d$block == "interior", ]
  floor_v <- mean(int_b$mean, na.rm = TRUE)
  pk_i  <- which.max(tss_b$mean)
  pk_bp <- if (length(pk_i)) tss_b$bp[pk_i] else NA_real_
  pk_v  <- if (length(pk_i)) tss_b$mean[pk_i] else NA_real_
  half  <- floor_v + (pk_v - floor_v) / 2
  dn    <- tss_b[!is.na(tss_b$bp) & tss_b$bp > pk_bp, ]
  hd    <- if (nrow(dn) && any(dn$mean <= half, na.rm = TRUE))
             dn$bp[which(dn$mean <= half)[1]] else NA_real_
  m <- pmax(d$mean, 0); m[!is.finite(m)] <- 0; tot <- sum(m)
  tes_i <- which.max(tes_b$mean)
  data.frame(
    tss_peak_bp       = pk_bp,
    tss_peak_value    = pk_v,
    tss_half_decay_bp = hd,          # bp past the TSS peak at which the shoulder
                                     # falls halfway to the interior floor
    interior_mean     = floor_v,
    tes_peak_bp       = if (length(tes_i)) tes_b$bp[tes_i] else NA_real_,
    tes_peak_value    = if (length(tes_i)) tes_b$mean[tes_i] else NA_real_,
    tes_over_interior = if (length(tes_i)) tes_b$mean[tes_i] / floor_v else NA_real_,
    window_mean       = mean(d$mean, na.rm = TRUE),
    frac_tss_block    = if (tot > 0) sum(m[d$block == "tss"]) / tot else NA_real_,
    frac_interior     = if (tot > 0) sum(m[d$block == "interior"]) / tot else NA_real_,
    frac_tes_block    = if (tot > 0) sum(m[d$block == "tes"]) / tot else NA_real_)
}
stats_by <- function(df, keys)
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) %>%
    dplyr::group_modify(~ cbind(data.frame(n_genes = .x$n_genes[1]),
                                composite_stats(.x))) %>%
    dplyr::ungroup()

shape_stats <- stats_by(profiles, c("assay", "arm", "genotype", "expr_group"))
readr::write_csv(shape_stats, "results/tables/gene_body_expression_shape_stats.csv")
overlay_stats <- stats_by(overlay, c("assay", "genotype"))
readr::write_csv(overlay_stats, "results/tables/gene_body_genotype_shape_stats.csv")

# --- 7. Figures --------------------------------------------------------------
sub_base <- sprintf(paste0("anchors +/-%s kb in %d bp bins (REAL bp); interior ",
                           "scaled to %d bins (median %d bp per interior bin, range %d-%d)"),
                    format(W / 1000, trim = TRUE), abin, n_bbins,
                    round(stats::median(ilen) / n_bbins),
                    round(min(ilen) / n_bbins), round(max(ilen) / n_bbins))
scale_note <- paste0("x is linear in bp WITHIN each anchor block; only the shaded ",
                     "interior is scaled, and no distance may be read off it")
floor_note <- sprintf(paste("%s genes (>= %s bp, so the two anchor windows do not",
                            "overlap); the length floor drops %s -- not evenly, so",
                            "compare groups with that in mind"),
                      format(length(genes), big.mark = ",", trim = TRUE),
                      format(min_len, big.mark = ",", trim = TRUE), attr_txt)
sem_note <- paste("ribbon = SEM across genes (gene-to-gene spread),",
                  "NOT a test of the quartile or genotype difference")

lty_vals <- c(setNames(rep("solid", n_q), q_levels),
              silent = as.character(mp$expr_silent_linetype %||% "22"))

# The composite axis, shared by every figure in this rule (and by rule 53) so
# they cannot drift on where the TSS and TES are drawn. The shaded band IS the
# scaled stretch: the warning is visual, not only written in the subtitle.
gene_axis_layers <- function()
  list(annotate("rect", xmin = axis_spec$interior[1], xmax = axis_spec$interior[2],
                ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.45),
       geom_vline(xintercept = axis_spec$dividers, colour = "grey60", linewidth = 0.3),
       geom_vline(xintercept = axis_spec$anchors, linetype = "dashed", colour = "grey30"),
       scale_x_continuous(breaks = axis_spec$breaks, labels = axis_spec$labels,
                          limits = axis_spec$limits, expand = ggplot2::expansion(mult = 0.004)))

build_fig <- function(df, ycol, semcol, colour_col, cols, legend_title,
                      fname, title, subtitle, use_lty = FALSE, height = 7) {
  n_lab <- df %>%
    dplyr::distinct(genotype, .data[[colour_col]], n_genes) %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarise(label = paste(sprintf("%s %s", .data[[colour_col]],
                                    format(n_genes, big.mark = ",", trim = TRUE)),
                                   collapse = " | "), .groups = "drop")
  p <- ggplot(df, aes(x, .data[[ycol]],
                      colour = .data[[colour_col]], fill = .data[[colour_col]])) +
    gene_axis_layers() +
    geom_ribbon(aes(ymin = .data[[ycol]] - .data[[semcol]],
                    ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.15, colour = NA) +
    {if (use_lty) geom_line(aes(linetype = .data[[colour_col]]), linewidth = 0.7)
     else geom_line(linewidth = 0.7)} +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.02, vjust = 1.4, size = 2.3, colour = "grey25") +
    facet_wrap(~ genotype, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    scale_fill_manual(values = cols, guide = "none", drop = FALSE) +
    {if (use_lty) scale_linetype_manual(values = lty_vals, guide = "none", drop = FALSE)} +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom",
                        axis.text.x = element_text(size = 7))
  save_plot(p, fname, fig_subdir, fig_root, width = 12, height = height)
}

for (a in assay_keys) {
  tag <- assay_tag(a); nm <- assay_name(a); cols <- ramp_for(a)
  pa   <- profiles %>% dplyr::filter(assay == a)
  prim <- pa %>% dplyr::filter(arm == "per_genotype")

  build_fig(prim, "mean", "sem", "expr_group", cols, "Expression",
    paste0("gene_signal_by_expression", tag),
    sprintf("%s signal across the gene, TSS to TES, by RNA-seq expression quartile", nm),
    sprintf("%s; %s\n%s\n%s\nfree y per genotype panel: curves WITHIN a panel are comparable, heights BETWEEN panels are not; %s",
            sub_base, arm_labels[["per_genotype"]], scale_note, floor_note, sem_note),
    use_lty = TRUE)

  if (isTRUE(mp$expr_shape_normalise %||% TRUE))
    build_fig(prim, "shape", "shape_sem", "expr_group", cols, "Expression",
      paste0("gene_signal_by_expression", tag, "_shape"),
      sprintf("%s signal SHAPE across the gene, TSS to TES, by expression quartile", nm),
      sprintf("%s; each curve divided by its own window mean, so only POSITION differs\n%s\n%s\nread tss_half_decay_bp and tes_over_interior off gene_body_expression_shape_stats.csv",
              sub_base, scale_note, floor_note),
      use_lty = TRUE)

  if (fixed_arm)
    build_fig(pa %>% dplyr::filter(arm == "fixed"), "mean", "sem", "expr_group",
      cols, "Expression", paste0("gene_signal_by_expression", tag, "_fixedcuts"),
      sprintf("%s signal across the gene by expression quartile (fixed cuts)", nm),
      sprintf("%s; %s\n%s\n%s\nsensitivity arm: gene membership is IDENTICAL across panels, so any panel difference here is signal, not relabelling",
              sub_base, arm_labels[["fixed"]], scale_note, floor_note),
      use_lty = TRUE)

  # Length tertiles. Diagnostic, so a deliberately neutral grey ramp: this is
  # about whether the interior scaling is distorting the comparison, not biology.
  len_cols <- setNames(as.character(unlist(
    mp$expr_length_palette %||% c("#BDBDBD", "#737373", "#252525")))[seq_len(n_tert)],
    tert_lab)
  build_fig(lenprof %>% dplyr::filter(assay == a), "mean", "sem", "length_tertile",
    len_cols, "Gene length", paste0("gene_signal_by_length", tag),
    sprintf("%s signal across the gene, TSS to TES, by GENE LENGTH tertile", nm),
    sprintf("%s\n%s\nDIAGNOSTIC: the anchor blocks are real bp and should be near-identical across tertiles; any tertile difference in the INTERIOR is the scaling, since an interior bin is %d bp for the shortest genes and %d bp for the longest\n%s",
            sub_base, scale_note, round(min(ilen) / n_bbins),
            round(max(ilen) / n_bbins), floor_note))
}

# --- 7b. Genotype overlay: all genes, one panel, four curves -----------------
# The only figure in this rule that compares AMPLITUDE between genotypes, which
# is exactly where raw bigWig values are least defensible: nothing in the
# metaprofile rules scales a library by its depth or its background. Rules 39/40
# do normalise (cnt_norm: background_median) but that step re-reads every bigWig
# uncached and would multiply this rule's runtime. Both figures are drawn: the
# raw one carries the caveat, and the shape companion divides each genotype by
# its own window mean so POSITION can be compared without depending on depth.
geno_cols <- setNames(
  as.character(unlist(mp$expr_genotype_palette %||%
    c("#333333", "#0072B2", "#D55E00", "#CC79A7")))[seq_along(genotypes)],
  genotypes)
if (anyNA(geno_cols))
  stop("metaprofile.expr_genotype_palette has fewer colours (",
       length(unlist(mp$expr_genotype_palette)), ") than genotypes (",
       length(genotypes), ").")

build_overlay_fig <- function(df, ycol, semcol, fname, title, subtitle) {
  n_lab <- data.frame(label = sprintf("n = %s genes per curve",
                        format(max(df$n_genes), big.mark = ",", trim = TRUE)))
  p <- ggplot(df, aes(x, .data[[ycol]], colour = genotype, fill = genotype)) +
    gene_axis_layers() +
    geom_ribbon(aes(ymin = .data[[ycol]] - .data[[semcol]],
                    ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.85) +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.04, vjust = 1.6, size = 2.8, colour = "grey25") +
    scale_colour_manual(values = geno_cols, name = NULL, drop = FALSE) +
    scale_fill_manual(values = geno_cols, guide = "none", drop = FALSE) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom",
                        axis.text.x = element_text(size = 7))
  save_plot(p, fname, fig_subdir, fig_root, width = 9, height = 5.5)
}

if (isTRUE(mp$expr_gene_overlay %||% TRUE)) {
  n_all <- format(length(genes), big.mark = ",", trim = TRUE)
  for (a in assay_keys) {
    tag <- assay_tag(a); nm <- assay_name(a)
    oa  <- overlay %>% dplyr::filter(assay == a)
    build_overlay_fig(oa, "mean", "sem", paste0("gene_signal_by_genotype", tag),
      sprintf("%s signal across the gene, TSS to TES, by genotype", nm),
      sprintf("all %s profiled genes, NO expression split; %s\n%s\nCURVE HEIGHTS ARE RAW bigWig VALUES, NOT DEPTH-NORMALISED BETWEEN GENOTYPES -- for a genotype comparison read the _shape companion; %s",
              n_all, sub_base, scale_note, sem_note))
    build_overlay_fig(oa, "shape", "shape_sem",
      paste0("gene_signal_by_genotype", tag, "_shape"),
      sprintf("%s signal SHAPE across the gene, TSS to TES, by genotype", nm),
      sprintf("all %s profiled genes; each genotype divided by its OWN window mean, so only POSITION differs and library depth cancels\n%s",
              n_all, scale_note))
  }
}

# unknown_expression is excluded from the bars (it has no place on the expression
# ramp) but named in the subtitle, so it is visible rather than gone. One figure,
# not one per assay: the gene set and its expression groups are assay-independent.
n_unknown <- counts$n_genes[counts$arm == "per_genotype" &
                            counts$expr_group == "unknown_expression"][1]
p_n <- ggplot(counts %>% dplyr::filter(arm == "per_genotype",
                                       expr_group != "unknown_expression"),
              aes(expr_group, n_genes, fill = expr_group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = format(n_genes, big.mark = ",", trim = TRUE)),
            vjust = -0.3, size = 2.6) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_manual(values = ramp_for(cfg$assay), guide = "none", drop = FALSE) +
  labs(x = NULL, y = "Genes (n)",
       title = "Genes per expression group (composite TSS-TES window)",
       subtitle = sprintf(paste("Also profiled but not shown: %d genes of UNKNOWN",
                                "expression (symbol absent from the RNA-seq",
                                "vocabulary -- not silent).\nExcluded: %d with no",
                                "single-strand TxDb gene, %d shorter than %d bp,",
                                "%d longer than the cap, %d within %d bp of a",
                                "chromosome end.\nLength floor removed, per group: %s"),
                          n_unknown, n_no_gene, n_short, min_len, n_long, n_edge,
                          W, attr_txt)) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_n, "gene_expression_group_sizes", fig_subdir, fig_root,
          width = 11, height = 4.5)

# --- 8. Report ---------------------------------------------------------------
message("\nAssays profiled: ", paste(assay_keys, collapse = ", "))
message("\nGenes per (arm, genotype, expression group):")
print(as.data.frame(counts))
message("\nComposite shape statistics (per_genotype arm).")
message("tss_half_decay_bp = bp past the TSS peak at which the shoulder falls ",
        "halfway to the interior floor -- the number the v1 axis could not give:")
print(as.data.frame(shape_stats %>% dplyr::filter(arm == "per_genotype") %>%
        dplyr::select(assay, genotype, expr_group, n_genes, tss_peak_bp,
                      tss_half_decay_bp, interior_mean, tes_peak_bp, tes_over_interior)))
message("\nGenotype overlay (all genes, no expression split):")
print(as.data.frame(overlay_stats))
message("\nDone. Rule 53 reads ", cache_file(cfg$assay), " and opens no bigWigs.")
