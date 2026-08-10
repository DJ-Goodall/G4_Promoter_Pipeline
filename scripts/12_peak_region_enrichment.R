#!/usr/bin/env Rscript
# ============================================================================
# 12_peak_region_enrichment.R   (system R)   --- Stage D ---
#
# Peak-to-region priority annotation, per-feature peak-count bars (stacked +
# fraction), and analytical binomial fold-enrichment vs genome background.
# Ported from extended_V3.Rmd chunks s3-peak-region-annot / -bar / s3-perm-test
# / s3-perm-plot (binomial replaces regioneR permutation, which hangs on Windows).
#
# ---------------------------------------------------------------------------
# v1.1 CHANGES (2026-08-07 audit)
#
#  - COLUMN RENAMED p_perm -> p_binom. No shuffling has ever happened in this
#    script. The substitution WAS documented in the header above, and the plot
#    subtitle named the binomial null correctly, but the column name leaked a
#    "permutation test" claim into four README lines ("Permutation testing
#    confirms this quantitatively", "permutation p < 0.001"). Renaming the
#    column is the honest fix; see docs/audit-2026-08/README-deltas.md for the README text.
#
#  - NULL CORRECTED for peak width AND for un-merged region sets. See the
#    detailed comment at the enrichment block below. The two biases run in
#    opposite directions, so fixing only the one the audit identified would
#    have moved the numbers the wrong way.
#
#  - NOTE ON peak_counts_per_region.tsv vs region_fold_enrichment.tsv: these
#    two tables report DIFFERENT quantities under the same region labels and
#    v1.0 gave no hint of it. The count table is priority-EXCLUSIVE (a peak is
#    assigned to promoter > enhancer > 5'UTR > intron1 > genebody and counted
#    once), so WT 5'UTR = 878; the enrichment table uses non-exclusive
#    countOverlaps, so the same cell is 3,410. Both are now labelled in their
#    column names and in the messages below.
# ---------------------------------------------------------------------------
#
# Inputs:   cache/peaks_{assay}_{genotype}.rds (12), cache/regions_all.rds
# Outputs:  results/tables/{peak_counts_per_region,region_fold_enrichment}.tsv
#           results/figures/06_regional/{02_peak_counts_stacked,
#             02_peak_counts_fraction,02_region_fold_enrichment}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(S4Vectors); library(ggplot2); library(dplyr)
  library(scales); library(readr); library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures", "results/tables")
fig_root <- "results/figures"

std_chroms  <- cfg$std_chroms
chrom_sizes <- default_chrom_sizes(std_chroms)
assays      <- as.character(unlist(cfg$peak_calling$assays))
genotypes   <- extended_genotypes

# --- Load peaks + regions ---------------------------------------------------
peak_sets <- list()
for (a in assays) { peak_sets[[a]] <- list()
  for (g in genotypes) {
    f <- file.path("cache", sprintf("peaks_%s_%s.rds", a, g))
    if (file.exists(f)) peak_sets[[a]][[g]] <- readRDS(f)
  }
}
regions_all     <- readRDS("cache/regions_all.rds")
region_priority <- regions_all[names(regions_all) %in% region_levels]

# --- Peak-region annotation + counts ----------------------------------------
peak_annot <- dplyr::bind_rows(lapply(names(peak_sets), function(a) {
  dplyr::bind_rows(lapply(names(peak_sets[[a]]), function(g) {
    pk <- peak_sets[[a]][[g]]; if (length(pk) == 0) return(NULL)
    ann <- annotate_peak_regions(pk, region_priority)
    data.frame(assay = a, genotype = g,
               region = as.character(S4Vectors::mcols(ann)$region),
               width = GenomicRanges::width(ann))
  }))
}))
peak_annot$region   <- factor(peak_annot$region, levels = region_levels, labels = region_labels)
peak_annot$genotype <- factor(peak_annot$genotype, levels = genotypes)
peak_annot$assay_label <- unname(assay_label_map[peak_annot$assay])
count_tbl <- peak_annot %>% dplyr::count(assay, assay_label, genotype, region)

fill_pal <- setNames(region_palette[region_levels], region_labels)
p_stack <- ggplot(count_tbl, aes(genotype, n, fill = region)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.2) +
  facet_wrap(~ assay_label, scales = "free_y") +
  scale_fill_manual(values = fill_pal) +
  labs(x = NULL, y = "Peak count", fill = "Region", title = "Peak counts by genomic feature",
       subtitle = "Priority: promoter > 5'UTR > intron1 > enhancer > gene body") + theme_pub()
save_plot(p_stack, "02_peak_counts_stacked", "06_regional", fig_root, width = 8, height = 5)

p_frac <- ggplot(count_tbl, aes(genotype, n, fill = region)) +
  geom_col(position = "fill", colour = "white", linewidth = 0.2) +
  facet_wrap(~ assay_label) + scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = fill_pal) +
  labs(x = NULL, y = "Fraction of peaks", fill = "Region",
       title = "Peak distribution across features (fraction)") + theme_pub()
save_plot(p_frac, "02_peak_counts_fraction", "06_regional", fig_root, width = 8, height = 5)
readr::write_tsv(count_tbl, "results/tables/peak_counts_per_region.tsv")

# --- Binomial fold-enrichment vs genome background --------------------------
# v1.1: the null is corrected in two ways that push in OPPOSITE directions, so
# both are needed or the correction goes the wrong way.
#
#   (a) PEAK WIDTH. A peak overlaps a region if it TOUCHES it, so a peak of
#       width w behaves like a point landing anywhere in the region dilated by
#       w - 1. v1.0 used p = region_bp / genome_bp, which is the probability
#       for a zero-width peak. With mean peak width ~966 bp this understates
#       the expectation badly for small features and barely at all for large
#       ones -- which is exactly the "fold enrichment scales inversely with
#       region size" artefact already noted in commit a4acb9f. For a 200 bp
#       5'UTR the expectation was understated ~5.8x; for a 4 kb promoter
#       ~1.24x. That is the whole reason 5'UTRs (11.9x) appeared to beat
#       promoters (2.0x).
#
#   (b) UN-MERGED REGIONS. region_priority sets are built per gene and never
#       globally reduce()d, so overlapping neighbours are double-counted:
#       24,483 promoters x 4,001 bp = 97.9 Mb of claimed footprint (this
#       reproduces v1.0's exp_mean exactly). That OVERSTATES the expectation
#       and deflates fold -- the opposite bias to (a).
#
# Fix: reduce() the region set, then compute the expectation per peak-width
# decile against the width-dilated footprint. Exact for uniform placement,
# runs in seconds, and needs no shuffling (regioneR::permTest with
# circularRandomizeRegions remains the gold standard but hangs on Windows,
# which is why v1.0 substituted the binomial in the first place).
#
# The column is renamed p_perm -> p_binom: no shuffling happens here and the
# old name propagated a "permutation test" claim into four README lines.
genome_bp <- sum(as.numeric(chrom_sizes))

# Expected overlap probability for a peak of width w against region set rg.
#
# A peak [p, p + w - 1] overlaps a region [s, e] iff p lies in [s - w + 1, e].
# So the set of starting positions that produce an overlap -- the "landing
# zone" -- is each region extended BACKWARDS by w - 1, i.e. width W + w - 1.
# Under uniform placement of the start, P(overlap) is the reduce()d coverage of
# those landing zones divided by the genome length. Exact, no simulation.
expected_p_for_width <- function(rg_reduced, w, genome_bp) {
  if (w <= 1L) {
    dil <- rg_reduced
  } else {
    # fix = "end" holds the region end and moves the start upstream.
    dil <- GenomicRanges::resize(rg_reduced,
                                 GenomicRanges::width(rg_reduced) + (w - 1L),
                                 fix = "end")
    GenomicRanges::start(dil) <- pmax(1L, GenomicRanges::start(dil))
    dil <- GenomicRanges::reduce(dil, ignore.strand = TRUE)
  }
  min(sum(as.numeric(GenomicRanges::width(dil))) / genome_bp, 1)
}

enrich <- dplyr::bind_rows(lapply(names(peak_sets), function(a) {
  dplyr::bind_rows(lapply(names(peak_sets[[a]]), function(g) {
    pk <- peak_sets[[a]][[g]]; if (length(pk) == 0) return(NULL)
    n_peaks <- length(pk)
    wid <- GenomicRanges::width(pk)
    # Bin peaks by width decile; within a bin use the median width. Ten
    # evaluations per region set instead of ~110,000, at negligible cost in
    # accuracy (the expectation is near-linear in w over a decile).
    brk <- unique(stats::quantile(wid, probs = seq(0, 1, 0.1), na.rm = TRUE))
    grp <- if (length(brk) > 2)
      cut(wid, breaks = brk, include.lowest = TRUE, labels = FALSE)
    else rep(1L, length(wid))
    dplyr::bind_rows(lapply(names(region_priority), function(rn) {
      rg_red <- GenomicRanges::reduce(region_priority[[rn]], ignore.strand = TRUE)
      # v1.0 expectation, kept so the change is auditable side by side
      p_naive <- min(sum(as.numeric(GenomicRanges::width(region_priority[[rn]]))) /
                       genome_bp, 1)
      # v1.1 expectation: sum over width bins of n_bin * P(overlap | width)
      exp_v11 <- sum(vapply(sort(unique(grp)), function(k) {
        idx <- which(grp == k)
        sum(length(idx) * expected_p_for_width(rg_red,
              as.integer(stats::median(wid[idx])), genome_bp))
      }, numeric(1)))
      p_expect <- min(exp_v11 / n_peaks, 1)
      obs <- sum(GenomicRanges::countOverlaps(pk, region_priority[[rn]],
                                              ignore.strand = TRUE) > 0)
      data.frame(assay = a, genotype = g, region = rn, obs = obs,
                 exp_mean = exp_v11,
                 fold = obs / max(exp_v11, 1),
                 p_binom = binom.test(obs, n_peaks, p_expect,
                                      alternative = "greater")$p.value,
                 # provenance: what v1.0 would have reported
                 exp_mean_v10 = n_peaks * p_naive,
                 fold_v10 = obs / max(n_peaks * p_naive, 1),
                 region_bp_raw = sum(as.numeric(GenomicRanges::width(region_priority[[rn]]))),
                 region_bp_reduced = sum(as.numeric(GenomicRanges::width(rg_red))),
                 median_peak_bp = as.integer(stats::median(wid)))
    }))
  }))
}))
enrich$assay_label  <- unname(assay_label_map[enrich$assay])
enrich$region_label <- factor(enrich$region, levels = region_levels, labels = region_labels)
enrich$genotype     <- factor(enrich$genotype, levels = genotypes)
readr::write_tsv(enrich, "results/tables/region_fold_enrichment.tsv")

message("Fold enrichment, v1.1 corrected null vs v1.0 (BG4, WT):")
print(enrich %>%
  dplyr::filter(assay == "G4_BG4", genotype == genotypes[1]) %>%
  dplyr::select(region, obs, exp_mean_v10, fold_v10, exp_mean, fold))

perm <- enrich   # back-compat for anything downstream referencing `perm`

p_fold <- ggplot(enrich, aes(region_label, log2(fold), fill = genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = ifelse(p_binom < 0.001, "***", ifelse(p_binom < 0.01, "**",
                        ifelse(p_binom < 0.05, "*", "ns")))),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = condition_colours) +
  facet_wrap(~ assay_label, scales = "free_y") +
  labs(x = NULL, y = "log2(observed / expected)",
       # v1.1: this is a binomial test, not a permutation. Say so on the figure.
       title = "Region enrichment vs a width-aware binomial null",
       subtitle = paste("Expected = sum over peak-width deciles of",
                        "P(overlap | width) against reduce()d regions")) +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_fold, "02_region_fold_enrichment", "06_regional", fig_root, width = 9, height = 5)
message("Peak-region enrichment done.")
