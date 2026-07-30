#!/usr/bin/env Rscript
# ============================================================================
# 12_peak_region_enrichment.R   (system R)   --- Stage D ---
#
# Peak-to-region priority annotation, per-feature peak-count bars (stacked +
# fraction), and analytical binomial fold-enrichment vs genome background.
# Ported from extended_V3.Rmd chunks s3-peak-region-annot / -bar / s3-perm-test
# / s3-perm-plot (binomial replaces regioneR permutation, which hangs on Windows).
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
genome_bp <- sum(as.numeric(chrom_sizes))
perm <- dplyr::bind_rows(lapply(names(peak_sets), function(a) {
  dplyr::bind_rows(lapply(names(peak_sets[[a]]), function(g) {
    pk <- peak_sets[[a]][[g]]; if (length(pk) == 0) return(NULL)
    n_peaks <- length(pk)
    dplyr::bind_rows(lapply(names(region_priority), function(rn) {
      rg <- region_priority[[rn]]
      p_expect <- min(sum(as.numeric(GenomicRanges::width(rg))) / genome_bp, 1)
      obs <- sum(GenomicRanges::countOverlaps(pk, rg) > 0)
      data.frame(assay = a, genotype = g, region = rn, obs = obs,
                 exp_mean = n_peaks * p_expect, fold = obs / max(n_peaks * p_expect, 1),
                 p_perm = binom.test(obs, n_peaks, p_expect, alternative = "greater")$p.value)
    }))
  }))
}))
perm$assay_label  <- unname(assay_label_map[perm$assay])
perm$region_label <- factor(perm$region, levels = region_levels, labels = region_labels)
perm$genotype     <- factor(perm$genotype, levels = genotypes)
readr::write_tsv(perm, "results/tables/region_fold_enrichment.tsv")

p_fold <- ggplot(perm, aes(region_label, log2(fold), fill = genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = ifelse(p_perm < 0.001, "***", ifelse(p_perm < 0.01, "**",
                        ifelse(p_perm < 0.05, "*", "ns")))),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = condition_colours) +
  facet_wrap(~ assay_label, scales = "free_y") +
  labs(x = NULL, y = "log2(observed / expected)",
       title = "Region enrichment relative to random placement",
       subtitle = "Binomial null: expected = n_peaks x (region_bp / genome_bp)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_fold, "02_region_fold_enrichment", "06_regional", fig_root, width = 9, height = 5)
message("Peak-region enrichment done.")
