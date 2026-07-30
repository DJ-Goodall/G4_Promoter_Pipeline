#!/usr/bin/env Rscript
# ============================================================================
# 10_peak_qc.R   (system R)   --- Stage B ---
#
# Peak QC figures: peak-width distribution (per assay x genotype) and region
# distribution by chromosome. Ported from 20260522_g4_gloop_extended_V3.Rmd
# chunks s2-peak-widths / s2-region-chr-bar.
#
# Inputs:   cache/peaks_{assay}_{genotype}.rds (12 per-genotype sets),
#           cache/regions_all.rds
# Outputs:  results/figures/04_peaks/{01_peak_width_distribution,
#           01_regions_by_chromosome}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(S4Vectors); library(ggplot2)
  library(dplyr); library(scales); library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures")
fig_root <- "results/figures"

std_chroms <- cfg$std_chroms
assays     <- as.character(unlist(cfg$peak_calling$assays))
genotypes  <- as.character(unlist(cfg$peak_calling$genotypes))

# --- Peak-width distribution ------------------------------------------------
combos <- expand.grid(assay = assays, genotype = genotypes, stringsAsFactors = FALSE)
width_df <- dplyr::bind_rows(lapply(seq_len(nrow(combos)), function(i) {
  a <- combos$assay[i]; g <- combos$genotype[i]
  f <- file.path("cache", sprintf("peaks_%s_%s.rds", a, g))
  if (!file.exists(f)) return(NULL)
  pk <- readRDS(f); if (length(pk) == 0) return(NULL)
  data.frame(assay = a, genotype = g, width = GenomicRanges::width(pk))
}))
width_df$genotype <- factor(width_df$genotype, levels = genotypes)

p_width <- ggplot(width_df, aes(x = width, fill = genotype)) +
  geom_histogram(bins = 50, colour = "white", linewidth = 0.1, alpha = 0.85) +
  scale_fill_manual(values = condition_colours) +
  scale_x_log10(labels = scales::comma) +
  facet_grid(assay ~ genotype, labeller = labeller(assay = assay_label_map)) +
  labs(title = "Peak-width distribution", x = "Peak width (bp, log10)", y = "Peak count") +
  theme_pub() + theme(legend.position = "none")
save_plot(p_width, "01_peak_width_distribution", "04_peaks", fig_root, width = 9, height = 5)

# --- Region distribution by chromosome --------------------------------------
regions_all <- readRDS(file.path("cache", "regions_all.rds"))
region_palette <- setNames(RColorBrewer::brewer.pal(length(regions_all), "Set2"),
                           names(regions_all))
chr_tbl <- dplyr::bind_rows(lapply(names(regions_all), function(nm) {
  data.frame(region = nm, chrom = as.character(GenomicRanges::seqnames(regions_all[[nm]])))
})) %>% dplyr::filter(chrom %in% std_chroms) %>% dplyr::count(region, chrom)

p_chr <- ggplot(chr_tbl, aes(x = factor(chrom, levels = std_chroms), y = n, fill = region)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = region_palette) +
  labs(x = "Chromosome", y = "Region count", title = "Region distribution by chromosome") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_chr, "01_regions_by_chromosome", "04_peaks", fig_root, width = 10, height = 4)

message("Peak QC figures written to results/figures/04_peaks/")
