#!/usr/bin/env Rscript
# ============================================================================
# 01_rnaseq_qc.R   (system R)   --- Stage A ---
#
# Load the RNA-seq count table, build the DESeqDataSet, run DESeq2 + rlog, and
# produce QC figures (library size, PCA, sample-correlation, gene detection).
# The dds / rld / sample_info / gene_meta objects are cached for scripts 02-04.
# Ported from 20260521_g4_gloop_analysis_V2.Rmd (chunks load-counts .. qc-*).
#
# Inputs:   config rnaseq.count_table (GSE269081_count_table.tsv)
# Outputs:  cache/{dds,rld,sample_info,gene_meta}.rds
#           results/figures/01_qc/{library_size_barplot,pca_plot,
#             sample_correlation_heatmap,gene_detection_barplot}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(cowplot); library(ggrepel)
  library(pheatmap); library(RColorBrewer); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures")
fig_root <- "results/figures"

rc <- cfg$rnaseq
rl <- load_rnaseq_counts(rc$count_table, rc$sample_map, cfg$genotypes,
                         ref_genotype = cfg$ref_genotype, min_count = rc$min_count)
counts    <- rl$counts                       # rownames = gene_id
coldata   <- rl$coldata                      # $genotype (WT first)
gene_meta <- rl$gene_meta                    # gene_id, gene_name, label
message(sprintf("Genes after filter: %d | samples: %d", nrow(counts), ncol(counts)))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                      design = ~ genotype)
dds <- DESeq2::DESeq(dds)
rld <- DESeq2::rlog(dds, blind = TRUE)

sample_info <- data.frame(sample_id = rownames(coldata),
                          condition = coldata$genotype,
                          row.names = rownames(coldata))

saveRDS(dds,         "cache/dds.rds")
saveRDS(rld,         "cache/rld.rds")
saveRDS(sample_info, "cache/sample_info.rds")
saveRDS(gene_meta,   "cache/gene_meta.rds")

# --- 2a. Library size --------------------------------------------------------
lib_df <- data.frame(sample = colnames(counts), reads = colSums(counts),
                     condition = sample_info$condition)
p_lib <- ggplot(lib_df, aes(reorder(sample, -reads), reads / 1e6, fill = condition)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = condition_colours) +
  labs(title = "Library size per sample", y = "Total counts (millions)", x = "") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_lib, "library_size_barplot", "01_qc", fig_root, width = 10, height = 6)

# --- 2b. PCA -----------------------------------------------------------------
pca_data <- DESeq2::plotPCA(rld, intgroup = "genotype", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = genotype)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
  scale_colour_manual(values = condition_colours) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA — rlog-transformed counts") + theme_cowplot()
save_plot(p_pca, "pca_plot", "01_qc", fig_root, width = 8, height = 6)

# --- 2c. Sample-to-sample correlation ---------------------------------------
cor_mat <- cor(SummarizedExperiment::assay(rld), method = "pearson")
ann_col <- data.frame(Condition = sample_info$condition, row.names = rownames(sample_info))
ann_colours <- list(Condition = condition_colours)
hm_fn <- function() pheatmap::pheatmap(cor_mat, annotation_col = ann_col,
  annotation_colors = ann_colours,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100),
  display_numbers = TRUE, number_format = "%.3f", fontsize_number = 7,
  main = "Sample-to-sample Pearson correlation (rlog)")
save_base_plot(hm_fn, "sample_correlation_heatmap", "01_qc", fig_root, width = 10, height = 8)

# --- 2d. Gene detection (>1 CPM) --------------------------------------------
cpm_mat  <- t(t(counts) / colSums(counts) * 1e6)
det_df   <- data.frame(sample = colnames(counts), n_genes = colSums(cpm_mat > 1),
                       condition = sample_info$condition)
p_det <- ggplot(det_df, aes(reorder(sample, -n_genes), n_genes, fill = condition)) +
  geom_col(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = condition_colours) +
  labs(title = "Genes detected per sample (>1 CPM)", y = "Number of genes", x = "") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_det, "gene_detection_barplot", "01_qc", fig_root, width = 10, height = 6)

message("RNA-seq QC done.")
