#!/usr/bin/env Rscript
# ============================================================================
# 02_deseq2.R   (system R)   --- Stage A ---
#
# Differential expression for all 5 contrasts (ashr shrinkage), from the shared
# dds cached by script 01. Writes per-contrast full tables (analysis_V2 format)
# AND the canonical long deseq2_KO_vs_WT.csv consumed unchanged by the DeepG4
# downstream rules (29, 06, 38, 31, 32, 46, 47) — schema matches the original
# DeepG4 13a_deseq2.R exactly: gene_id, gene_name, ko, baseMean, log2FC, padj.
#
# Inputs:   cache/dds.rds, cache/gene_meta.rds
# Outputs:  cache/deseq2_res_list.rds
#           results/tables/deseq2_results_<contrast>.tsv  (x5)
#           results/tables/deseq2_KO_vs_WT.csv            (3 KO-vs-WT, long)
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/tables")

dds       <- readRDS("cache/dds.rds")
gene_meta <- readRDS("cache/gene_meta.rds")
use_ashr  <- isTRUE(cfg$rnaseq$lfc_shrink) && requireNamespace("ashr", quietly = TRUE)

# config contrasts: list of [numerator, denominator]
contrasts <- lapply(cfg$rnaseq$contrasts, function(x) as.character(unlist(x)))
names(contrasts) <- vapply(contrasts, function(x) paste0(x[1], "_vs_", x[2]), character(1))

gm <- gene_meta[, c("gene_id", "gene_name")]
res_list <- lapply(names(contrasts), function(nm) {
  pr <- contrasts[[nm]]
  ct <- c("genotype", pr[1], pr[2])
  res <- DESeq2::results(dds, contrast = ct, alpha = cfg$rnaseq$padj)
  if (use_ashr) res <- DESeq2::lfcShrink(dds, contrast = ct, res = res, type = "ashr")
  df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                   log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                   pvalue = res$pvalue, padj = res$padj, stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, gm, by = "gene_id")
  df[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
})
names(res_list) <- names(contrasts)
saveRDS(res_list, "cache/deseq2_res_list.rds")

# --- per-contrast full tables ----------------------------------------------
for (nm in names(res_list)) {
  df <- res_list[[nm]] %>% dplyr::arrange(padj)
  readr::write_tsv(df, file.path("results/tables", paste0("deseq2_results_", nm, ".tsv")))
}

# --- canonical long KO_vs_WT.csv (DeepG4 13a schema) ------------------------
ref  <- cfg$ref_genotype
ko_contrasts <- Filter(function(x) x[2] == ref, contrasts)
long <- dplyr::bind_rows(lapply(names(ko_contrasts), function(nm) {
  ko <- ko_contrasts[[nm]][1]
  res_list[[nm]] %>%
    transmute(gene_id, gene_name, ko = ko, baseMean, log2FC = log2FoldChange, padj)
}))
readr::write_csv(long, "results/tables/deseq2_KO_vs_WT.csv")

message("DESeq2 done. DEG counts (padj<", cfg$rnaseq$padj, ", |log2FC|>", cfg$rnaseq$lfc, "):")
for (nm in names(res_list)) {
  df <- res_list[[nm]]
  n_up <- sum(df$padj < cfg$rnaseq$padj & df$log2FoldChange >  cfg$rnaseq$lfc, na.rm = TRUE)
  n_dn <- sum(df$padj < cfg$rnaseq$padj & df$log2FoldChange < -cfg$rnaseq$lfc, na.rm = TRUE)
  message(sprintf("  %-16s up=%d down=%d", nm, n_up, n_dn))
}
