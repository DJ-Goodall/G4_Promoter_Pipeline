#!/usr/bin/env Rscript
# ============================================================================
# 04_enrichment.R   (system R)   --- Stage A ---
#
# Functional enrichment on the DESeq2 results: GO:BP (up/down per 3 vs-WT
# contrasts), GSEA against MSigDB Hallmark, and targeted C5 GO:BP pathway tests
# for dKO. Ported from analysis_V2.Rmd chunks gene-id-mapping / go-enrichment /
# gsea-hallmark / gsea-specific-dko. dplyr::select is qualified (org.Mm.eg.db S4
# `select` masks it).
#
# Inputs:   cache/deseq2_res_list.rds
# Outputs:  results/tables/go_enrichment_summary.tsv
#           results/figures/03_enrichment/{GO_dotplot_<c>_{up,down},
#             gsea_hallmark_<c>, gsea_specific_pathways_dKO}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(fgsea)
  library(msigdbr); library(ggplot2); library(cowplot); library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/figures", "results/tables")
fig_root <- "results/figures"
rc <- cfg$rnaseq

res_list <- readRDS("cache/deseq2_res_list.rds")
vs_wt    <- names(res_list)[grepl("_vs_WT$", names(res_list))][1:3]

# SYMBOL -> ENTREZ over all gene_names in the results
all_symbols <- unique(unlist(lapply(res_list, function(d) d$gene_name)))
gene_map <- clusterProfiler::bitr(all_symbols, fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Mm.eg.db)
bg_entrez <- gene_map$ENTREZID

# --- 4a. GO:BP enrichment (up/down) -----------------------------------------
go_results <- list()
for (nm in vs_wt) {
  df <- res_list[[nm]]
  sig <- df[!is.na(df$padj) & df$padj < rc$padj, ]
  up  <- sig$gene_name[sig$log2FoldChange >  rc$lfc]
  dn  <- sig$gene_name[sig$log2FoldChange < -rc$lfc]
  for (dir in c("up", "down")) {
    genes <- if (dir == "up") up else dn
    ent <- gene_map$ENTREZID[gene_map$SYMBOL %in% genes]
    if (length(ent) < 5) next
    ego <- clusterProfiler::enrichGO(gene = ent, universe = bg_entrez,
      OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = rc$go_pvalue, readable = TRUE)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    go_results[[paste0(nm, "_", dir)]] <- ego
    p <- clusterProfiler::dotplot(ego, showCategory = rc$go_show) +
      ggtitle(paste(nm, "—", dir, "GO:BP"))
    save_plot(p, paste0("GO_dotplot_", nm, "_", dir), "03_enrichment", fig_root,
              width = 10, height = 10)
  }
}
go_summary <- dplyr::bind_rows(lapply(names(go_results), function(n) {
  d <- as.data.frame(go_results[[n]]); if (nrow(d) == 0) return(NULL)
  d$contrast_direction <- n; d
}))
if (!is.null(go_summary) && nrow(go_summary) > 0)
  readr::write_tsv(go_summary, "results/tables/go_enrichment_summary.tsv")

# --- helper: ranked entrez stat vector for a contrast -----------------------
build_ranks <- function(df) {
  df <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]
  df$entrez <- gene_map$ENTREZID[match(df$gene_name, gene_map$SYMBOL)]
  df <- df[!is.na(df$entrez), ]
  r <- setNames(-log10(df$pvalue) * sign(df$log2FoldChange), df$entrez)
  r <- r[order(-abs(r))]; r <- r[!duplicated(names(r))]
  sort(r, decreasing = TRUE)
}

# --- 4b. GSEA Hallmark ------------------------------------------------------
hallmark <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- lapply(split(hallmark$entrez_gene, hallmark$gs_name), as.character)
for (nm in vs_wt) {
  ranks <- build_ranks(res_list[[nm]])
  fr <- fgsea::fgsea(pathways = hallmark_list, stats = ranks,
                     minSize = rc$gsea_minsize, maxSize = rc$gsea_maxsize)
  fr <- fr[order(fr$pval), ]
  top_pos <- head(fr[fr$NES > 0, ][order(fr[fr$NES > 0, ]$padj), ], 10)
  top_neg <- head(fr[fr$NES < 0, ][order(fr[fr$NES < 0, ]$padj), ], 10)
  tp <- dplyr::bind_rows(top_pos, top_neg)
  if (nrow(tp) == 0) next
  tp$pathway_short <- gsub("_", " ", gsub("^HALLMARK_", "", tp$pathway))
  tp$highlight <- grepl("DNA|REPAIR|REPLIC|DAMAGE|STRESS", tp$pathway_short, ignore.case = TRUE)
  p <- ggplot(tp, aes(reorder(pathway_short, NES), NES, fill = padj < 0.05)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey70"),
      labels = c("TRUE" = "padj < 0.05", "FALSE" = "NS"), name = "Significance") +
    coord_flip() + labs(title = paste("GSEA Hallmark —", nm), x = "",
      y = "Normalised Enrichment Score (NES)") +
    theme_cowplot() +
    theme(axis.text.y = element_text(face = ifelse(tp$highlight[order(tp$NES)], "bold", "plain")))
  save_plot(p, paste0("gsea_hallmark_", nm), "03_enrichment", fig_root, width = 12, height = 8)
}

# --- 4c. Targeted C5 GO:BP pathways for dKO ---------------------------------
c5 <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
c5_list <- lapply(split(c5$entrez_gene, c5$gs_name), function(x) as.character(unique(x)))
target_patterns <- c(
  DNA_REPAIR             = "DNA_REPAIR",
  DNA_DAMAGE             = "RESPONSE_TO_DNA_DAMAGE|DNA_DAMAGE_RESPONSE",
  REPLICATION_STRESS     = "DNA_REPLICATION|REPLICATION_FORK|STALLED_REPLICATION",
  ATR_ATM_SIGNALLING     = "ATR|ATM|SIGNAL_TRANSDUCTION.*DAMAGE",
  G_QUADRUPLEX           = "QUADRUPLEX|G_QUADRUPLEX",
  R_LOOP_RNA_DNA_HYBRID  = "R_LOOP|RNA_DNA_HYBRID")
ranks_dko <- build_ranks(res_list[["dKO_vs_WT"]])
target_results <- lapply(names(target_patterns), function(pn) {
  matching <- grep(target_patterns[[pn]], names(c5_list), value = TRUE, ignore.case = TRUE)
  if (length(matching) == 0) return(NULL)
  fr <- fgsea::fgsea(pathways = c5_list[matching], stats = ranks_dko,
                     minSize = rc$gsea_targeted_minsize, maxSize = rc$gsea_targeted_maxsize)
  fr$category <- pn; fr
})
target_df <- dplyr::bind_rows(target_results)
if (!is.null(target_df) && nrow(target_df) > 0) {
  target_df$pathway_short <- gsub("_", " ", gsub("^GOBP_", "", target_df$pathway))
  top <- head(target_df[order(target_df$padj), ], 20)
  p <- ggplot(top, aes(reorder(pathway_short, NES), NES, fill = category)) +
    geom_col() + coord_flip() +
    labs(title = "Targeted pathway enrichment — dKO vs WT", x = "", y = "NES", fill = "Category") +
    theme_cowplot() + theme(axis.text.y = element_text(size = 8))
  save_plot(p, "gsea_specific_pathways_dKO", "03_enrichment", fig_root, width = 10, height = 8)
}
message("Enrichment done.")
