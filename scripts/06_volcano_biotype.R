#!/usr/bin/env Rscript
# ============================================================================
# 22_volcano_plots.R   (env: system R 4.5.1)   --- Phase 9 ---
#
# RNA-seq volcano plots, one pair per KO-vs-WT contrast (DHX36KO, FANCJKO, dKO),
# with up/down DEGs DISTINCTLY coloured (up = #E41A1C, down = #377EB8, NS grey).
# Two flavours per contrast:
#   (A) <ko>_top50  -- all points coloured by direction; label the N most up +
#                      N most down significant genes by |log2FC| (config
#                      lncrna.volcano_label_n, default 50).
#   (B) <ko>_lncRNA -- non-lncRNA points greyed out, lncRNA points coloured by
#                      direction; label significant lncRNAs (top N by |log2FC|
#                      per direction). Focuses the figure on lncRNA regulation.
#
# Significance = padj < deg_padj AND |log2FC| > deg_lfc (project standard).
# Biotype (lncRNA vs other) joined from rule 21 by unversioned ENSEMBL gene_id.
#
# Inputs:   results/tables/deseq2_KO_vs_WT.csv (rule 13a),
#           results/tables/gene_biotype.csv    (rule 21)
# Deps:     dplyr, readr, ggplot2, ggrepel
# Outputs:  results/figures/16_volcano/volcano_<ko>_{top50,lncRNA}.{pdf,png}
#           results/tables/volcano_label_genes.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/16_volcano")
fig_root <- "results/figures"

deg_padj <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc  <- cfg$metaprofile$deg_lfc  %||% 0.5
label_n  <- cfg$lncrna$volcano_label_n %||% 50
lnc_biotypes <- as.character(unlist(cfg$lncrna$biotypes %||% cfg$lncrna$biotype %||% "lncRNA"))
dir_cols <- c(up = "#E41A1C", down = "#377EB8", NS = "grey80")
strip_ver <- function(x) sub("\\..*$", "", as.character(x))

# --- Load DESeq2 + biotype; classify direction ------------------------------
de <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
bt <- readr::read_csv("results/tables/gene_biotype.csv", show_col_types = FALSE)
de$gene_id <- strip_ver(de$gene_id)
bt$gene_id <- strip_ver(bt$gene_id)
de <- de %>% dplyr::left_join(bt %>% dplyr::select(gene_id, gene_type), by = "gene_id")
de$is_lncRNA <- !is.na(de$gene_type) & de$gene_type %in% lnc_biotypes
matched <- mean(de$gene_id %in% bt$gene_id)
message(sprintf("Biotype match: %.1f%% of DESeq2 gene_ids; %d lncRNA rows",
                100 * matched, sum(de$is_lncRNA, na.rm = TRUE)))

de <- de %>% dplyr::filter(!is.na(padj), !is.na(log2FC))
de$direction <- dplyr::case_when(
  de$padj < deg_padj & de$log2FC >  deg_lfc ~ "up",
  de$padj < deg_padj & de$log2FC < -deg_lfc ~ "down",
  TRUE ~ "NS")
# -log10(padj) with a finite cap so padj==0 genes still plot.
finite_max <- max(-log10(de$padj[de$padj > 0]), na.rm = TRUE)
de$negLog10 <- ifelse(de$padj <= 0, finite_max * 1.05, -log10(de$padj))

kos <- intersect(c("DHX36KO", "FANCJKO", "dKO"), unique(as.character(de$ko)))
hline <- -log10(deg_padj)

base_volcano <- function(df, title, subtitle, point_aes_col) {
  ggplot(df, aes(x = log2FC, y = negLog10)) +
    geom_vline(xintercept = c(-deg_lfc, deg_lfc), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = hline, linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = .data[[point_aes_col]]), size = 1.1, alpha = 0.7) +
    scale_colour_manual(values = dir_cols,
                        labels = c(up = "Up", down = "Down", NS = "n.s."),
                        breaks = c("up", "down", "NS"), name = NULL,
                        na.value = "grey85") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(p)),
         title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "top")
}

label_rows <- list()
for (k in kos) {
  dk <- de %>% dplyr::filter(ko == k)
  sig <- dk %>% dplyr::filter(direction %in% c("up", "down"))

  # ---- (A) normal volcano: top N up + N down by |log2FC| -------------------
  topA <- sig %>% dplyr::group_by(direction) %>%
    dplyr::slice_max(abs(log2FC), n = label_n, with_ties = FALSE) %>% dplyr::ungroup()
  pA <- base_volcano(
    dk, sprintf("%s vs WT", k),
    sprintf("padj<%g & |log2FC|>%g; %d up / %d down (top %d/dir labelled)",
            deg_padj, deg_lfc, sum(dk$direction == "up"), sum(dk$direction == "down"), label_n),
    "direction") +
    ggrepel::geom_text_repel(data = topA, aes(label = gene_name), size = 2.6,
                             max.overlaps = 25, segment.size = 0.2, colour = "grey15")
  save_plot(pA, sprintf("volcano_%s_top50", k), "16_volcano", fig_root, width = 7.5, height = 7)

  # ---- (B) lncRNA-focused volcano ------------------------------------------
  dk_lnc <- dk %>% dplyr::mutate(
    lnc_colour = ifelse(is_lncRNA & direction != "NS", direction, "NS"))
  sig_lnc <- dk %>% dplyr::filter(is_lncRNA, direction %in% c("up", "down"))
  topB <- sig_lnc %>% dplyr::group_by(direction) %>%
    dplyr::slice_max(abs(log2FC), n = label_n, with_ties = FALSE) %>% dplyr::ungroup()
  pB <- base_volcano(
    dk_lnc, sprintf("%s vs WT  -  lncRNAs highlighted", k),
    sprintf("non-lncRNA genes greyed; %d lncRNA up / %d lncRNA down (padj<%g & |log2FC|>%g)",
            sum(sig_lnc$direction == "up"), sum(sig_lnc$direction == "down"), deg_padj, deg_lfc),
    "lnc_colour") +
    ggrepel::geom_text_repel(data = topB, aes(label = gene_name), size = 2.6,
                             max.overlaps = 30, segment.size = 0.2, colour = "grey15")
  save_plot(pB, sprintf("volcano_%s_lncRNA", k), "16_volcano", fig_root, width = 7.5, height = 7)

  if (nrow(topA) > 0) label_rows[[length(label_rows) + 1]] <-
    topA %>% dplyr::transmute(ko = k, gene_name, gene_id, direction, log2FC, padj,
                              is_lncRNA, plot = "top50")
  if (nrow(topB) > 0) label_rows[[length(label_rows) + 1]] <-
    topB %>% dplyr::transmute(ko = k, gene_name, gene_id, direction, log2FC, padj,
                              is_lncRNA, plot = "lncRNA")
  message(sprintf("%s: %d up / %d down DEGs; %d lncRNA DEGs labelled",
                  k, sum(dk$direction == "up"), sum(dk$direction == "down"), nrow(topB)))
}

labels_out <- dplyr::bind_rows(label_rows)
readr::write_csv(labels_out, "results/tables/volcano_label_genes.csv")
message("Done. Volcano plots for ", length(kos), " contrasts (top50 + lncRNA each).")
