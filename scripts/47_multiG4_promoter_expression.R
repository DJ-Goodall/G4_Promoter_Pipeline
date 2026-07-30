#!/usr/bin/env Rscript
# ============================================================================
# 32_multiG4_promoter_expression.R  (env: r_g4)  --- Phase 13, Analysis 4 ---
#
# How many promoters carry MULTIPLE G4 motifs within TSS +/- 2 kb, and how does the
# G4 count relate to gene activity? Two RNA-seq proxies:
#   (1) BASELINE  : DESeq2 baseMean (steady-state expression) vs G4-count bucket.
#   (2) KO RESPONSE: |log2FC| and DEG fraction (KO vs WT) vs G4-count bucket, per KO.
# Plus the topology composition of multi-G4 promoters.
#
# Inputs:   cache_v2 regions_promoter.rds (Entrez names), cache/motifs_all.rds,
#           results/tables/deseq2_KO_vs_WT.csv
# Outputs:  results/tables/multiG4_promoter_counts.csv, multiG4_expression_stats.csv
#           results/figures/21_multiG4_expression/{promoter_g4_count,
#             g4count_vs_baseline_expression,g4count_vs_ko_response,
#             g4count_deg_fraction,g4count_topology_composition}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges); library(S4Vectors)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/21_multiG4_expression")
fig_root <- "results/figures"

kos       <- setdiff(cfg$genotypes, cfg$ref_genotype)
hw        <- cfg$motif_analysis$promoter_half_width %||% 2000
deg_padj  <- cfg$metaprofile$deg_padj %||% 0.05
deg_lfc   <- cfg$metaprofile$deg_lfc %||% 0.5
definite  <- c("parallel", "antiparallel", "hybrid")
bucket_of <- function(n) cut(n, breaks = c(-1, 0, 1, 2, 3, Inf),
                             labels = c("0", "1", "2", "3", "4+"))

# --- Count motifs per promoter (TSS +/- hw) ---------------------------------
prom <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$promoter))
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), cfg$std_chroms), pruning.mode = "coarse")
prom_entrez <- names(prom); if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)
tss      <- GenomicRanges::resize(prom, width = 1, fix = "center")
prom_win <- GenomicRanges::resize(tss, width = 2 * hw, fix = "center")

motifs_gr <- readRDS("cache/motifs_all.rds")
motif_cen <- as.integer(S4Vectors::mcols(motifs_gr)$motif_center)
motif_topo <- as.character(S4Vectors::mcols(motifs_gr)$topology)
motif_pts <- GenomicRanges::GRanges(GenomicRanges::seqnames(motifs_gr),
                                    IRanges::IRanges(motif_cen, width = 1))

ov <- GenomicRanges::findOverlaps(motif_pts, prom_win, ignore.strand = TRUE)
n_g4 <- tabulate(S4Vectors::subjectHits(ov), nbins = length(prom_win))

prom_df <- data.frame(entrez = prom_entrez, n_g4 = as.integer(n_g4), stringsAsFactors = FALSE)
prom_df$bucket <- bucket_of(prom_df$n_g4)
prom_df$symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_df$entrez,
                                        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# --- Join expression (Entrez -> SYMBOL -> DESeq2 gene_name) ------------------
de <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
baseline <- de %>% group_by(gene_name) %>%
  summarise(baseMean = mean(baseMean, na.rm = TRUE), .groups = "drop")

prom_expr <- prom_df %>% filter(!is.na(symbol)) %>%
  left_join(baseline, by = c("symbol" = "gene_name"))
readr::write_csv(prom_expr, "results/tables/multiG4_promoter_counts.csv")

# --- Promoter G4-count distribution ------------------------------------------
bucket_n <- prom_df %>% count(bucket, name = "n_promoters")
p_cnt <- ggplot(bucket_n, aes(bucket, n_promoters)) +
  geom_col(fill = "#1B9E77", colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n_promoters), vjust = -0.3, size = 3) +
  labs(x = "G4 motifs in promoter (TSS +/- 2 kb)", y = "Promoters",
       title = "Promoters by number of G4 motifs",
       subtitle = sprintf("%.1f%% of promoters carry >=2 G4 motifs",
                          100 * mean(prom_df$n_g4 >= 2))) +
  theme_pub()
save_plot(p_cnt, "promoter_g4_count", "21_multiG4_expression", fig_root, width = 7, height = 5)

# --- (1) Baseline expression vs G4 count ------------------------------------
be <- prom_expr %>% filter(is.finite(baseMean))
sp <- suppressWarnings(cor.test(be$n_g4, log10(be$baseMean + 1), method = "spearman"))
p_base <- ggplot(be, aes(bucket, log10(baseMean + 1))) +
  geom_boxplot(outlier.shape = NA, fill = "#4DAF4A", alpha = 0.8) +
  labs(x = "G4 motifs in promoter", y = "log10(baseMean + 1)",
       title = "Baseline expression rises with promoter G4 count",
       subtitle = sprintf("Spearman rho = %.3f (p = %.2g); DESeq2 baseMean as activity proxy",
                          unname(sp$estimate), sp$p.value)) +
  theme_pub()
save_plot(p_base, "g4count_vs_baseline_expression", "21_multiG4_expression", fig_root,
          width = 7, height = 5)

# --- (2) KO response vs G4 count --------------------------------------------
resp <- prom_expr %>% dplyr::select(symbol, n_g4, bucket) %>%
  inner_join(de, by = c("symbol" = "gene_name")) %>%
  mutate(abs_lfc = abs(log2FC),
         is_deg = !is.na(padj) & padj < deg_padj & abs(log2FC) > deg_lfc,
         ko = factor(ko, levels = kos))

p_resp <- ggplot(resp %>% filter(is.finite(abs_lfc)), aes(bucket, abs_lfc)) +
  geom_boxplot(outlier.shape = NA, fill = "#E41A1C", alpha = 0.75) +
  facet_wrap(~ ko, nrow = 1) +
  coord_cartesian(ylim = c(0, quantile(resp$abs_lfc, 0.97, na.rm = TRUE))) +
  labs(x = "G4 motifs in promoter", y = "|log2 fold-change| (KO vs WT)",
       title = "KO transcriptional response vs promoter G4 count",
       subtitle = "Are multi-G4 promoters more dysregulated on helicase loss?") +
  theme_pub()
save_plot(p_resp, "g4count_vs_ko_response", "21_multiG4_expression", fig_root, width = 11, height = 4.5)

deg_frac <- resp %>% group_by(ko, bucket) %>%
  summarise(n = dplyr::n(), deg_frac = mean(is_deg), .groups = "drop")
p_frac <- ggplot(deg_frac, aes(bucket, deg_frac, fill = ko)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = condition_colours, name = "KO") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "G4 motifs in promoter", y = "DEG fraction",
       title = "Fraction of genes called DEG vs promoter G4 count",
       subtitle = sprintf("DEG = padj<%g & |log2FC|>%g", deg_padj, deg_lfc)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_frac, "g4count_deg_fraction", "21_multiG4_expression", fig_root, width = 8, height = 5)

# --- (3) Topology composition by G4-count bucket ----------------------------
pair_df <- data.frame(prom_i = S4Vectors::subjectHits(ov),
                      topology = motif_topo[S4Vectors::queryHits(ov)],
                      stringsAsFactors = FALSE)
pair_df$bucket <- prom_df$bucket[pair_df$prom_i]
comp <- pair_df %>% filter(topology %in% definite, bucket != "0") %>%
  count(bucket, topology, name = "n") %>%
  group_by(bucket) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
comp$topology <- factor(comp$topology, levels = definite)
p_comp <- ggplot(comp, aes(bucket, pct, fill = topology)) +
  geom_col(colour = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = topology_palette,
                    labels = setNames(topology_labels[match(definite, topology_levels)], definite),
                    name = "Topology") +
  labs(x = "G4 motifs in promoter", y = "Motif topology share (%)",
       title = "Topology mix of motifs in multi-G4 promoters",
       subtitle = "Does topology composition shift as promoters accumulate more G4s?") +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_comp, "g4count_topology_composition", "21_multiG4_expression", fig_root,
          width = 7, height = 5)

# --- Stats table ------------------------------------------------------------
stat_rows <- list(data.frame(
  metric = "spearman_nG4_vs_log_baseMean",
  ko = NA_character_, statistic = unname(sp$estimate), p_value = sp$p.value,
  stringsAsFactors = FALSE))
for (k in kos) {
  d <- resp %>% filter(ko == k, is.finite(abs_lfc))
  kw <- tryCatch(kruskal.test(abs_lfc ~ bucket, data = d), error = function(e) NULL)
  stat_rows[[length(stat_rows) + 1]] <- data.frame(
    metric = "kruskal_absLFC_by_bucket", ko = k,
    statistic = if (is.null(kw)) NA_real_ else unname(kw$statistic),
    p_value = if (is.null(kw)) NA_real_ else kw$p.value, stringsAsFactors = FALSE)
}
stats <- bind_rows(stat_rows)
readr::write_csv(stats, "results/tables/multiG4_expression_stats.csv")

message("Done. Analysis 4 (multi-G4 promoter expression). Promoters per bucket:")
print(bucket_n)
print(stats)
