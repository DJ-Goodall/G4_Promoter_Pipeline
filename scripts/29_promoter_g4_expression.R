#!/usr/bin/env Rscript
# ============================================================================
# 13b_promoter_g4_expression.R   (env: r_g4)   --- Q1d ---
#
# Link promoter-G4 topology to gene expression and KO-induced change. Headline:
# are genes with a PARALLEL promoter G4 preferentially dysregulated in DHX36KO?
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           cache_v2 regions_promoter.rds (names = Entrez gene_id),
#           results/tables/deseq2_KO_vs_WT.csv
# Outputs:  results/tables/promoter_g4_expression.csv,
#           results/tables/promoter_g4_expression_tests.csv,
#           results/figures/14_expression/promoter_g4_expression.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/14_expression")
fig_root <- "results/figures"

kos      <- setdiff(cfg$genotypes, cfg$ref_genotype)
definite <- c("parallel", "antiparallel", "hybrid")

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
has_maxz <- "max_z" %in% colnames(S4Vectors::mcols(union))

prom <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$promoter))
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom_entrez <- names(prom)
if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)

# --- Assign each promoter the topology of its strongest overlapping G4 ------
ov <- GenomicRanges::findOverlaps(prom, union, ignore.strand = TRUE)
ov_df <- data.frame(prom_i = S4Vectors::queryHits(ov),
                    peak_i = S4Vectors::subjectHits(ov),
                    topology = union$topology[S4Vectors::subjectHits(ov)],
                    strength = if (has_maxz) S4Vectors::mcols(union)$max_z[S4Vectors::subjectHits(ov)]
                               else GenomicRanges::width(union)[S4Vectors::subjectHits(ov)],
                    stringsAsFactors = FALSE) %>%
  filter(!is.na(topology)) %>%
  group_by(prom_i) %>% slice_max(strength, n = 1, with_ties = FALSE) %>% ungroup()

prom_df <- data.frame(entrez = prom_entrez, g4_topology = "no_G4", stringsAsFactors = FALSE)
prom_df$g4_topology[ov_df$prom_i] <- ov_df$topology
prom_df$g4_topology[prom_df$g4_topology == "no_canonical_PQS"] <- "no_PQS_G4"

# --- Map Entrez -> SYMBOL and join DESeq2 ----------------------------------
prom_df$symbol <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = prom_df$entrez,
                                        column = "SYMBOL", keytype = "ENTREZID",
                                        multiVals = "first")
prom_df <- prom_df %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE)

de <- readr::read_csv("results/tables/deseq2_KO_vs_WT.csv", show_col_types = FALSE)
expr <- de %>% inner_join(prom_df, by = c("gene_name" = "symbol")) %>%
  mutate(abs_log2FC = abs(log2FC),
         g4_topology = factor(g4_topology,
                              levels = c("parallel", "antiparallel", "hybrid",
                                         "no_PQS_G4", "no_G4")),
         ko = factor(ko, levels = kos))
readr::write_csv(expr, "results/tables/promoter_g4_expression.csv")

# --- Tests: dysregulation by promoter-G4 topology, per KO -------------------
test_rows <- list()
for (k in kos) {
  d <- expr %>% filter(ko == k, is.finite(abs_log2FC))
  kw <- tryCatch(kruskal.test(abs_log2FC ~ g4_topology, data = d), error = function(e) NULL)
  # focused: parallel-G4 genes vs no-G4 genes
  wp <- tryCatch(wilcox.test(d$abs_log2FC[d$g4_topology == "parallel"],
                             d$abs_log2FC[d$g4_topology == "no_G4"])$p.value,
                 error = function(e) NA_real_)
  test_rows[[length(test_rows) + 1]] <- data.frame(
    ko = k,
    kruskal_p = if (is.null(kw)) NA_real_ else kw$p.value,
    wilcox_parallel_vs_noG4_p = wp,
    median_absLFC_parallel = median(d$abs_log2FC[d$g4_topology == "parallel"], na.rm = TRUE),
    median_absLFC_noG4 = median(d$abs_log2FC[d$g4_topology == "no_G4"], na.rm = TRUE),
    stringsAsFactors = FALSE)
}
tests <- dplyr::bind_rows(test_rows)
readr::write_csv(tests, "results/tables/promoter_g4_expression_tests.csv")

p_expr <- ggplot(expr %>% filter(is.finite(abs_log2FC)),
                 aes(g4_topology, abs_log2FC, fill = g4_topology)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  facet_wrap(~ ko) +
  coord_cartesian(ylim = c(0, quantile(expr$abs_log2FC, 0.97, na.rm = TRUE))) +
  scale_fill_manual(values = c(parallel = "#1B9E77", antiparallel = "#D95F02",
                               hybrid = "#7570B3", no_PQS_G4 = "grey75", no_G4 = "grey55"),
                    guide = "none") +
  labs(x = NULL, y = "|log2 fold-change| (KO vs WT)",
       title = "Gene dysregulation by promoter-G4 topology",
       subtitle = "Headline: are parallel-promoter-G4 genes more dysregulated, esp. in DHX36KO?") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_expr, "promoter_g4_expression", "14_expression", fig_root, width = 11, height = 5)

message("Done. Promoter-G4 expression tests:")
print(tests)
