#!/usr/bin/env Rscript
# ============================================================================
# 20_unclassified_peak_motifs.R   (env: system R 4.5.1)   --- Phase 8 ---
#
# Is there a TF-motif signature in the no-canonical-PQS ("unclassified") G4
# peaks, or are they sequence-random? Known-motif (JASPAR2020 CORE) enrichment
# via motifmatchr, tested as:
#   (a) no-PQS vs classified G4 peaks  [headline contrast the user asked for]
#   (b) no-PQS vs GC-matched background [GC control -- classified peaks are
#       GC-rich, so (a) alone is GC-confounded; a motif that is enriched in (a)
#       only and NOT in (b) is a low-GC artifact, whereas one enriched in BOTH
#       is a genuinely specific motif].
# Each motif is annotated with its own GC so GC-driven hits are visible.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           cache/peak_seqs_201bp.fa, cache/background_seqs.fa
# Deps (system R): motifmatchr, TFBSTools, JASPAR2020 (+ Matrix, ggrepel)
# Outputs:  results/tables/unclassified_motif_enrichment.csv
#           results/figures/15_unclassified/unclassified_peak_motifs.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(Biostrings); library(GenomicRanges); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2); library(cowplot)
  library(ggrepel); library(Matrix)
  library(TFBSTools); library(motifmatchr); library(JASPAR2020)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/15_unclassified")
fig_root <- "results/figures"

# --- Peaks -> class; sequence sets -----------------------------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
union$class2 <- ifelse(union$topology == "no_canonical_PQS", "unclassified", "classified")

peak_seqs <- Biostrings::readDNAStringSet("cache/peak_seqs_201bp.fa")
bg_seqs   <- Biostrings::readDNAStringSet("cache/background_seqs.fa")
cl <- union$class2[match(names(peak_seqs), names(union))]

set.seed(7)
cap <- function(x, n) if (length(x) > n) x[sample(length(x), n)] else x
sub <- cfg$unclassified$subsample %||% 16980
unc <- cap(peak_seqs[which(cl == "unclassified")], sub)
cla <- cap(peak_seqs[which(cl == "classified")], length(unc))   # equal n to no-PQS
bg  <- cap(bg_seqs, length(unc))
allseqs <- c(unc, cla, bg)
grp <- rep(c("unclassified", "classified", "background"),
           c(length(unc), length(cla), length(bg)))
message("Sequence sets: ", length(unc), " no-PQS / ", length(cla),
        " classified / ", length(bg), " background")

# --- JASPAR2020 CORE motifs + scan -----------------------------------------
pfms <- TFBSTools::getMatrixSet(
  JASPAR2020::JASPAR2020,
  list(collection = "CORE", tax_group = cfg$unclassified$jaspar_tax %||% "vertebrates"))
pcut <- cfg$unclassified$motif_pcut %||% 5e-5
message("Scanning ", length(allseqs), " sequences x ", length(pfms),
        " JASPAR motifs (p.cutoff ", pcut, ") ...")
mm   <- motifmatchr::matchMotifs(pfms, allseqs, out = "matches", p.cutoff = pcut, bg = "even")
hits <- motifmatchr::motifMatches(mm)        # sparse lgCMatrix: rows = seqs, cols = motifs

idx <- list(unclassified = which(grp == "unclassified"),
            classified   = which(grp == "classified"),
            background   = which(grp == "background"))
cnt <- sapply(idx, function(i) Matrix::colSums(hits[i, , drop = FALSE]))  # n_motif x 3
nn  <- vapply(idx, length, integer(1))

ids <- vapply(pfms, TFBSTools::ID, character(1))
tfs <- vapply(pfms, TFBSTools::name, character(1))
pfm_gc <- vapply(pfms, function(m) {
  M <- TFBSTools::Matrix(m); cs <- colSums(M); mean((M["C", ] + M["G", ]) / cs)
}, numeric(1))

# Per-motif 2x2 Fisher: no-PQS hit-fraction vs a reference set's hit-fraction.
# as.numeric() strips names that mapply propagates from the named count vectors
# (otherwise c(log2OR=...) builds compound names like "log2OR.MA0001").
fisher_or <- function(a, n1, c, n2) {
  a <- as.numeric(a); c <- as.numeric(c); n1 <- as.numeric(n1); n2 <- as.numeric(n2)
  b <- n1 - a; d <- n2 - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  c(log2OR = log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))),
    p = as.numeric(ft$p.value))
}
res_cla <- t(mapply(function(a, c) fisher_or(a, nn["unclassified"], c, nn["classified"]),
                    cnt[, "unclassified"], cnt[, "classified"], USE.NAMES = FALSE))
res_bg  <- t(mapply(function(a, c) fisher_or(a, nn["unclassified"], c, nn["background"]),
                    cnt[, "unclassified"], cnt[, "background"], USE.NAMES = FALSE))

out <- data.frame(
  motif_id = ids, tf = tfs, motif_gc = round(pfm_gc, 3),
  frac_unclassified = cnt[, "unclassified"] / nn["unclassified"],
  frac_classified   = cnt[, "classified"]   / nn["classified"],
  frac_background   = cnt[, "background"]    / nn["background"],
  log2OR_vs_classified = res_cla[, "log2OR"], p_vs_classified = res_cla[, "p"],
  log2OR_vs_background = res_bg[, "log2OR"],  p_vs_background  = res_bg[, "p"],
  stringsAsFactors = FALSE)
out$q_vs_classified <- p.adjust(out$p_vs_classified, "BH")
out$q_vs_background <- p.adjust(out$p_vs_background, "BH")
out <- out[order(-out$log2OR_vs_background), ]
readr::write_csv(out, "results/tables/unclassified_motif_enrichment.csv")
message("Top motifs enriched in no-PQS vs GC-matched background:")
print(head(out[, c("tf", "motif_gc", "frac_unclassified", "frac_background",
                   "log2OR_vs_background", "q_vs_background")], 12))

# --- Figure: specificity scatter + top GC-controlled motifs -----------------
out$specific <- out$q_vs_background < 0.05 & out$log2OR_vs_background > 0
p_sc <- ggplot(out, aes(log2OR_vs_classified, log2OR_vs_background, colour = specific)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(alpha = 0.6, size = 1.3) +
  ggrepel::geom_text_repel(
    data = head(out, 12), aes(label = tf), size = 2.8, max.overlaps = 20, colour = "black") +
  scale_colour_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                      labels = c("TRUE" = "q<0.05 vs GC-bg", "FALSE" = "ns"), name = NULL) +
  labs(x = "log2 OR (no-PQS vs classified)", y = "log2 OR (no-PQS vs GC-matched bg)",
       title = "TF-motif enrichment in no-canonical-PQS peaks",
       subtitle = "Upper-right = specific (enriched vs both); lower-right = GC artifact (vs classified only)") +
  theme_pub() + theme(legend.position = "top")

top <- head(out, 15)
top$lab <- factor(sprintf("%s (%s)", top$tf, top$motif_id),
                  levels = rev(sprintf("%s (%s)", top$tf, top$motif_id)))
p_bar <- ggplot(top, aes(log2OR_vs_background, lab, fill = q_vs_background < 0.05)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey60"),
                    labels = c("TRUE" = "q<0.05", "FALSE" = "ns"), name = NULL) +
  labs(x = "log2 OR vs GC-matched background", y = NULL,
       title = "Top 15 motifs (GC-controlled)") +
  theme_pub() + theme(legend.position = "top")
save_plot(cowplot::plot_grid(p_sc, p_bar, ncol = 1, rel_heights = c(1.3, 1)),
          "unclassified_peak_motifs", "15_unclassified", fig_root, width = 9, height = 11)

message("Done. Known-motif (JASPAR) enrichment for unclassified peaks.")
