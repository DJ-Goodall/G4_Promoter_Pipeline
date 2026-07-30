#!/usr/bin/env Rscript
# ============================================================================
# 48_strand_asymmetry.R   (system R)   --- Stage L ---
#
# G4 strand-asymmetry analysis. Ported from g4_gloop_strand_analysis.Rmd
# (chunks s1-s8). For each G4 peak in promoter/5'UTR regions: extract flanking
# sequence at each window size, score G4Hunter on both strands (custom impl),
# assign the peak to the higher-scoring strand, then test genomic-strand bias
# (binomial), gene-relative bias, G4-score distributions (Wilcoxon), detection
# rate (chi-square), and build strand + template-vs-coding metaprofiles.
# gene_name is populated (fixes the blank-label bug used by the Gviz rule 49).
#
# Inputs:   cache/peaks_G4_BG4_{WT,DHX36KO,FANCJKO,dKO}.rds,
#           cache/regions_{promoter,5UTR}.rds, config bigwig_dir
# Outputs:  cache/sequences_scored.rds
#           results/tables/{strand_bias_summary,strand_bias_relative,
#             g4_detection_by_strand,g4_score_stats}.tsv
#           results/figures/22_strand/{strand_bias_comprehensive,
#             metaprofile_by_strand,mean_signal_by_strand,
#             metaprofile_template_vs_coding}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(BSgenome.Mmusculus.UCSC.mm10); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db); library(AnnotationDbi); library(rtracklayer)
  library(IRanges); library(ggplot2); library(cowplot); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache", "results/figures", "results/tables")
fig_root <- "results/figures"
sc <- cfg$strand
window_sizes   <- as.integer(unlist(sc$window_sizes))
g4_threshold   <- sc$g4_threshold
genotype_order <- c("WT", "DHX36KO", "FANCJKO", "dKO")
region_types   <- c("promoter", "5'UTR")
n_bins <- sc$n_bins; half_width <- sc$half_width

genome <- BSgenome.Mmusculus.UCSC.mm10
txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene

# --- G4Hunter (custom, faithful to the RMD) ---------------------------------
g4_hunter_score <- function(seq_str) {
  if (is.na(seq_str) || nchar(seq_str) == 0) return(NA_real_)
  chars <- strsplit(toupper(seq_str), "")[[1]]
  r <- rle(chars)
  val <- ifelse(r$values == "G",  pmin(r$lengths, 4L),
         ifelse(r$values == "C", -pmin(r$lengths, 4L), 0L))
  s <- rep(val, r$lengths); m <- mean(s, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

# --- Load peaks + assign region_type ----------------------------------------
g4_peaks_list <- lapply(genotype_order, function(g) {
  f <- file.path("cache", sprintf("peaks_G4_BG4_%s.rds", g))
  pk <- readRDS(f); pk$genotype <- g; pk$assay <- "G4_BG4"; pk })
g4_peaks <- Reduce(c, g4_peaks_list)

reg <- list(promoter = readRDS("cache/regions_promoter.rds"),
            `5'UTR`   = readRDS("cache/regions_5UTR.rds"))
g4_peaks$region_type <- NA_character_
for (rt in region_types) {
  ov <- GenomicRanges::findOverlaps(g4_peaks, reg[[rt]])
  idx <- unique(S4Vectors::queryHits(ov))
  idx <- idx[is.na(g4_peaks$region_type[idx])]
  g4_peaks$region_type[idx] <- rt
}
g4_peaks <- g4_peaks[!is.na(g4_peaks$region_type)]

# subsample 2000/genotype
set.seed(sc$seed)
g4_peaks <- do.call(c, lapply(genotype_order, function(g) {
  pg <- g4_peaks[g4_peaks$genotype == g]
  if (length(pg) > sc$subsample) pg[sample(length(pg), sc$subsample)] else pg }))

# gene_name per peak (bug fix for Gviz): overlapping gene -> SYMBOL
genes_gr <- GenomicFeatures::genes(txdb)
ov <- GenomicRanges::findOverlaps(g4_peaks, genes_gr, select = "first")
entrez <- names(genes_gr)[ov]
# NA keys (peaks not overlapping any gene) make mapIds() return a *list*, which
# then explodes into per-element columns inside data.frame() below (-> "differing
# number of rows: 1, 0"). Map only the non-NA keys and keep a plain char vector.
sym <- rep(NA_character_, length(g4_peaks))
valid <- !is.na(entrez)
if (any(valid)) {
  m <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = entrez[valid],
        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
  sym[valid] <- unname(as.character(m))
}
g4_peaks$gene_name <- sym

# --- Extract sequences (per region_type x window) ---------------------------
scored_cache <- "cache/sequences_scored.rds"
sequences_df <- cache_or_build(scored_cache, {
  extract_seqs <- function(gr, window_bp) {
    ex <- GenomicRanges::trim(GenomicRanges::resize(gr, width = 2 * window_bp + 1, fix = "center"))
    seqs <- BSgenome::getSeq(genome, ex)
    data.frame(seqname = as.character(GenomicRanges::seqnames(gr)),
      peak_start = GenomicRanges::start(gr), peak_end = GenomicRanges::end(gr),
      strand = as.character(GenomicRanges::strand(gr)), genotype = gr$genotype,
      gene_name = gr$gene_name, region_type = gr$region_type, window_size = window_bp,
      sequence = as.character(seqs), stringsAsFactors = FALSE)
  }
  parts <- list()
  for (rt in region_types) for (w in window_sizes) {
    pk <- g4_peaks[g4_peaks$region_type == rt]
    if (length(pk) == 0) next
    parts[[paste(rt, w)]] <- extract_seqs(pk, w)
  }
  sdf <- do.call(rbind, parts); rownames(sdf) <- NULL

  # both-strand G4Hunter scoring
  fwd <- vapply(sdf$sequence, g4_hunter_score, numeric(1), USE.NAMES = FALSE)
  rev <- vapply(sdf$sequence, function(s) g4_hunter_score(
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))),
    numeric(1), USE.NAMES = FALSE)
  det <- ifelse(is.na(fwd) & is.na(rev), "neither",
         ifelse(is.na(rev) | (!is.na(fwd) & fwd >= rev), "+", "-"))
  sdf$g4_score_fwd <- fwd; sdf$g4_score_rev <- rev
  sdf$g4_detected_strand <- det
  sdf$has_g4_motif <- pmax(abs(fwd), abs(rev), na.rm = TRUE) >= g4_threshold
  sdf$g4_score <- ifelse(det == "+", fwd, rev)
  sdf
})

# --- Gene-relative strand + template/coding ---------------------------------
gr_peaks <- GenomicRanges::GRanges(sequences_df$seqname,
  IRanges::IRanges(sequences_df$peak_start, sequences_df$peak_end),
  strand = ifelse(sequences_df$g4_detected_strand %in% c("+", "-"),
                  sequences_df$g4_detected_strand, "*"))
gs_ov <- GenomicRanges::findOverlaps(gr_peaks, genes_gr, select = "first")
sequences_df$gene_strand <- as.character(GenomicRanges::strand(genes_gr))[gs_ov]
sequences_df$strand_class <- with(sequences_df, ifelse(
  is.na(gene_strand) | g4_detected_strand == "neither", NA_character_,
  ifelse(g4_detected_strand == gene_strand, "with_strand", "against_strand")))
sequences_df$g4_strand_class <- with(sequences_df, ifelse(
  is.na(gene_strand) | g4_detected_strand == "neither", "unmapped",
  ifelse(g4_detected_strand == gene_strand, "coding", "template")))

# --- Statistics tables ------------------------------------------------------
grid <- expand.grid(genotype = genotype_order, region_type = region_types,
                    window_size = window_sizes, stringsAsFactors = FALSE)
strand_bias_summary <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i], ]
  cp <- sum(s$g4_detected_strand == "+"); cm <- sum(s$g4_detected_strand == "-")
  nt <- cp + cm
  data.frame(grid[i, ], count_plus = cp, count_minus = cm, n_total = nt,
    pct_plus = if (nt > 0) 100 * cp / nt else NA,
    p_value_binom = if (nt > 0) binom.test(cp, nt, 0.5)$p.value else NA,
    fold_plus = if (nt > 0) (cp / nt) / 0.5 else NA)
}))
strand_bias_summary$padj <- p.adjust(strand_bias_summary$p_value_binom, method = "BH")

strand_bias_relative <- dplyr::bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  s <- sequences_df[sequences_df$genotype == grid$genotype[i] &
    sequences_df$region_type == grid$region_type[i] &
    sequences_df$window_size == grid$window_size[i] & !is.na(sequences_df$strand_class), ]
  cw <- sum(s$strand_class == "with_strand"); ca <- sum(s$strand_class == "against_strand")
  nt <- cw + ca
  data.frame(grid[i, ], count_with = cw, count_against = ca, n_total = nt,
    pct_with = if (nt > 0) 100 * cw / nt else NA,
    p_value_binom = if (nt > 0) binom.test(cw, nt, 0.5)$p.value else NA)
}))
strand_bias_relative$padj <- p.adjust(strand_bias_relative$p_value_binom, method = "BH")

g4_score_stats <- sequences_df %>%
  dplyr::filter(g4_detected_strand %in% c("+", "-"), !is.na(g4_score)) %>%
  dplyr::group_by(genotype, region_type, window_size, strand = g4_detected_strand) %>%
  dplyr::summarise(n = dplyr::n(), mean_g4 = mean(g4_score), median_g4 = median(g4_score),
    sd_g4 = sd(g4_score), q25 = quantile(g4_score, .25), q75 = quantile(g4_score, .75),
    .groups = "drop")

g4_detection_by_strand <- sequences_df %>%
  dplyr::filter(g4_detected_strand %in% c("+", "-")) %>%
  dplyr::group_by(genotype, region_type, window_size, strand = g4_detected_strand) %>%
  dplyr::summarise(n_total = dplyr::n(), n_with_motif = sum(has_g4_motif, na.rm = TRUE),
    pct_with_motif = 100 * mean(has_g4_motif, na.rm = TRUE), .groups = "drop")

readr::write_tsv(strand_bias_summary,   "results/tables/strand_bias_summary.tsv")
readr::write_tsv(strand_bias_relative,  "results/tables/strand_bias_relative.tsv")
readr::write_tsv(g4_detection_by_strand,"results/tables/g4_detection_by_strand.tsv")
readr::write_tsv(g4_score_stats,        "results/tables/g4_score_stats.tsv")

# --- Comprehensive strand-bias figure (multi-panel) -------------------------
gl <- function(df) { df$genotype <- factor(df$genotype, levels = genotype_order); df }
pA <- ggplot(gl(strand_bias_summary), aes(genotype, pct_plus, fill = genotype)) +
  geom_col() + geom_hline(yintercept = 50, linetype = "dashed") +
  geom_text(aes(label = ifelse(padj < 0.05, "*", "")), vjust = -0.2) +
  facet_grid(region_type ~ window_size) + scale_fill_manual(values = condition_colours) +
  labs(title = "A. Genomic-strand G4 bias (% + strand)", y = "% + strand", x = NULL) +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
pB <- ggplot(gl(strand_bias_relative), aes(genotype, pct_with, fill = genotype)) +
  geom_col() + geom_hline(yintercept = 50, linetype = "dashed") +
  facet_grid(region_type ~ window_size) + scale_fill_manual(values = condition_colours) +
  labs(title = "B. Gene-relative bias (% with gene strand)", y = "% with-strand", x = NULL) +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
pC <- ggplot(gl(dplyr::filter(sequences_df, g4_detected_strand %in% c("+", "-"), window_size == 250)),
             aes(genotype, g4_score, fill = g4_detected_strand)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + facet_wrap(~ region_type) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "C. G4Hunter score by strand (250 bp)", y = "G4Hunter score",
       x = NULL, fill = "strand") + theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pD <- ggplot(gl(g4_detection_by_strand), aes(genotype, pct_with_motif, fill = strand)) +
  geom_col(position = "dodge") + facet_grid(region_type ~ window_size) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "D. G4 detection rate by strand", y = "% with motif", x = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_comp <- cowplot::plot_grid(pA, pB, pC, pD, ncol = 2, labels = NULL)
save_plot(p_comp, "strand_bias_comprehensive", "22_strand", fig_root, width = 16, height = 12)

# --- Strand + template/coding metaprofiles ----------------------------------
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotype_order, assay = "G4_BG4")
bw_by_genotype <- split(bw_meta$filepath, as.character(bw_meta$genotype))

# start-anchored, strand-flipping profile (faithful to strand RMD s7)
prof_matrix <- function(regions, bw_fp) {
  cov <- read_bigwig_rle(bw_fp)
  shared <- intersect(as.character(GenomicRanges::seqnames(regions)), names(cov))
  if (!length(shared)) return(matrix(NA_real_, length(regions), n_bins))
  bm <- matrix(NA_real_, length(regions), n_bins)
  for (chr in shared) {
    ci <- which(as.character(GenomicRanges::seqnames(regions)) == chr)
    if (!length(ci)) next
    chr_len <- length(cov[[chr]])
    st <- pmax(1L, GenomicRanges::start(regions[ci]) - half_width)
    en <- pmin(chr_len, st + 2L * half_width)
    vw <- IRanges::Views(cov[[chr]], start = st, end = en)
    sm <- do.call(rbind, lapply(vw, function(v) { x <- as.numeric(v)
      if (length(x) < 2 * half_width + 1) x <- c(x, rep(NA_real_, 2 * half_width + 1 - length(x))); x }))
    bi <- cut(seq_len(ncol(sm)), breaks = n_bins, labels = FALSE)
    for (b in seq_len(n_bins)) bm[ci, b] <- rowMeans(sm[, which(bi == b), drop = FALSE], na.rm = TRUE)
    minus <- as.character(GenomicRanges::strand(regions[ci])) == "-"
    if (any(minus)) bm[ci[minus], ] <- bm[ci[minus], seq(n_bins, 1L), drop = FALSE]
  }
  bm
}
prof_ci <- function(regions, bws) {
  if (length(regions) == 0) return(list(mean = rep(NA, n_bins), lower = rep(NA, n_bins), upper = rep(NA, n_bins)))
  rm <- do.call(rbind, lapply(bws, function(fp) colMeans(prof_matrix(regions, fp), na.rm = TRUE)))
  m <- colMeans(rm, na.rm = TRUE)
  sem <- apply(rm, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  z <- qnorm(0.975); list(mean = m, lower = m - z * sem, upper = m + z * sem)
}
mk_gr <- function(df) GenomicRanges::GRanges(df$seqname,
  IRanges::IRanges(df$peak_start, df$peak_end),
  strand = ifelse(df$g4_detected_strand %in% c("+", "-"), df$g4_detected_strand, "*"))

collect <- function(split_col, levels) {
  out <- list()
  for (g in genotype_order) for (rt in region_types) {
    sub <- sequences_df[sequences_df$genotype == g & sequences_df$region_type == rt &
                        sequences_df$window_size == 250, ]
    for (lv in levels) {
      d <- sub[sub[[split_col]] == lv, ]
      pr <- prof_ci(mk_gr(d), bw_by_genotype[[g]])
      out[[length(out) + 1]] <- data.frame(position_bin = seq_len(n_bins), genotype = g,
        region_type = rt, level = lv, signal = pr$mean, lower = pr$lower, upper = pr$upper,
        n = nrow(d))
    }
    evict_bigwigs(bw_by_genotype[[g]])
  }
  dplyr::bind_rows(out)
}
strand_prof <- collect("g4_detected_strand", c("+", "-"))
strand_prof$genotype <- factor(strand_prof$genotype, levels = genotype_order)
p_meta <- ggplot(dplyr::filter(strand_prof, !is.na(signal)),
    aes(position_bin, signal, colour = level, fill = level, linetype = level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + facet_grid(region_type ~ genotype) +
  scale_colour_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "G4 CUT&Tag signal by detected strand (95% CI)",
       x = "Position (bin, +/-2kb)", y = "Mean signal", colour = "strand",
       fill = "strand", linetype = "strand") + theme_pub()
save_plot(p_meta, "metaprofile_by_strand", "22_strand", fig_root, width = 14, height = 8)

p_meanbar <- strand_prof %>% dplyr::filter(!is.na(signal)) %>%
  dplyr::group_by(genotype, region_type, level) %>%
  dplyr::summarise(mean_signal = mean(signal, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(level, mean_signal, fill = level)) + geom_col() +
  facet_grid(region_type ~ genotype) +
  scale_fill_manual(values = c("+" = "#E41A1C", "-" = "#377EB8")) +
  labs(title = "Mean G4 signal by strand", x = "strand", y = "Mean signal") + theme_pub()
save_plot(p_meanbar, "mean_signal_by_strand", "22_strand", fig_root, width = 12, height = 8)

tc_prof <- collect("g4_strand_class", c("template", "coding"))
tc_prof$genotype <- factor(tc_prof$genotype, levels = genotype_order)
p_tc <- ggplot(dplyr::filter(tc_prof, !is.na(signal)),
    aes(position_bin, signal, colour = level, fill = level)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.8) + facet_grid(region_type ~ genotype) +
  scale_colour_manual(values = c(template = "#1B9E77", coding = "#D95F02")) +
  scale_fill_manual(values = c(template = "#1B9E77", coding = "#D95F02")) +
  labs(title = "G4 signal: template vs coding strand (95% CI)",
       x = "Position (bin, +/-2kb)", y = "Mean signal", colour = "strand type",
       fill = "strand type") + theme_pub()
save_plot(p_tc, "metaprofile_template_vs_coding", "22_strand", fig_root, width = 14, height = 10)
message("Strand asymmetry analysis done.")
