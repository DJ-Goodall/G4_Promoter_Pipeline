#!/usr/bin/env Rscript
# ============================================================================
# 19_unclassified_peak_context.R   (env: r_g4 / system R)   --- Phase 8 ---
#
# WHERE are the no-canonical-PQS ("unclassified") G4 peaks, and WHAT KIND of DNA
# are they? Direct no-PQS-vs-classified contrasts:
#   (1) Genomic feature (promoter / 5'UTR / intron1 / gene-body / enhancer /
#       intergenic): fraction by class + log2 odds ratio + per-feature Fisher.
#   (2) RepeatMasker class (SINE/LINE/LTR/DNA/Simple_repeat/Satellite/...):
#       same enrichment test. PURPOSE -- intergenic/low-complexity DNA is where
#       CUT&Tag most often gives non-G4 signal (multi-mapping at high-copy
#       repeats, pA-Tn5/antibody stickiness). If the no-PQS peaks are dominated
#       by a repeat class, they are best explained as repeat/accessibility
#       artifacts rather than genuine (non-canonical) G4 binding; if not, they
#       point to a real alternative entity.
#   (3) Low-complexity: dinucleotide Shannon entropy + longest homopolymer run,
#       no-PQS vs classified vs GC-matched background.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv,
#           cache_v2 feature region RDS (config feature.region_files),
#           cache/peak_seqs_201bp.fa, cache/background_seqs.fa,
#           mm10 RepeatMasker (downloaded -> cache/rmsk_mm10.rds)
# Outputs:  results/tables/unclassified_{feature,repeat}_enrichment.csv,
#           results/tables/unclassified_lowcomplexity.csv,
#           results/figures/15_unclassified/unclassified_peak_{feature,repeat,lowcomplexity}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(Biostrings); library(dplyr); library(readr); library(ggplot2)
  library(stringr); library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/15_unclassified", "cache", "data")
fig_root <- "results/figures"

class_cols <- c("No canonical PQS" = "#E41A1C", "Classified (P/AP/H)" = "#377EB8")
group_cols <- c("No canonical PQS" = "#E41A1C", "Classified (P/AP/H)" = "#377EB8",
                "GC-matched background" = "grey55")

# --- Peaks -> 2-way class ---------------------------------------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
union$class2 <- ifelse(union$topology == "no_canonical_PQS", "unclassified", "classified")
GenomeInfoDb::seqlevelsStyle(union) <- "UCSC"
union <- GenomeInfoDb::keepSeqlevels(
  union, intersect(GenomeInfoDb::seqlevels(union), cfg$std_chroms), pruning.mode = "coarse")
message("Peak classes:"); print(table(union$class2))

# --- Enrichment helper: per-category 2x2 Fisher (no-PQS vs classified) ------
enrich_by <- function(cat_vec, class_vec) {
  cats <- sort(unique(cat_vec))
  rows <- lapply(cats, function(ct) {
    a <- sum(class_vec == "unclassified" & cat_vec == ct)
    b <- sum(class_vec == "unclassified" & cat_vec != ct)
    cc <- sum(class_vec == "classified" & cat_vec == ct)
    d <- sum(class_vec == "classified" & cat_vec != ct)
    ft <- fisher.test(matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE))
    data.frame(category = ct, n_unclassified = a, n_classified = cc,
               frac_unclassified = a / (a + b), frac_classified = cc / (cc + d),
               log2_OR = log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (cc + 0.5))),
               p = ft$p.value, stringsAsFactors = FALSE)
  })
  out <- dplyr::bind_rows(rows)
  out$padj <- p.adjust(out$p, "BH")
  out
}

# --- Two-panel enrichment figure (fraction-by-class + log2 OR) --------------
plot_enrich <- function(df, title, order_by_or = TRUE) {
  if (order_by_or) df <- df[order(df$log2_OR), ]
  df$category <- factor(df$category, levels = df$category)
  long <- dplyr::bind_rows(
    data.frame(category = df$category, class = "No canonical PQS",    frac = df$frac_unclassified),
    data.frame(category = df$category, class = "Classified (P/AP/H)", frac = df$frac_classified))
  long$class <- factor(long$class, levels = names(class_cols))
  p1 <- ggplot(long, aes(category, frac, fill = class)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = class_cols, name = NULL) +
    labs(x = NULL, y = "Fraction of peaks", title = title) +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                        legend.position = "top")
  df$sig <- ifelse(df$padj < 0.05, "*", "")
  p2 <- ggplot(df, aes(category, log2_OR, fill = log2_OR > 0)) +
    geom_col(width = 0.7) + geom_hline(yintercept = 0, colour = "grey40") +
    geom_text(aes(label = sprintf("%.2f%s", log2_OR, sig)),
              vjust = ifelse(df$log2_OR > 0, -0.3, 1.2), size = 2.5) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), guide = "none") +
    labs(x = NULL, y = "log2 OR (no-PQS vs classified)") +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 0.95))
}

# ============================================================================
# (1) Genomic feature enrichment
# ============================================================================
region_files <- cfg$feature$region_files
region_list <- lapply(names(region_files), function(nm) {
  gr <- readRDS(file.path(cfg$paths$cache_v2, region_files[[nm]]))
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"; gr
})
names(region_list) <- names(region_files)
union <- annotate_peak_regions(union, region_list, intergenic_label = "intergenic")
feat <- enrich_by(as.character(S4Vectors::mcols(union)$region), union$class2)
readr::write_csv(feat, "results/tables/unclassified_feature_enrichment.csv")
save_plot(plot_enrich(feat, "Genomic feature: no canonical PQS vs classified"),
          "unclassified_peak_feature", "15_unclassified", fig_root, width = 8, height = 8)
message("Feature enrichment:"); print(feat[, c("category", "frac_unclassified", "frac_classified", "log2_OR", "padj")])

# ============================================================================
# (2) RepeatMasker class enrichment
# ============================================================================
rmsk_rds <- "cache/rmsk_mm10.rds"
if (!file.exists(rmsk_rds)) {
  rmsk_gz <- "data/rmsk_mm10.txt.gz"
  if (!file.exists(rmsk_gz)) {
    # Try the configured URL first, then current UCSC mirrors (the legacy
    # hgdownload.soc.ucsc.edu host no longer resolves on some networks).
    urls <- unique(c(cfg$unclassified$rmsk_url,
      "https://hgdownload.gi.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz",
      "https://hgdownload2.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz",
      "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"))
    urls <- urls[!is.na(urls) & nzchar(urls)]
    ok <- FALSE
    for (url in urls) {
      message("Downloading mm10 RepeatMasker (~140 MB, one-time): ", url)
      ok <- tryCatch({
        utils::download.file(url, rmsk_gz, mode = "wb")
        file.exists(rmsk_gz) && file.info(rmsk_gz)$size > 1e6
      }, error = function(e) { message("  failed: ", conditionMessage(e)); FALSE })
      if (ok) break
    }
    if (!ok) stop("Could not download mm10 RepeatMasker from any mirror; ",
                  "download it manually to ", rmsk_gz)
  }
  message("Parsing RepeatMasker ...")
  rm_cols <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName",
               "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass",
               "repFamily", "repStart", "repEnd", "repLeft", "id")
  rm <- readr::read_tsv(rmsk_gz, col_names = rm_cols, show_col_types = FALSE,
    col_types = readr::cols_only(genoName = "c", genoStart = "i", genoEnd = "i",
                                 repClass = "c", repFamily = "c"))
  rmg <- GenomicRanges::GRanges(rm$genoName,
           IRanges::IRanges(rm$genoStart + 1L, rm$genoEnd),   # UCSC 0-based -> 1-based
           repClass = rm$repClass, repFamily = rm$repFamily)
  GenomeInfoDb::seqlevelsStyle(rmg) <- "UCSC"
  rmg <- GenomeInfoDb::keepSeqlevels(
    rmg, intersect(GenomeInfoDb::seqlevels(rmg), cfg$std_chroms), pruning.mode = "coarse")
  saveRDS(rmg, rmsk_rds)
} else rmg <- readRDS(rmsk_rds)

# Per peak: repeat class with the largest overlap (else "none").
rc <- rep("none", length(union))
ov <- GenomicRanges::findOverlaps(union, rmg, ignore.strand = TRUE)
if (length(ov) > 0) {
  ow <- BiocGenerics::width(GenomicRanges::pintersect(
    union[S4Vectors::queryHits(ov)], rmg[S4Vectors::subjectHits(ov)], ignore.strand = TRUE))
  o  <- order(S4Vectors::queryHits(ov), -ow)
  qo <- S4Vectors::queryHits(ov)[o]; so <- S4Vectors::subjectHits(ov)[o]
  f  <- !duplicated(qo)
  rc[qo[f]] <- rmg$repClass[so[f]]
}
# Lump rare classes for a readable figure (full detail stays via the test set).
main_rep <- c("none", "SINE", "LINE", "LTR", "DNA", "Simple_repeat",
              "Low_complexity", "Satellite")
rc_lumped <- ifelse(rc %in% main_rep, rc, "Other")
union$repClass <- rc_lumped
rep_enr <- enrich_by(rc_lumped, union$class2)
readr::write_csv(rep_enr, "results/tables/unclassified_repeat_enrichment.csv")
save_plot(plot_enrich(rep_enr, "RepeatMasker class: no canonical PQS vs classified"),
          "unclassified_peak_repeat", "15_unclassified", fig_root, width = 9, height = 8)
message("Repeat enrichment:"); print(rep_enr[, c("category", "frac_unclassified", "frac_classified", "log2_OR", "padj")])

# ============================================================================
# (3) Low-complexity: dinucleotide entropy + longest homopolymer run
# ============================================================================
peak_seqs <- Biostrings::readDNAStringSet("cache/peak_seqs_201bp.fa")
bg_seqs   <- Biostrings::readDNAStringSet("cache/background_seqs.fa")
peak_cl   <- union$class2[match(names(peak_seqs), names(union))]

lc_metrics <- function(dss) {
  of  <- Biostrings::oligonucleotideFrequency(dss, width = 2)
  p   <- of / pmax(rowSums(of), 1); p[p == 0] <- NA
  ent <- -rowSums(p * log2(p), na.rm = TRUE)
  s   <- as.character(dss)
  maxhp <- vapply(stringr::str_extract_all(s, "A+|C+|G+|T+"),
                  function(r) if (length(r)) max(nchar(r)) else 0L, integer(1))
  data.frame(dinuc_entropy = ent, max_homopolymer = maxhp)
}
set.seed(7)
cap <- function(idx, n) if (length(idx) > n) sample(idx, n) else idx
lc <- dplyr::bind_rows(
  cbind(group = "No canonical PQS",
        lc_metrics(peak_seqs[cap(which(peak_cl == "unclassified"), 8000)])),
  cbind(group = "Classified (P/AP/H)",
        lc_metrics(peak_seqs[cap(which(peak_cl == "classified"), 8000)])),
  cbind(group = "GC-matched background",
        lc_metrics(bg_seqs[cap(seq_along(bg_seqs), 8000)])))
lc$group <- factor(lc$group, levels = names(group_cols))

lc_summ <- lc %>% dplyr::group_by(group) %>% dplyr::summarise(
  n = dplyr::n(),
  dinuc_entropy_med   = round(median(dinuc_entropy), 3),
  max_homopolymer_med = median(max_homopolymer), .groups = "drop")
readr::write_csv(lc_summ, "results/tables/unclassified_lowcomplexity.csv")
message("Low-complexity summary:"); print(as.data.frame(lc_summ))

lc_long <- dplyr::bind_rows(
  data.frame(group = lc$group, metric = "Dinucleotide entropy (bits)", value = lc$dinuc_entropy),
  data.frame(group = lc$group, metric = "Longest homopolymer run (bp)", value = lc$max_homopolymer))
p_lc <- ggplot(lc_long, aes(group, value, fill = group)) +
  geom_violin(scale = "width", alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(x = NULL, y = "Value",
       title = "Low-complexity content: are unclassified peaks repetitive?",
       subtitle = "Lower dinucleotide entropy / longer homopolymers = more low-complexity") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_lc, "unclassified_peak_lowcomplexity", "15_unclassified", fig_root,
          width = 9, height = 4.5)

message("Done. Unclassified-peak genomic context (feature + repeats + low-complexity).")
