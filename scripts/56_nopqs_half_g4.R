#!/usr/bin/env Rscript
# ============================================================================
# 56_nopqs_half_g4.R   (Stage O, rule nopqs_half_g4)
#
# Are the no-canonical-PQS BG4 peaks hybrid / intermolecular G4s?
#
# A DNA:RNA hybrid G4 needs only two G-tracts on the NON-TEMPLATE DNA strand;
# the transcript supplies the rest (Zheng et al. 2013 NAR). A bimolecular DNA
# G4 needs two G-tracts per strand. So if these peaks are real G4 sites that
# pqsfinder cannot see because the quadruplex is only HALF encoded in the
# local DNA, they should carry half sites -- and those half sites should behave
# like hybrid-G4 substrates. This rule tests that, in order:
#
#   A. Relaxed pqsfinder: are there weak or long-loop full G4s that the
#      pipeline's score >= 20 / max_len 50 simply missed?
#   B. Half-site rates (scripts/_shared/half_g4.R): T3 / T2 / T2gg per kb and
#      per peak, against the peaks' OWN dinucleotide shuffle (the primary null:
#      same composition), a width+GC-matched random genome, and PQS peaks.
#   C. Strand vs transcription: in genic peaks, are half sites on the
#      non-template (coding) strand, and only in expressed genes?
#   D. R-loop co-occupancy: do half-site-bearing peaks overlap S9.6 peaks more?
#   E. Signal strength: do the STRONGEST no-PQS peaks carry more half sites?
#   F. Genotype: do the helicase KOs expose half-site peaks WT does not have?
#
# Every rate is reported with a peak-bootstrap CI, and no-PQS vs its shuffle is
# a PAIRED bootstrap (each peak is resampled together with its own shuffled
# twin). The evidence summary states each answer as a number with a CI and a
# mechanical verdict, so a null result reads as a result.
#
# No bigWig is opened. The DeepG4 control is rule 55, deliberately separate so a
# TensorFlow failure cannot take this rule down with it.
#
# Inputs:   cache/nopqs_sets.rds (rule 54), cache/peaks.rds,
#           cache/peaks_Rloop_S96_{union,<genotype>}.rds,
#           results/tables/tss_expression_gene_groups.csv (rule 50)
# Outputs:  results/tables/nopqs_{relaxed_pqs,half_sites,half_site_rates,
#             half_site_strand,half_site_rloop,half_site_strength,
#             half_site_genotype,evidence_summary}.csv
#           results/figures/26_nopqs_half_g4/nopqs_{relaxed_pqs,half_site_rates,
#             half_site_score,half_site_strand,half_site_rloop,
#             half_site_strength,half_site_genotype}.{pdf,png}
# Dev knob: NOPQS_ROOT (see script 54).
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges); library(S4Vectors)
  library(Biostrings); library(pqsfinder); library(dplyr); library(readr)
  library(ggplot2); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})
# NB: the TxDb package attaches AnnotationDbi, which masks dplyr's select().
# Every select() below is therefore written dplyr::select().

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
source(file.path("scripts", "_shared", "half_g4.R"))
cfg <- load_config(config_path)
p   <- cfg$nopqs_half_g4

root     <- Sys.getenv("NOPQS_ROOT", ".")
out      <- function(x) file.path(root, x)
tab      <- function(x) out(file.path("results/tables", x))
fig_root <- out("results/figures")
fig_dir  <- "26_nopqs_half_g4"
ensure_dirs(out("cache"), out("results/tables"), file.path(fig_root, fig_dir))

genotypes <- cfg$genotypes
B         <- as.integer(p$n_boot)
seed      <- p$seed
pqs_min   <- cfg$pqsfinder$min_score          # the pipeline's own PQS threshold

sets_file <- out("cache/nopqs_sets.rds")
sets <- readRDS(sets_file)
meta <- sets$meta; seqs <- sets$seqs
meta$set <- as.character(meta$set)
meta$kb  <- meta$width / 1000
set_levels <- c("noPQS", "noPQS_shuffled", "PQS_sample", "random_gc")
set_labels <- c(noPQS = "No-PQS peaks", noPQS_shuffled = "No-PQS, dinucleotide-shuffled",
                PQS_sample = "PQS peaks (width-matched)", random_gc = "Random genome (width+GC)")
set_cols   <- c(noPQS = "#B2182B", noPQS_shuffled = "grey60",
                PQS_sample = "#2166AC", random_gc = "grey25")

# paired design: noPQS[i] and noPQS_shuffled[i] share a parent. Assert it rather
# than assume it -- every paired CI below depends on it.
i_nop  <- which(meta$set == "noPQS"); i_shuf <- which(meta$set == "noPQS_shuffled")
stopifnot(identical(meta$parent_id[i_nop], meta$parent_id[i_shuf]))
message(sprintf("Sets: %s", paste(sprintf("%s %d", set_levels,
                                          table(factor(meta$set, set_levels))), collapse = " | ")))

# cache keys hash the sets file, so a rebuilt rule-54 set (new seed, new GC
# bins, dev subsample) can never be served a stale scan. Rule 40's cooc_prof_*
# keys omit their parameters and do exactly that; this is the fix, not a copy.
sets_key <- substr(unname(tools::md5sum(sets_file)), 1, 10)

# ============================================================================
# A. Relaxed pqsfinder
# ============================================================================
pqs_cache <- out(sprintf("cache/nopqs_relaxed_pqs_%s_ms%s_ml%s.rds", sets_key,
                         p$relaxed_min_score, p$relaxed_max_len))
pqs_hits <- cache_or_build(pqs_cache, {
  message("Relaxed pqsfinder over ", length(seqs), " sequences (min_score ",
          p$relaxed_min_score, ", max_len ", p$relaxed_max_len, ")...")
  res <- vector("list", length(seqs))
  # pqsfinder writes four progress lines per call to the C-level stderr, which
  # neither capture.output() nor suppressMessages() catches (tried); 9k calls
  # would bury the log. Sinking R's message stream does catch it. The sink is
  # lifted per call so our own progress messages still reach the log.
  quiet <- file(nullfile(), open = "wt")
  for (i in seq_along(seqs)) {
    if (i %% 1000 == 0) message("  pqsfinder ", i, "/", length(seqs))
    sink(quiet, type = "message")
    pv <- tryCatch(
      pqsfinder::pqsfinder(seqs[[i]], strand = "*", min_score = p$relaxed_min_score,
                           max_len = as.integer(p$relaxed_max_len)),
      error = function(e) NULL)
    sink(type = "message")
    if (is.null(pv) || length(pv) == 0) next
    md <- S4Vectors::elementMetadata(pv)
    res[[i]] <- data.frame(seq_id = names(seqs)[i], start = start(pv), end = end(pv),
                           strand = as.character(md$strand), score = md$score,
                           nt = md$nt, stringsAsFactors = FALSE)
  }
  close(quiet)
  hits <- dplyr::bind_rows(res)
  if (!nrow(hits))   # keep the schema so every join below still works
    hits <- data.frame(seq_id = character(0), start = integer(0), end = integer(0),
                       strand = character(0), score = numeric(0), nt = integer(0))
  hits
})
thr_pqs <- as.numeric(unlist(p$pqs_thresholds))

best <- pqs_hits %>% group_by(seq_id) %>%
  summarise(best_score = max(score), best_len = (end - start + 1L)[which.max(score)],
            n_pqs_ge_pipeline = sum(score >= pqs_min), .groups = "drop")
relaxed <- meta %>% dplyr::select(seq_id, set, parent_id) %>%
  left_join(best, by = "seq_id") %>%
  mutate(n_pqs_ge_pipeline = dplyr::coalesce(n_pqs_ge_pipeline, 0L),
         # a hit >= the pipeline threshold that is longer than pqsfinder's
         # default max_len (50) is one the pipeline could not have seen
         long_loop_ge_pipeline = !is.na(best_score) & best_score >= pqs_min & best_len > 50)
readr::write_csv(relaxed, tab("nopqs_relaxed_pqs.csv"))

relaxed_summ <- dplyr::bind_rows(lapply(thr_pqs, function(t) {
  relaxed %>% group_by(set) %>% summarise(
    threshold = t, n = n(), k = sum(!is.na(best_score) & best_score >= t), .groups = "drop")
})) %>% rowwise() %>%
  mutate(frac = k / n,
         lo = binom.test(k, n)$conf.int[1], hi = binom.test(k, n)$conf.int[2]) %>%
  ungroup()
message("Relaxed pqsfinder, fraction of peaks with best PQS >= threshold:")
print(as.data.frame(relaxed_summ %>% dplyr::select(set, threshold, n, k, frac)))

# ============================================================================
# B. Half-site scan
# ============================================================================
scan_one_variant <- function(loop_max) {
  message("Half-site scan, loop_max = ", loop_max, "...")
  s <- half_g4_scan(seqs, min_run = as.integer(p$min_run), gg_min_run = as.integer(p$gg_min_run),
                    loop_min = as.integer(p$loop_min), loop_max = as.integer(loop_max),
                    include_gg = TRUE)
  s$loop_max <- as.integer(loop_max)
  s
}
scan_cache <- out(sprintf("cache/nopqs_half_sites_%s_r%s_gg%s_l%s-%s-%s.rds", sets_key,
                          p$min_run, p$gg_min_run, p$loop_min, p$loop_max,
                          p$loop_max_sensitivity))
sites_all <- cache_or_build(scan_cache, dplyr::bind_rows(
  scan_one_variant(p$loop_max), scan_one_variant(p$loop_max_sensitivity)))
sites_all$set <- meta$set[match(sites_all$seq_id, meta$seq_id)]
sites_all$tier <- factor(sites_all$tier, levels = half_g4_tiers)

# in_pqs: the site overlaps (same strand) a relaxed-scan PQS scoring >= the
# PIPELINE threshold. In noPQS that can only be a long-loop PQS the pipeline's
# max_len hid; in PQS_sample it separates free-standing half sites from mere
# fragments of a full G4 -- the "orphan" class.
if (nrow(pqs_hits)) {
  ph <- pqs_hits[pqs_hits$score >= pqs_min, ]
  hit_gr  <- GRanges(ph$seq_id, IRanges(ph$start, ph$end), strand = ph$strand)
  site_gr <- GRanges(sites_all$seq_id, IRanges(sites_all$start, sites_all$end),
                     strand = sites_all$strand)
  sites_all$in_pqs <- overlapsAny(site_gr, hit_gr)
} else sites_all$in_pqs <- FALSE

# High-confidence threshold per tier, MEASURED on the shuffled null at the main
# loop setting (and separately at the sensitivity setting).
fam <- function(t) ifelse(t %in% c("T4gg", "T2gg"), "GG", "G3")
null_thr <- dplyr::bind_rows(lapply(unique(sites_all$loop_max), function(lm) {
  nul <- sites_all[sites_all$set == "noPQS_shuffled" & sites_all$loop_max == lm, ]
  dplyr::bind_rows(lapply(half_g4_tiers, function(t) {
    x <- nul$half_score[nul$tier == t]; src <- "tier"
    if (length(x) < p$null_min_sites) {
      x <- nul$half_score[fam(as.character(nul$tier)) == fam(t)]
      src <- sprintf("pooled %s family (tier had %d null sites)", fam(t),
                     sum(nul$tier == t))
    }
    data.frame(loop_max = lm, tier = t, n_null_sites = length(x),
               threshold = if (length(x)) unname(quantile(x, p$null_quantile)) else NA_real_,
               threshold_source = src, stringsAsFactors = FALSE)
  }))
}))
sites_all <- sites_all %>%
  left_join(null_thr %>% dplyr::select(loop_max, tier, threshold), by = c("loop_max", "tier")) %>%
  mutate(high_conf = !is.na(threshold) & half_score >= threshold)
message("High-confidence thresholds (", p$null_quantile, " quantile of the shuffled null):")
print(null_thr)

# site catalogue (main loop setting), with genomic coordinates where they exist
cat_sites <- sites_all %>% filter(loop_max == p$loop_max) %>%
  left_join(meta %>% dplyr::select(seq_id, parent_id, chr, pk_start = start, genomic),
            by = "seq_id") %>%
  mutate(g_start = ifelse(genomic, pk_start + start - 1L, NA_integer_),
         g_end   = ifelse(genomic, pk_start + end - 1L, NA_integer_)) %>%
  dplyr::select(set, parent_id, chr, g_start, g_end, strand, tier, n_runs, tetrads,
         run_lens, loops, half_score, threshold, high_conf, in_pqs, offset_start = start)
readr::write_csv(cat_sites, tab("nopqs_half_sites.csv"))

# --- per-peak count columns -------------------------------------------------
count_by <- function(sel) as.integer(table(factor(sites_all$seq_id[sel], levels = meta$seq_id)))
main <- sites_all$loop_max == p$loop_max
combos <- list()
for (t in half_g4_tiers) {
  combos[[paste(t, "all", "main")]]       <- list(t, "all", p$loop_max, main & sites_all$tier == t)
  combos[[paste(t, "high_conf", "main")]] <- list(t, "high_conf", p$loop_max,
                                                  main & sites_all$tier == t & sites_all$high_conf)
  combos[[paste(t, "orphan", "main")]]    <- list(t, "orphan", p$loop_max,
                                                  main & sites_all$tier == t & !sites_all$in_pqs)
  combos[[paste(t, "all", "sens")]]       <- list(t, "all", p$loop_max_sensitivity,
                                                  !main & sites_all$tier == t)
}
cnt <- sapply(combos, function(cb) count_by(cb[[4]]))
# the two families every downstream test uses: G3 half sites (T3 + T2 -- the
# hypothesis; T4 excluded, it is a full G4) and GG half sites (T2gg)
meta$n_g3 <- count_by(main & sites_all$tier %in% c("T3", "T2"))
meta$n_gg <- count_by(main & sites_all$tier == "T2gg")

# --- bootstrap helpers --------------------------------------------------------
boot_mat <- function(n, s) { set.seed(s); matrix(sample.int(n, n * B, replace = TRUE), nrow = n) }
rate_boot <- function(k, kb, M) colSums(matrix(k[M], nrow(M))) / colSums(matrix(kb[M], nrow(M)))
ci <- function(x) unname(quantile(x[is.finite(x)], c(0.025, 0.975), names = FALSE))
safe_ci <- function(x) if (sum(is.finite(x)) < 10) c(NA_real_, NA_real_) else ci(x)

# rate of `k` (per peak) with CI; plus a PAIRED ratio when `k_ref` is given on
# the SAME peaks (noPQS vs its shuffle) and an UNPAIRED ratio for `k_other`.
rate_row <- function(k, kb, s, k_ref = NULL, k_other = NULL, kb_other = NULL) {
  n <- length(k)
  r <- data.frame(n_peaks = n, total_kb = sum(kb), n_sites = sum(k),
                  per_kb = if (n) sum(k) / sum(kb) else NA_real_,
                  per_kb_lo = NA_real_, per_kb_hi = NA_real_,
                  frac_peaks = if (n) mean(k > 0) else NA_real_,
                  frac_lo = NA_real_, frac_hi = NA_real_,
                  ratio_vs_shuffled = NA_real_, ratio_shuf_lo = NA_real_, ratio_shuf_hi = NA_real_,
                  ratio_vs_random_gc = NA_real_, ratio_rnd_lo = NA_real_, ratio_rnd_hi = NA_real_)
  if (n == 0) return(r)
  M  <- boot_mat(n, s)
  bk <- rate_boot(k, kb, M)
  r[c("per_kb_lo", "per_kb_hi")] <- safe_ci(bk)
  bt <- binom.test(sum(k > 0), n)$conf.int
  r[c("frac_lo", "frac_hi")] <- bt
  if (!is.null(k_ref)) {
    r$ratio_vs_shuffled <- (sum(k) / sum(kb)) / (sum(k_ref) / sum(kb))
    r[c("ratio_shuf_lo", "ratio_shuf_hi")] <- safe_ci(bk / rate_boot(k_ref, kb, M))
  }
  if (!is.null(k_other) && length(k_other)) {
    r$ratio_vs_random_gc <- (sum(k) / sum(kb)) / (sum(k_other) / sum(kb_other))
    Mo <- boot_mat(length(k_other), s + 1L)
    r[c("ratio_rnd_lo", "ratio_rnd_hi")] <- safe_ci(bk / rate_boot(k_other, kb_other, Mo))
  }
  r
}

# ============================================================================
# B. rates table
# ============================================================================
i_rnd <- which(meta$set == "random_gc")
rates <- dplyr::bind_rows(lapply(names(combos), function(nm) {
  cb <- combos[[nm]]; k <- cnt[, nm]
  dplyr::bind_rows(lapply(set_levels, function(st) {
    ii <- which(meta$set == st)
    rr <- rate_row(k[ii], meta$kb[ii], seed,
                   k_ref   = if (st == "noPQS") k[i_shuf] else NULL,
                   k_other = if (st %in% c("noPQS", "PQS_sample")) k[i_rnd] else NULL,
                   kb_other = meta$kb[i_rnd])
    cbind(data.frame(set = st, tier = cb[[1]], class = cb[[2]], loop_max = cb[[3]],
                     stringsAsFactors = FALSE), rr)
  }))
})) %>%
  left_join(null_thr %>% dplyr::select(loop_max, tier, threshold, threshold_source),
            by = c("loop_max", "tier")) %>%
  mutate(threshold = ifelse(class == "high_conf", threshold, NA_real_),
         threshold_source = ifelse(class == "high_conf", threshold_source, NA_character_))
readr::write_csv(rates, tab("nopqs_half_site_rates.csv"))
message("Half-site rates (main loop setting, all sites):")
print(as.data.frame(rates %>% filter(class == "all", loop_max == p$loop_max) %>%
                      dplyr::select(set, tier, n_sites, per_kb, frac_peaks, ratio_vs_shuffled,
                             ratio_shuf_lo, ratio_shuf_hi)), digits = 3)

# ============================================================================
# C. strand vs transcription
# ============================================================================
# Gene orientation is assigned at the PEAK (~1 kb): a peak inside genes on one
# strand only gets that strand; both strands -> excluded ("ambiguous"); none ->
# intergenic. Shuffled sequences take their parent's gene, which is the point:
# a shuffle keeps the parent's G/C skew, so it asks "coding-strand bias beyond
# what composition alone produces?".
genes <- suppressMessages(GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.knownGene))
grp <- readr::read_csv("results/tables/tss_expression_gene_groups.csv",
                       col_types = readr::cols(entrez = readr::col_character(),
                                               .default = readr::col_guess()))
gcol <- paste0("group_", p$expression_genotype)
if (!gcol %in% names(grp)) stop("tss_expression_gene_groups.csv has no column ", gcol)
expr_lab <- setNames(ifelse(grepl("^Q[0-9]+$", grp[[gcol]]), "expressed",
                            ifelse(grp[[gcol]] %in% "silent", "silent", "unknown")),
                     grp$entrez)

peak_gr <- GRanges(meta$chr, IRanges(meta$start, meta$end))
ov <- findOverlaps(peak_gr, genes)
ov_df <- data.frame(q = queryHits(ov), gstrand = as.character(strand(genes))[subjectHits(ov)],
                    entrez = names(genes)[subjectHits(ov)], stringsAsFactors = FALSE)
ov_df$expr <- unname(expr_lab[ov_df$entrez]); ov_df$expr[is.na(ov_df$expr)] <- "unknown"
peak_gene <- ov_df %>% group_by(q) %>% summarise(
  gene_strand = if (n_distinct(gstrand) == 1) gstrand[1] else "ambiguous",
  # any expressed -> expressed; else all silent -> silent; else unknown.
  # silent and unknown are never merged (CLAUDE.md: the bridge is lossy).
  expr = if (any(expr == "expressed")) "expressed" else if (all(expr == "silent")) "silent" else "unknown",
  .groups = "drop")
meta$gene_strand <- "intergenic"; meta$expr <- NA_character_
meta$gene_strand[peak_gene$q] <- peak_gene$gene_strand
meta$expr[peak_gene$q] <- peak_gene$expr

g3gg <- sites_all %>% filter(main, tier %in% c("T3", "T2", "T2gg")) %>%
  mutate(family = ifelse(tier == "T2gg", "GG (T2gg)", "G3 (T3 + T2)")) %>%
  left_join(meta %>% dplyr::select(seq_id, gene_strand, expr), by = "seq_id") %>%
  filter(gene_strand %in% c("+", "-")) %>%
  mutate(coding = strand == gene_strand)

strand_tab <- g3gg %>% group_by(set, family, expr) %>%
  summarise(n_sites = n(), n_coding = sum(coding), .groups = "drop") %>%
  rowwise() %>%
  mutate(frac_coding = n_coding / n_sites,
         lo = binom.test(n_coding, n_sites)$conf.int[1],
         hi = binom.test(n_coding, n_sites)$conf.int[2],
         p_vs_half = binom.test(n_coding, n_sites, 0.5)$p.value) %>% ungroup()
# noPQS vs each null, same family x expression class (Fisher on 2x2)
vs <- function(ref) strand_tab %>% filter(set == "noPQS") %>%
  left_join(strand_tab %>% filter(set == ref) %>% dplyr::select(family, expr, r_n = n_sites, r_k = n_coding),
            by = c("family", "expr")) %>% rowwise() %>%
  mutate(p = if (is.na(r_n)) NA_real_ else
           fisher.test(matrix(c(n_coding, n_sites - n_coding, r_k, r_n - r_k), 2))$p.value) %>%
  ungroup() %>% dplyr::select(family, expr, p)
strand_tab <- strand_tab %>%
  left_join(vs("noPQS_shuffled") %>% rename(p_noPQS_vs_shuffled = p), by = c("family", "expr")) %>%
  left_join(vs("random_gc") %>% rename(p_noPQS_vs_random_gc = p), by = c("family", "expr")) %>%
  mutate(p_noPQS_vs_shuffled = ifelse(set == "noPQS", p_noPQS_vs_shuffled, NA_real_),
         p_noPQS_vs_random_gc = ifelse(set == "noPQS", p_noPQS_vs_random_gc, NA_real_))
peak_ctx <- meta %>% count(set, gene_strand, expr, name = "n_peaks")
readr::write_csv(strand_tab, tab("nopqs_half_site_strand.csv"))
message("Strand vs transcription (fraction of half sites on the non-template strand):")
print(as.data.frame(strand_tab %>% dplyr::select(set, family, expr, n_sites, frac_coding, p_vs_half)), digits = 3)

# ============================================================================
# D. R-loop co-occupancy
# ============================================================================
rl_union <- readRDS("cache/peaks_Rloop_S96_union.rds")
rl_geno  <- lapply(setNames(genotypes, genotypes),
                   function(g) readRDS(sprintf("cache/peaks_Rloop_S96_%s.rds", g)))
meta$rloop_union <- overlapsAny(peak_gr, rl_union)
for (g in genotypes) meta[[paste0("rloop_", g)]] <- overlapsAny(peak_gr, rl_geno[[g]])
meta$has_g3 <- meta$n_g3 > 0

rl_row <- function(d, scope) {
  if (!nrow(d)) return(NULL)
  rl <- d[[if (scope == "union") "rloop_union" else paste0("rloop_", scope)]]
  tb <- table(factor(d$has_g3, c(TRUE, FALSE)), factor(rl, c(TRUE, FALSE)))
  ft <- if (all(rowSums(tb) > 0)) fisher.test(tb) else list(estimate = NA, p.value = NA,
                                                            conf.int = c(NA, NA))
  data.frame(scope = scope,
             n_half = sum(d$has_g3), n_half_rloop = sum(d$has_g3 & rl),
             n_nohalf = sum(!d$has_g3), n_nohalf_rloop = sum(!d$has_g3 & rl),
             frac_rloop_half = if (sum(d$has_g3)) mean(rl[d$has_g3]) else NA_real_,
             frac_rloop_nohalf = if (sum(!d$has_g3)) mean(rl[!d$has_g3]) else NA_real_,
             odds_ratio = unname(ft$estimate), or_lo = ft$conf.int[1], or_hi = ft$conf.int[2],
             p = ft$p.value)
}
rloop_tab <- dplyr::bind_rows(lapply(set_levels, function(st) {
  d <- meta[meta$set == st, ]
  rows <- list(rl_row(d, "union"))
  # per-genotype scope only where the peak has a genotype call (random has none)
  if (st != "random_gc") for (g in genotypes) rows[[g]] <- rl_row(d[d[[paste0("in_", g)]] %in% TRUE, ], g)
  cbind(set = st, dplyr::bind_rows(rows))
}))
readr::write_csv(rloop_tab, tab("nopqs_half_site_rloop.csv"))
message("R-loop co-occupancy (G3 half site + vs -):")
print(as.data.frame(rloop_tab %>% filter(scope == "union")), digits = 3)

# ============================================================================
# E. signal strength
# ============================================================================
# Quartile cuts from the WHOLE union's max_z, so "Q4" means strong among all BG4
# peaks, not merely strong among the weaker no-PQS peaks.
union <- readRDS("cache/peaks.rds")$union
qcut  <- unname(quantile(union$max_z, c(0.25, 0.5, 0.75), na.rm = TRUE))
q_lab <- c("Q1 (weakest)", "Q2", "Q3", "Q4 (strongest)")
meta$maxz_q <- factor(q_lab[findInterval(meta$max_z, qcut) + 1L], levels = q_lab)

strength_tab <- dplyr::bind_rows(lapply(c("n_g3", "n_gg"), function(col) {
  dplyr::bind_rows(lapply(q_lab, function(q) {
    dplyr::bind_rows(lapply(c("noPQS", "noPQS_shuffled", "PQS_sample"), function(st) {
      ii <- which(meta$set == st & meta$maxz_q %in% q)
      k_ref <- if (st == "noPQS") meta[[col]][i_shuf[match(meta$parent_id[ii],
                                                           meta$parent_id[i_shuf])]] else NULL
      cbind(data.frame(family = ifelse(col == "n_g3", "G3 (T3 + T2)", "GG (T2gg)"),
                       maxz_quartile = q, set = st, stringsAsFactors = FALSE),
            rate_row(meta[[col]][ii], meta$kb[ii], seed, k_ref = k_ref)[
              , c("n_peaks", "n_sites", "per_kb", "per_kb_lo", "per_kb_hi", "frac_peaks",
                  "ratio_vs_shuffled", "ratio_shuf_lo", "ratio_shuf_hi")])
    }))
  }))
}))
readr::write_csv(strength_tab, tab("nopqs_half_site_strength.csv"))
# within noPQS: does the half-site EXCESS over the peak's own shuffle track
# signal? Excess, not raw count, so an AT-content trend with signal cannot
# masquerade as a half-site trend.
nopd <- meta[i_nop, ]
excess <- (nopd$n_g3 - meta$n_g3[i_shuf]) / nopd$kb
sp_raw <- suppressWarnings(cor.test(nopd$n_g3 / nopd$kb, nopd$max_z, method = "spearman"))
sp_exc <- suppressWarnings(cor.test(excess, nopd$max_z, method = "spearman"))

# ============================================================================
# F. genotype
# ============================================================================
ko <- setdiff(genotypes, cfg$ref_genotype)
in_mat <- as.matrix(nopd[, paste0("in_", genotypes)]); in_mat[is.na(in_mat)] <- FALSE
subsets <- c(setNames(lapply(genotypes, function(g) in_mat[, paste0("in_", g)]),
                      paste("called in", genotypes)),
             list("KO-only (any KO, not WT)" = rowSums(in_mat[, paste0("in_", ko), drop = FALSE]) > 0 &
                    !in_mat[, paste0("in_", cfg$ref_genotype)],
                  "WT-only" = in_mat[, paste0("in_", cfg$ref_genotype)] &
                    rowSums(in_mat[, paste0("in_", ko), drop = FALSE]) == 0,
                  "all four genotypes" = rowSums(in_mat) == length(genotypes)))
geno_tab <- dplyr::bind_rows(lapply(c("n_g3", "n_gg"), function(col) {
  dplyr::bind_rows(lapply(names(subsets), function(sn) {
    w <- which(subsets[[sn]])
    # n recorded before any guard: an empty subset still reports n_peaks = 0
    cbind(data.frame(family = ifelse(col == "n_g3", "G3 (T3 + T2)", "GG (T2gg)"),
                     subset = sn, stringsAsFactors = FALSE),
          rate_row(nopd[[col]][w], nopd$kb[w], seed, k_ref = meta[[col]][i_shuf][w])[
            , c("n_peaks", "n_sites", "per_kb", "per_kb_lo", "per_kb_hi", "frac_peaks",
                "ratio_vs_shuffled", "ratio_shuf_lo", "ratio_shuf_hi")])
  }))
}))
geno_tab$subset <- factor(geno_tab$subset, levels = names(subsets))
readr::write_csv(geno_tab, tab("nopqs_half_site_genotype.csv"))

# ============================================================================
# Evidence summary -- one row per question, a number, a CI and a verdict
# ============================================================================
verdict_ratio <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "not estimable",
                               ifelse(lo > 1, "enriched", ifelse(hi < 1, "depleted", "no difference")))
# With ZERO observed sites the bootstrap numerator is always 0, the CI collapses
# to [0, 0] and verdict_ratio() would print "depleted" with false precision.
# Say what was actually seen instead.
verdict_n <- function(n_obs, lo, hi) if (isTRUE(n_obs == 0)) "no sites observed" else verdict_ratio(lo, hi)
ev <- list()
add <- function(question, metric, value, lo = NA, hi = NA, p = NA, n = NA, verdict = NA)
  ev[[length(ev) + 1L]] <<- data.frame(question = question, metric = metric, value = value,
                                       lo = lo, hi = hi, p = p, n = n, verdict = verdict,
                                       stringsAsFactors = FALSE)

add("0. Sequence-less peaks", "noPQS peaks dropped: > max_n_frac N (inside mm10 assembly gaps)",
    sets$params$n_gap_dropped %||% NA, n = sets$params$n_gap_dropped %||% NA,
    verdict = "artefact: no sequence to carry a PQS")
for (t in c(pqs_min, min(thr_pqs))) {
  r0 <- relaxed_summ %>% filter(set == "noPQS", threshold == t)
  r1 <- relaxed_summ %>% filter(set == "noPQS_shuffled", threshold == t)
  if (!nrow(r0)) next
  ft <- fisher.test(matrix(c(r0$k, r0$n - r0$k, r1$k, r1$n - r1$k), 2))
  add("A. Missed full G4s (relaxed pqsfinder)",
      sprintf("frac noPQS peaks with a PQS >= %s (max_len %s); shuffled: %.3f",
              t, p$relaxed_max_len, r1$frac),
      r0$frac, r0$lo, r0$hi, ft$p.value, r0$n,
      if (ft$p.value < 0.05) (if (r0$frac > r1$frac) "more than composition" else "fewer than composition")
      else "as composition predicts")
}
add("A. Missed full G4s (relaxed pqsfinder)",
    sprintf("frac noPQS peaks whose PQS >= %s is longer than 50 nt (long loop)", pqs_min),
    mean(relaxed$long_loop_ge_pipeline[relaxed$set == "noPQS"]),
    n = sum(relaxed$set == "noPQS"),
    verdict = sprintf("%d peaks carry a G4 only the pipeline's max_len 50 hid",
                      sum(relaxed$long_loop_ge_pipeline[relaxed$set == "noPQS"])))
for (t in c("T3", "T2", "T2gg")) for (cl in c("all", "high_conf")) {
  r <- rates %>% filter(set == "noPQS", tier == t, class == cl, loop_max == p$loop_max)
  add(sprintf("B. Half-site enrichment (%s, %s)", t, cl),
      "noPQS / own-shuffle, per kb (paired bootstrap)",
      r$ratio_vs_shuffled, r$ratio_shuf_lo, r$ratio_shuf_hi, n = r$n_sites,
      verdict = verdict_n(r$n_sites, r$ratio_shuf_lo, r$ratio_shuf_hi))
  add(sprintf("B. Half-site enrichment (%s, %s)", t, cl),
      "noPQS / random genome (width+GC matched), per kb",
      r$ratio_vs_random_gc, r$ratio_rnd_lo, r$ratio_rnd_hi, n = r$n_sites,
      verdict = verdict_n(r$n_sites, r$ratio_rnd_lo, r$ratio_rnd_hi))
}
r <- rates %>% filter(set == "noPQS", tier == "T2", class == "all", loop_max == p$loop_max_sensitivity)
add("B. Half-site enrichment (T2, all)",
    sprintf("noPQS / own-shuffle, per kb, loop_max %s (sensitivity)", p$loop_max_sensitivity),
    r$ratio_vs_shuffled, r$ratio_shuf_lo, r$ratio_shuf_hi, n = r$n_sites,
    verdict = verdict_ratio(r$ratio_shuf_lo, r$ratio_shuf_hi))
for (fm in c("G3 (T3 + T2)", "GG (T2gg)")) for (ex in c("expressed", "silent")) {
  s <- strand_tab %>% filter(set == "noPQS", family == fm, expr == ex)
  if (!nrow(s)) { add(sprintf("C. Non-template-strand bias (%s, %s genes)", fm, ex),
                      "frac of genic half sites on the coding strand", NA, n = 0,
                      verdict = "no sites"); next }
  add(sprintf("C. Non-template-strand bias (%s, %s genes)", fm, ex),
      sprintf("frac on coding strand (p vs 0.5; vs shuffle p = %.3g)", s$p_noPQS_vs_shuffled),
      s$frac_coding, s$lo, s$hi, s$p_vs_half, s$n_sites,
      if (s$p_vs_half >= 0.05) "no strand bias" else if (s$frac_coding > 0.5)
        "non-template (hybrid-compatible)" else "template strand")
}
rr <- rloop_tab %>% filter(set == "noPQS", scope == "union")
add("D. R-loop co-occupancy", "odds ratio: S9.6 peak overlap, G3 half site + vs -",
    rr$odds_ratio, rr$or_lo, rr$or_hi, rr$p, rr$n_half + rr$n_nohalf,
    verdict_ratio(rr$or_lo, rr$or_hi))
add("E. Signal strength", "Spearman: G3 half sites per kb vs max_z (noPQS)",
    unname(sp_raw$estimate), p = sp_raw$p.value, n = nrow(nopd),
    verdict = if (sp_raw$p.value < 0.05) (if (sp_raw$estimate > 0) "rises with signal" else "falls with signal") else "no trend")
add("E. Signal strength", "Spearman: half-site EXCESS over own shuffle vs max_z (noPQS)",
    unname(sp_exc$estimate), p = sp_exc$p.value, n = nrow(nopd),
    verdict = if (sp_exc$p.value < 0.05) (if (sp_exc$estimate > 0) "rises with signal" else "falls with signal") else "no trend")
s4 <- strength_tab %>% filter(set == "noPQS", family == "G3 (T3 + T2)", maxz_quartile == q_lab[4])
add("E. Signal strength", "strongest quartile (union max_z Q4): noPQS / own-shuffle, G3",
    s4$ratio_vs_shuffled, s4$ratio_shuf_lo, s4$ratio_shuf_hi, n = s4$n_peaks,
    verdict = verdict_ratio(s4$ratio_shuf_lo, s4$ratio_shuf_hi))
gk <- geno_tab %>% filter(family == "G3 (T3 + T2)", subset == "KO-only (any KO, not WT)")
add("F. Genotype", "KO-only noPQS peaks: / own-shuffle, G3",
    gk$ratio_vs_shuffled, gk$ratio_shuf_lo, gk$ratio_shuf_hi, n = gk$n_peaks,
    verdict = verdict_ratio(gk$ratio_shuf_lo, gk$ratio_shuf_hi))
evidence <- dplyr::bind_rows(ev)
readr::write_csv(evidence, tab("nopqs_evidence_summary.csv"))
message("Evidence summary:"); print(as.data.frame(evidence %>% dplyr::select(question, value, lo, hi, p, n, verdict)), digits = 3)

# ============================================================================
# Figures -- no point clouds anywhere, so every PDF stays editable
# ============================================================================
n_sub <- sprintf("n = %s peaks per set (%s assembly-gap peaks dropped)",
                 format(sum(meta$set == "noPQS"), big.mark = ","),
                 sets$params$n_gap_dropped %||% "?")
fset <- function(d) mutate(d, set = factor(set, levels = set_levels))

# A. relaxed pqsfinder
p_rel <- ggplot(fset(relaxed_summ), aes(factor(threshold), frac, fill = set)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.25) +
  scale_fill_manual(values = set_cols, labels = set_labels, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = sprintf("Best relaxed pqsfinder score >= (max_len %s)", p$relaxed_max_len),
       y = "Peaks with a PQS at or above threshold",
       title = "Control A: are there weak or long-loop G4s the pipeline missed?",
       subtitle = sprintf("Pipeline definition: score >= %s at max_len 50. %s", pqs_min, n_sub)) +
  theme_pub() + theme(legend.position = "bottom", legend.direction = "vertical")
save_plot(p_rel, "nopqs_relaxed_pqs", fig_dir, fig_root, width = 8, height = 6)

# B. rates -- the headline
tier_show <- c("T3", "T2", "T2gg", "T4", "T4gg")
rp <- rates %>% filter(loop_max == p$loop_max, tier %in% tier_show) %>% fset() %>%
  mutate(tier = factor(tier, levels = tier_show),
         class = factor(class, levels = c("all", "high_conf", "orphan"),
                        labels = c("All sites", "High-confidence (> null 95th pct)",
                                   "Orphan (not inside a PQS >= 20)")))
# No shuffle ratio on the ORPHAN row: shuffling CREATES full PQSs (48% of
# shuffled noPQS peaks reach score >= 20, measured), which swallow the
# shuffle's own half sites and remove them from its orphan count -- the ratio
# would read as enrichment when it is the null being depleted. Orphan is for
# PQS peaks vs random only. The evidence summary never uses it.
rlab <- rp %>% filter(set == "noPQS") %>%
  mutate(lab = ifelse(grepl("^Orphan", class), "no shuffle ratio (see script)",
               ifelse(is.na(ratio_vs_shuffled), "vs shuffle: n/a",
                      sprintf("vs shuffle %.2fx [%.2f-%.2f]", ratio_vs_shuffled,
                              ratio_shuf_lo, ratio_shuf_hi))))
p_rate <- ggplot(rp, aes(set, per_kb, colour = set)) +
  geom_pointrange(aes(ymin = per_kb_lo, ymax = per_kb_hi), size = 0.35) +
  geom_text(data = rlab, aes(x = 2.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.5, size = 2.6, colour = "grey20") +
  # facet_wrap, not facet_grid: grid frees y per ROW, and the tiers in a row
  # differ ~50x in rate (T2gg vs T3), which would flatten T3 to the baseline
  facet_wrap(class ~ tier, ncol = length(tier_show), scales = "free_y",
             labeller = labeller(tier = half_g4_tier_labels, class = label_wrap_gen(40))) +
  scale_colour_manual(values = set_cols, labels = set_labels, name = NULL) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "Half sites per kb (95% peak-bootstrap CI)",
       title = "Do no-PQS peaks carry half-G4 sites beyond what their composition predicts?",
       subtitle = sprintf("G-runs >= %s (GG tiers >= %s), loops %s-%s nt, both strands. %s",
                          p$min_run, p$gg_min_run, p$loop_min, p$loop_max, n_sub)) +
  theme_pub(10) + theme(legend.position = "bottom", axis.ticks.x = element_blank(),
                        strip.text = element_text(size = 7.5))
save_plot(p_rate, "nopqs_half_site_rates", fig_dir, fig_root, width = 15, height = 9)

# half_score distributions vs the null
sd_ <- sites_all %>% filter(main, tier %in% c("T3", "T2", "T2gg")) %>% fset() %>%
  mutate(tier = factor(as.character(tier), levels = c("T3", "T2", "T2gg")))
thr_df <- null_thr %>% filter(loop_max == p$loop_max, tier %in% c("T3", "T2", "T2gg")) %>%
  mutate(tier = factor(tier, levels = c("T3", "T2", "T2gg")))
p_sc <- ggplot(sd_, aes(half_score, colour = set)) +
  geom_freqpoly(aes(y = after_stat(density)), binwidth = 0.1, linewidth = 0.7) +
  geom_vline(data = thr_df, aes(xintercept = threshold), linetype = "dashed", colour = "grey30") +
  facet_wrap(~ tier, scales = "free_y", labeller = labeller(tier = half_g4_tier_labels)) +
  scale_colour_manual(values = set_cols, labels = set_labels, name = NULL) +
  labs(x = "half_score (G4Hunter over the site, oriented to its G-rich strand)",
       y = "Density of sites",
       title = "Half-site quality: are no-PQS sites better than chance sites?",
       subtitle = sprintf("Dashed = high-confidence threshold, the %s quantile of the shuffled null",
                          p$null_quantile)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_sc, "nopqs_half_site_score", fig_dir, fig_root, width = 12, height = 5)

# C. strand
st_p <- strand_tab %>% fset() %>%
  mutate(expr = factor(expr, levels = c("expressed", "silent", "unknown")))
p_st <- ggplot(st_p, aes(expr, frac_coding, colour = set)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.6), size = 0.35) +
  geom_text(aes(y = 0.02, label = n_sites), position = position_dodge(0.6),
            size = 2.4, show.legend = FALSE) +
  facet_wrap(~ family) +
  scale_colour_manual(values = set_cols, labels = set_labels, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = sprintf("Overlapping gene's expression (%s)", p$expression_genotype),
       y = "Half sites G-rich on the non-template strand",
       title = "Test C: DNA:RNA hybrid G4s need G-tracts on the NON-TEMPLATE strand",
       subtitle = "Genic peaks only (genes on one strand). Numbers = sites. Hybrid prediction: > 50% in expressed genes only") +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_st, "nopqs_half_site_strand", fig_dir, fig_root, width = 10, height = 5.5)

# D. R-loop
rl_p <- dplyr::bind_rows(
  rloop_tab %>% transmute(set, scope, half = "G3 half site present", frac = frac_rloop_half,
                          k = n_half_rloop, n = n_half),
  rloop_tab %>% transmute(set, scope, half = "no G3 half site", frac = frac_rloop_nohalf,
                          k = n_nohalf_rloop, n = n_nohalf)) %>%
  filter(n > 0) %>% rowwise() %>%
  mutate(lo = binom.test(k, n)$conf.int[1], hi = binom.test(k, n)$conf.int[2]) %>% ungroup() %>%
  fset() %>% mutate(scope = factor(scope, levels = c("union", genotypes)))
p_rl <- ggplot(rl_p, aes(set, frac, fill = half)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~ scope, nrow = 1) +
  scale_fill_manual(values = c("G3 half site present" = "#B2182B", "no G3 half site" = "grey70"),
                    name = NULL) +
  scale_x_discrete(labels = c(noPQS = "no-PQS", noPQS_shuffled = "shuffled",
                              PQS_sample = "PQS", random_gc = "random")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Peaks overlapping an S9.6 (R-loop) peak",
       title = "Test D: do half-site-bearing peaks co-occupy R-loops?",
       subtitle = "Union = union R-loop peaks; genotype panels = peaks called in that genotype vs its own R-loop peaks") +
  theme_pub() + theme(legend.position = "bottom",
                      axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_rl, "nopqs_half_site_rloop", fig_dir, fig_root, width = 13, height = 5)

# E. strength
sp_p <- strength_tab %>% fset() %>% mutate(maxz_quartile = factor(maxz_quartile, levels = q_lab))
p_sp <- ggplot(sp_p, aes(maxz_quartile, per_kb, colour = set, group = set)) +
  geom_line(linewidth = 0.7) +
  geom_pointrange(aes(ymin = per_kb_lo, ymax = per_kb_hi), size = 0.3) +
  geom_text(data = sp_p %>% filter(set == "noPQS"), aes(y = 0, label = paste0("n=", n_peaks)),
            vjust = -0.3, size = 2.5, show.legend = FALSE) +
  facet_wrap(~ family, scales = "free_y") +
  scale_colour_manual(values = set_cols, labels = set_labels, name = NULL) +
  labs(x = "Peak strength: quartile of max_z across ALL union BG4 peaks",
       y = "Half sites per kb (95% CI)",
       title = "Test E: do the strongest no-PQS peaks carry more half sites?",
       subtitle = sprintf("Within noPQS, Spearman(G3 excess over own shuffle, max_z) = %.3f, p = %.2g",
                          unname(sp_exc$estimate), sp_exc$p.value)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_sp, "nopqs_half_site_strength", fig_dir, fig_root, width = 10, height = 5)

# F. genotype
gt_p <- dplyr::bind_rows(
  geno_tab %>% transmute(family, subset, set = "noPQS", per_kb, lo = per_kb_lo, hi = per_kb_hi, n_peaks),
  dplyr::bind_rows(lapply(c("n_g3", "n_gg"), function(col) dplyr::bind_rows(lapply(names(subsets), function(sn) {
    w <- which(subsets[[sn]])
    r <- rate_row(meta[[col]][i_shuf][w], nopd$kb[w], seed)
    data.frame(family = ifelse(col == "n_g3", "G3 (T3 + T2)", "GG (T2gg)"), subset = sn,
               set = "noPQS_shuffled", per_kb = r$per_kb, lo = r$per_kb_lo, hi = r$per_kb_hi,
               n_peaks = r$n_peaks)
  }))))) %>% fset() %>% mutate(subset = factor(subset, levels = names(subsets)))
p_gt <- ggplot(gt_p, aes(subset, per_kb, colour = set)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.5), size = 0.3) +
  geom_text(data = gt_p %>% filter(set == "noPQS"), aes(y = 0, label = paste0("n=", n_peaks)),
            vjust = -0.3, size = 2.4, show.legend = FALSE) +
  facet_wrap(~ family, scales = "free_y") +
  scale_colour_manual(values = set_cols, labels = set_labels, name = NULL) +
  labs(x = NULL, y = "Half sites per kb (95% CI)",
       title = "Test F: does helicase loss expose half-site peaks?",
       subtitle = "no-PQS peaks by the genotype(s) they were called in, vs each peak's own shuffle") +
  theme_pub() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_gt, "nopqs_half_site_genotype", fig_dir, fig_root, width = 11, height = 5.5)

message("Done.")
