#!/usr/bin/env Rscript
# ============================================================================
# 54_nopqs_sets.R   (Stage O, rule nopqs_sets)
#
# Build the four sequence sets that rules 55 (DeepG4 scan) and 56 (half-G4
# scan) score, ONCE, so both rules see byte-identical sequences:
#
#   noPQS           the union BG4 peaks rule 21 labelled no_canonical_PQS
#                   (pqsfinder found nothing >= pqsfinder.min_score anywhere in
#                   the full peak). 2,405 of 62,913 in the 2026-08 run.
#   noPQS_shuffled  each noPQS peak dinucleotide-shuffled (universalmotif,
#                   Euler method). The PRIMARY null: same length, composition,
#                   CpG and G/C skew as its parent, no positional structure. A
#                   half site that is merely what AT-rich sequence produces by
#                   chance is equally frequent here.
#   PQS_sample      PQS-bearing union peaks, sampled to the noPQS WIDTH
#                   distribution (width by width), same n. Positive reference.
#   random_gc       random mm10 windows matched to noPQS on width AND GC
#                   (quantile bins, as rule 18 builds its DeepG4 background),
#                   excluding every union peak and any window with an N.
#
# Why these and not rule 18's background: rule 18's windows are 201 bp and
# matched to the GC of ALL peaks (median 0.49). The noPQS peaks are ~1 kb and
# AT-rich (median GC 0.36), so that background would be the wrong composition
# and the wrong length for a per-kb rate.
#
# Shuffled sequences inherit their parent's coordinates and metadata (so the
# strand-vs-gene test in rule 56 can ask "beyond composition?") but carry
# genomic = FALSE: their base positions are not genomic positions.
#
# Inputs:   cache/peaks.rds, results/tables/peak_topology.csv
# Outputs:  cache/nopqs_sets.rds  list(meta = data.frame, seqs = DNAStringSet,
#                                      params = list)
#           results/tables/nopqs_sets_summary.csv
#
# Dev knobs (environment variables; unset in the real run):
#   NOPQS_ROOT   write this stage's own outputs under this directory instead of
#                the repo root -- lets the rule be tested without creating
#                files Snakemake would mistake for finished outputs.
#   NOPQS_DEV_N  keep only this many noPQS peaks (comparators follow).
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges); library(S4Vectors)
  library(Biostrings); library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
p   <- cfg$nopqs_half_g4

root  <- Sys.getenv("NOPQS_ROOT", ".")
dev_n <- suppressWarnings(as.integer(Sys.getenv("NOPQS_DEV_N", "")))
out   <- function(x) file.path(root, x)
ensure_dirs(out("cache"), out("results/tables"))

genome <- BSgenome.Mmusculus.UCSC.mm10
std    <- cfg$std_chroms
seed   <- p$seed

# --- noPQS peaks + metadata --------------------------------------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
union <- GenomeInfoDb::keepSeqlevels(union, intersect(seqlevels(union), std),
                                     pruning.mode = "coarse")
strand(union) <- "*"

is_nop <- union$topology %in% "no_canonical_PQS"
nop    <- union[is_nop]
pqs    <- union[!is_nop & !is.na(union$topology)]
message(sprintf("Union peaks: %d | noPQS: %d | PQS-bearing: %d",
                length(union), length(nop), length(pqs)))
if (!is.na(dev_n) && dev_n < length(nop)) {
  set.seed(seed)
  nop <- nop[sort(sample(length(nop), dev_n))]
  message("DEV MODE: NOPQS_DEV_N = ", dev_n, " -> ", length(nop), " noPQS peaks")
}

# --- drop sequence-less peaks ------------------------------------------------
# Measured 2026-09-24: 148 of the 2,405 noPQS peaks contain N, and every one is
# >= 57% N -- 145 are >= 99% N, i.e. the whole "peak" lies inside an mm10
# assembly gap. They have no sequence, so "no PQS" is trivially true of them and
# they would dilute every per-kb rate towards zero. Only 3 of 60,508 PQS peaks
# are affected. The cut-off is DeepG4's own (N fraction > 0.1), so the set that
# survives is exactly the set rule 19 could score. They are counted in the
# summary table, never silently lost.
nop_nfrac <- as.numeric(letterFrequency(get_sequences(nop, genome), "N", as.prob = TRUE))
n_gap <- sum(nop_nfrac > p$max_n_frac)
message(sprintf("Dropping %d noPQS peaks with N fraction > %.2f (assembly gaps); %d remain",
                n_gap, p$max_n_frac, length(nop) - n_gap))
nop <- nop[nop_nfrac <= p$max_n_frac]
n <- length(nop)

peak_meta <- function(gr, set, genomic = TRUE, id = names(gr)) {
  data.frame(seq_id = paste0(set, "|", id), set = set, parent_id = names(gr),
             chr = as.character(seqnames(gr)), start = start(gr), end = end(gr),
             width = width(gr), genomic = genomic,
             max_z = gr$max_z %||% NA_real_, mean_signal = gr$mean_signal %||% NA_real_,
             in_WT = gr$in_WT %||% NA, in_DHX36KO = gr$in_DHX36KO %||% NA,
             in_FANCJKO = gr$in_FANCJKO %||% NA, in_dKO = gr$in_dKO %||% NA,
             group = gr$group %||% NA_character_, stringsAsFactors = FALSE)
}

# --- 1. noPQS ----------------------------------------------------------------
nop_seq <- get_sequences(nop, genome)
names(nop_seq) <- paste0("noPQS|", names(nop))
meta_nop <- peak_meta(nop, "noPQS")

# --- 2. dinucleotide shuffle of each noPQS peak -----------------------------
# universalmotif's Euler shuffle preserves every k-mer count exactly. It works
# per sequence, so each shuffled twin keeps its parent's length and GC.
message("Shuffling noPQS sequences (k = ", p$shuffle_k, ")...")
shuf_seq <- universalmotif::shuffle_sequences(nop_seq, k = p$shuffle_k,
                                              method = "euler", rng.seed = seed)
names(shuf_seq) <- paste0("noPQS_shuffled|", names(nop))
meta_shuf <- peak_meta(nop, "noPQS_shuffled", genomic = FALSE)
stopifnot(all(width(shuf_seq) == width(nop_seq)))

# --- 3. PQS-bearing peaks, width-matched ------------------------------------
# Width by width: for every width w that k noPQS peaks have, draw k PQS peaks of
# width exactly w. Peaks are multiples of the 500 bp calling bin, so every
# width present in noPQS is well supplied among 60k PQS peaks; if one ever is
# not, the nearest widths fill in and the shortfall is logged.
set.seed(seed + 1L)
w_tab <- table(width(nop))
pick <- unlist(lapply(names(w_tab), function(w) {
  k <- as.integer(w_tab[[w]]); w <- as.integer(w)
  pool <- which(width(pqs) == w)
  if (length(pool) >= k) return(sample(pool, k))
  message(sprintf("  PQS_sample: only %d PQS peaks of width %d for %d wanted; ",
                  length(pool), w, k), "filling from the nearest widths")
  near <- order(abs(width(pqs) - w))
  near <- setdiff(near, pool)[seq_len(k - length(pool))]
  c(pool, near)
}))
pqs_s <- pqs[sort(pick)]
pqs_seq <- get_sequences(pqs_s, genome)
names(pqs_seq) <- paste0("PQS_sample|", names(pqs_s))
meta_pqs <- peak_meta(pqs_s, "PQS_sample")

# --- 4. random genomic windows, width + GC matched --------------------------
# GC is computed over ACGT bases only. 148 noPQS peaks sit against assembly
# gaps (every one of them > 10% N); counting N in the denominator dragged their
# GC down and demanded random windows at a GC the genome barely has -- the
# first build came up 22 windows short for exactly that reason.
gc_acgt <- function(x) {
  lf <- letterFrequency(x, c("G", "C", "A", "T"))
  (lf[, "G"] + lf[, "C"]) / pmax(rowSums(lf), 1)
}
chrom_sizes <- default_chrom_sizes(std)
nop_gc <- gc_acgt(nop_seq)
brks <- unique(quantile(nop_gc, probs = seq(0, 1, length.out = p$gc_bins + 1)))
nb   <- length(brks) - 1L
nop_bin <- findInterval(nop_gc, brks, rightmost.closed = TRUE, all.inside = TRUE)

set.seed(seed + 2L)
rnd_list <- list(); match_rows <- list()
for (w in sort(unique(width(nop)))) {
  need_bin <- table(factor(nop_bin[width(nop) == w], levels = seq_len(nb)))
  # floor per width: rare widths (one 4.5 kb peak) otherwise get a 40-window
  # pool, and a < 0.30 GC window of 2-5 kb is rare -- the first build came up 22
  # short, every miss a 2-5 kb peak in the lowest GC bin.
  pool_n <- max(sum(need_bin) * p$random_pool_mult, p$random_pool_min)
  chr <- sample(names(chrom_sizes), pool_n, replace = TRUE,
                prob = as.numeric(chrom_sizes) / sum(as.numeric(chrom_sizes)))
  st  <- floor(runif(pool_n) * (as.numeric(chrom_sizes[chr]) - w - 1)) + 1
  pool <- GRanges(chr, IRanges(as.integer(st), width = w))
  pool <- pool[!overlapsAny(pool, union)]
  ps   <- getSeq(genome, pool)
  ok   <- as.numeric(letterFrequency(ps, "N")) == 0
  pool <- pool[ok]; ps <- ps[ok]
  pgc  <- gc_acgt(ps)
  # windows outside the noPQS GC range are not a match for any bin
  inr  <- pgc >= brks[1] & pgc <= brks[length(brks)]
  pbin <- findInterval(pgc, brks, rightmost.closed = TRUE, all.inside = TRUE)
  pbin[!inr] <- NA
  sel <- unlist(lapply(seq_len(nb), function(b) {
    avail <- which(pbin == b); k <- as.integer(need_bin[[b]])
    if (k == 0) return(integer(0))
    sample(avail, min(length(avail), k))
  }))
  rnd_list[[as.character(w)]] <- pool[sel]
  match_rows[[as.character(w)]] <- data.frame(
    width = w, gc_bin = seq_len(nb), target_n = as.integer(need_bin),
    realised_n = as.integer(table(factor(pbin[sel], levels = seq_len(nb)))))
}
rnd <- do.call(c, unname(rnd_list))
names(rnd) <- sprintf("rnd_%06d", seq_along(rnd))
rnd_seq <- getSeq(genome, rnd)
names(rnd_seq) <- paste0("random_gc|", names(rnd))
meta_rnd <- peak_meta(rnd, "random_gc")
gc_match <- do.call(rbind, match_rows)
short <- sum(gc_match$target_n - gc_match$realised_n)
if (short > 0)
  warning(sprintf("random_gc is %d windows short of %d requested (GC matching, ",
                  short, n), "see nopqs_sets_summary.csv)")

# --- assemble ---------------------------------------------------------------
seqs <- c(nop_seq, shuf_seq, pqs_seq, rnd_seq)
meta <- dplyr::bind_rows(meta_nop, meta_shuf, meta_pqs, meta_rnd)
stopifnot(identical(meta$seq_id, names(seqs)))
meta$gc     <- gc_acgt(seqs)
meta$n_frac <- as.numeric(letterFrequency(seqs, "N", as.prob = TRUE))
set_levels  <- c("noPQS", "noPQS_shuffled", "PQS_sample", "random_gc")
meta$set    <- factor(meta$set, levels = set_levels)

saveRDS(list(meta = meta, seqs = seqs,
             params = list(seed = seed, shuffle_k = p$shuffle_k, gc_bins = p$gc_bins,
                           gc_breaks = brks, gc_match = gc_match, dev_n = dev_n,
                           n_gap_dropped = n_gap, max_n_frac = p$max_n_frac)),
        out("cache/nopqs_sets.rds"))

summ <- meta %>% group_by(set) %>% summarise(
  n = n(), width_median = median(width), total_kb = sum(width) / 1000,
  gc_median = round(median(gc), 4), gc_q25 = round(quantile(gc, .25), 4),
  gc_q75 = round(quantile(gc, .75), 4),
  n_with_N = sum(n_frac > 0), n_N_over_10pct = sum(n_frac > 0.1),
  max_z_median = round(median(max_z, na.rm = TRUE), 3), .groups = "drop") %>%
  mutate(gc_match_shortfall = ifelse(set == "random_gc", short, NA_integer_),
         n_noPQS_dropped_assembly_gap = ifelse(set == "noPQS", n_gap, NA_integer_))
readr::write_csv(summ, out("results/tables/nopqs_sets_summary.csv"))
message("Sets built:"); print(as.data.frame(summ))
