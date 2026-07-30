#!/usr/bin/env Rscript
# ============================================================================
# 02_extract_sequences.R   (env: r_g4)
#
# (a) Detect the top pqsfinder PQS inside each peak -> the G4-forming subsequence
#     (G4SP input) + loop lengths (topology heuristic input) + its genomic centre.
# (b) Extract 201 bp sequences for DeepG4, CENTRED ON THE PQS (not the peak
#     midpoint): peaks here are 500 bp - 5.5 kb wide, so a midpoint window usually
#     misses the G4 motif. Fall back to the peak centre for peaks with no PQS.
# (c) Build a GC/length-matched random genomic background (DeepG4 AUROC negatives).
#
# Inputs:   cache/peaks.rds
# Outputs:  cache/peak_seqs_201bp.fa
#           cache/background_seqs.fa
#           cache/peak_pqs_n<N>.rds   (cache keyed on peak count, auto-invalidates)
#           results/tables/peak_pqs.csv   (peak_id, has_pqs, pqs_seq, score, nt,
#                                          ll1, ll2, ll3, strand, pqs_center, ...)
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges); library(S4Vectors)
  library(Biostrings); library(pqsfinder); library(dplyr); library(readr)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

genome      <- BSgenome.Mmusculus.UCSC.mm10
std_chroms  <- cfg$std_chroms
width201    <- cfg$sequences$deepg4_width
bg_n        <- cfg$sequences$background_n
bg_seed     <- cfg$sequences$background_seed
gc_bins     <- cfg$sequences$gc_bins
pqs_min     <- cfg$pqsfinder$min_score
pqs_strand  <- cfg$pqsfinder$strand
half201     <- width201 %/% 2L

ensure_dirs("cache", "results/tables")

chrom_sizes <- default_chrom_sizes(std_chroms)

peaks_obj <- readRDS("cache/peaks.rds")
union     <- peaks_obj$union

# --- (a) pqsfinder PQS detection per peak ----------------------------------
# Run pqsfinder on the FULL peak so the strongest G4 motif anywhere in the peak
# is captured. We record its genomic centre so the DeepG4 window (below) can be
# placed on the actual candidate G4 rather than the arbitrary peak midpoint.
# Cache is keyed on peak count: a different peak set (smoke-test subsample vs the
# full run) gets a different file, so a stale PQS table is never silently reused.
peak_full_seqs <- get_sequences(union, genome)
pk_start <- start(union)

pqs_cache <- sprintf("cache/peak_pqs_n%d.rds", length(union))
pqs_tbl <- cache_or_build(pqs_cache, {
  n <- length(peak_full_seqs)
  res <- vector("list", n)
  for (i in seq_len(n)) {
    if (i %% 5000 == 0) message("  pqsfinder ", i, "/", n)
    sq <- peak_full_seqs[[i]]
    pv <- tryCatch(
      pqsfinder::pqsfinder(sq, strand = pqs_strand, min_score = pqs_min),
      error = function(e) NULL)
    if (is.null(pv) || length(pv) == 0) {
      res[[i]] <- data.frame(peak_id = names(peak_full_seqs)[i], has_pqs = FALSE,
                             pqs_seq = NA_character_, score = NA_integer_, nt = NA_integer_,
                             ll1 = NA_integer_, ll2 = NA_integer_, ll3 = NA_integer_,
                             strand = NA_character_, pqs_center = NA_integer_,
                             stringsAsFactors = FALSE)
      next
    }
    md  <- S4Vectors::elementMetadata(pv)
    top <- which.max(md$score)
    strd <- tryCatch(as.character(md$strand[top]), error = function(e) "+")
    if (length(strd) == 0 || is.na(strd)) strd <- "+"
    sub  <- Biostrings::subseq(sq, start = start(pv)[top], width = width(pv)[top])
    if (!is.na(strd) && strd == "-") sub <- Biostrings::reverseComplement(sub)
    # Genomic centre of the top PQS (positions are on the + strand peak sequence).
    pqs_centre <- pk_start[i] + start(pv)[top] - 1L + (width(pv)[top] %/% 2L)
    res[[i]] <- data.frame(
      peak_id = names(peak_full_seqs)[i], has_pqs = TRUE,
      pqs_seq = as.character(sub), score = md$score[top], nt = md$nt[top],
      ll1 = md$ll1[top], ll2 = md$ll2[top], ll3 = md$ll3[top],
      strand = strd, pqs_center = as.integer(pqs_centre), stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(res)
})

pqs_tbl <- pqs_tbl %>%
  mutate(max_loop = pmax(ll1, ll2, ll3),
         min_loop = pmin(ll1, ll2, ll3),
         total_loop = ll1 + ll2 + ll3)
readr::write_csv(pqs_tbl, "results/tables/peak_pqs.csv")

message(sprintf("Done pqsfinder. Peaks with canonical PQS: %d / %d (%.1f%%)",
                sum(pqs_tbl$has_pqs), nrow(pqs_tbl),
                100 * mean(pqs_tbl$has_pqs)))

# --- (b) Peak 201 bp sequences for DeepG4, centred on the PQS ---------------
# Centre on the top PQS where present, else on the peak midpoint. Then drop any
# window that runs past a chromosome end (getSeq() errors on out-of-bounds).
idx        <- match(names(union), pqs_tbl$peak_id)
pqs_centre <- pqs_tbl$pqs_center[idx]
peak_mid   <- start(union) + (width(union) %/% 2L)
centre     <- ifelse(!is.na(pqs_centre), pqs_centre, peak_mid)
n_on_pqs   <- sum(!is.na(pqs_centre))

dg_win <- GRanges(seqnames = seqnames(union),
                  ranges = IRanges(start = as.integer(centre - half201), width = width201),
                  strand = "*")
names(dg_win) <- names(union)

GenomeInfoDb::seqlengths(dg_win) <- chrom_sizes[GenomeInfoDb::seqlevels(dg_win)]
in_bounds <- start(dg_win) >= 1 &
  end(dg_win) <= GenomeInfoDb::seqlengths(dg_win)[as.character(seqnames(dg_win))]
if (any(!in_bounds))
  message("Dropping ", sum(!in_bounds), " DeepG4 windows that exceed chromosome bounds.")
dg_win <- dg_win[in_bounds]

message("Extracting ", length(dg_win), " DeepG4 windows (", width201,
        " bp; ", n_on_pqs, " centred on PQS, ",
        length(union) - n_on_pqs, " on peak midpoint)...")
peak_seqs <- get_sequences(dg_win, genome)
write_fasta(peak_seqs, "cache/peak_seqs_201bp.fa")

# --- (c) GC/length-matched random background -------------------------------
# Matched to the GC distribution of the DeepG4 peak windows above (the set
# DeepG4 actually scores), so the AUROC compares like with like.
message("Building GC-matched background (target n = ", bg_n, ")...")
peak_gc <- as.numeric(Biostrings::letterFrequency(peak_seqs, "GC", as.prob = TRUE))

set.seed(bg_seed)
pool_n <- bg_n * 6L
chr_choice <- sample(names(chrom_sizes), pool_n, replace = TRUE,
                     prob = as.numeric(chrom_sizes) / sum(as.numeric(chrom_sizes)))
max_start  <- as.numeric(chrom_sizes[chr_choice]) - width201 - 1
starts     <- floor(runif(pool_n) * max_start) + 1
pool <- GRanges(seqnames = chr_choice,
                ranges = IRanges(start = as.integer(starts), width = width201),
                strand = "*")
pool <- pool[!overlapsAny(pool, union, ignore.strand = TRUE)]   # exclude real peaks
pool_seqs <- get_sequences(setNames(pool, sprintf("pool_%d", seq_along(pool))), genome)
# drop windows containing N
n_count   <- Biostrings::letterFrequency(pool_seqs, "N")
pool_ok   <- as.numeric(n_count) == 0
pool      <- pool[pool_ok]; pool_seqs <- pool_seqs[pool_ok]
pool_gc   <- as.numeric(Biostrings::letterFrequency(pool_seqs, "GC", as.prob = TRUE))

# GC-quantile bins from the peak distribution; sample background to match shape
brks <- unique(quantile(peak_gc, probs = seq(0, 1, length.out = gc_bins + 1),
                        na.rm = TRUE))
peak_bin <- findInterval(peak_gc, brks, rightmost.closed = TRUE, all.inside = TRUE)
pool_bin <- findInterval(pool_gc, brks, rightmost.closed = TRUE, all.inside = TRUE)
target   <- round(bg_n * (table(factor(peak_bin, levels = seq_len(length(brks) - 1))) /
                          length(peak_bin)))
sel_idx <- unlist(lapply(seq_len(length(brks) - 1), function(b) {
  avail <- which(pool_bin == b)
  take  <- min(length(avail), target[[as.character(b)]])
  if (take > 0) sample(avail, take) else integer(0)
}))
bg_seqs <- pool_seqs[sel_idx]
names(bg_seqs) <- sprintf("bg_%06d", seq_along(bg_seqs))
write_fasta(bg_seqs, "cache/background_seqs.fa")
message("Background sequences written: ", length(bg_seqs))
