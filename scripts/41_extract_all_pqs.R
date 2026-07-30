#!/usr/bin/env Rscript
# ============================================================================
# 27_extract_all_pqs.R   (env: r_g4)   --- Phase 13 (all-PQS), foundational ---
#
# Rules 01-12 keep only the single top-scoring pqsfinder PQS per peak. This script
# re-runs pqsfinder on the union G4 peaks and retains EVERY PQS, then COLLAPSES
# overlapping PQS (alternative registers of the same physical G4) down to one
# best-scoring motif per locus via a greedy, score-ordered, non-overlapping scan.
# Each surviving motif gets its genomic coordinates, loop-length heuristic topology,
# and sequence -- the substrate for the motif-resolution analyses (rules 29-32) and
# for per-motif G4ShapePredictor scoring (rule g4sp_topology_all -> 04 reused).
#
# Inputs:   cache/peaks.rds (union G4 peaks)
# Outputs:  cache/motifs_all_base.rds       strand-aware GRanges, one range per motif
#                                           (topology added by rule 28 -> motifs_all.rds)
#           results/tables/motif_catalog.csv
#           results/tables/motif_pqs_for_g4sp.csv  (peak_id = MOTIF id; for rule 04)
#           cache/motif_catalog_n<N>.rds    cache keyed on union size (auto-invalidates)
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
ensure_dirs("cache", "results/tables")

genome     <- BSgenome.Mmusculus.UCSC.mm10
mcfg       <- cfg$motif_analysis
min_score  <- mcfg$min_score %||% cfg$pqsfinder$min_score %||% 20
pqs_strand <- cfg$pqsfinder$strand %||% "*"
collapse   <- isTRUE(mcfg$collapse_overlaps %||% TRUE)
par_max    <- cfg$heuristic$parallel_max_loop %||% 3
anti_min   <- cfg$heuristic$antiparallel_min_long_loops %||% 2

union     <- readRDS("cache/peaks.rds")$union
pk_seqs   <- get_sequences(union, genome)      # + strand peak sequences
pk_chr    <- as.character(GenomicRanges::seqnames(union))
pk_start  <- GenomicRanges::start(union)
pk_id     <- names(union)

# Greedy collapse: keep highest-scoring PQS first, drop any later PQS overlapping an
# already-kept one (within-peak coordinates). Returns indices into the PQS object,
# best-score-first. With collapse = FALSE, returns all indices (score-desc order).
select_motifs <- function(starts, ends, scores, do_collapse) {
  o <- order(scores, decreasing = TRUE)
  if (!do_collapse || length(o) <= 1) return(o)
  kept <- integer(0)
  ks <- numeric(0); ke <- numeric(0)
  for (j in o) {
    if (length(kept) == 0 || !any(starts[j] <= ke & ends[j] >= ks)) {
      kept <- c(kept, j); ks <- c(ks, starts[j]); ke <- c(ke, ends[j])
    }
  }
  kept
}

# --- pqsfinder over every union peak, retaining all (collapsed) motifs ------
motif_cache <- sprintf("cache/motif_catalog_n%d.rds", length(union))
motif_df <- cache_or_build(motif_cache, {
  n <- length(pk_seqs)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    if (i %% 5000 == 0) message("  pqsfinder ", i, "/", n)
    sq <- pk_seqs[[i]]
    pv <- tryCatch(pqsfinder::pqsfinder(sq, strand = pqs_strand, min_score = min_score),
                   error = function(e) NULL)
    if (is.null(pv) || length(pv) == 0) next
    md <- S4Vectors::elementMetadata(pv)
    s_in <- as.integer(GenomicRanges::start(pv)); w_in <- as.integer(GenomicRanges::width(pv))
    e_in <- s_in + w_in - 1L
    sel  <- select_motifs(s_in, e_in, md$score, collapse)

    g_start <- pk_start[i] + s_in[sel] - 1L                 # genomic start of motif
    g_end   <- g_start + w_in[sel] - 1L
    g_cen   <- g_start + (w_in[sel] %/% 2L)
    strd    <- as.character(md$strand[sel]); strd[is.na(strd) | strd == "*"] <- "+"
    seqs    <- Biostrings::extractAt(sq, IRanges::IRanges(start = s_in[sel], width = w_in[sel]))
    rc      <- strd == "-"
    if (any(rc)) seqs[rc] <- Biostrings::reverseComplement(seqs[rc])

    out[[i]] <- data.frame(
      peak_id   = pk_id[i],
      motif_idx = seq_along(sel),
      chr       = pk_chr[i],
      motif_start = as.integer(g_start), motif_end = as.integer(g_end),
      motif_center = as.integer(g_cen), strand = strd,
      score = md$score[sel], nt = md$nt[sel],
      ll1 = md$ll1[sel], ll2 = md$ll2[sel], ll3 = md$ll3[sel],
      pqs_seq = as.character(seqs), stringsAsFactors = FALSE)
  }
  df <- dplyr::bind_rows(out)
  df$motif_id <- sprintf("%s_m%02d", df$peak_id, df$motif_idx)
  df
})

# Loop heuristic per motif + handy loop summaries.
motif_df <- motif_df %>%
  mutate(
    has_pqs = TRUE,
    topology_heuristic = mapply(heuristic_topology, has_pqs, ll1, ll2, ll3,
                               MoreArgs = list(par_max = par_max, anti_min = anti_min)),
    max_loop = pmax(ll1, ll2, ll3), min_loop = pmin(ll1, ll2, ll3),
    total_loop = ll1 + ll2 + ll3) %>%
  select(motif_id, peak_id, motif_idx, chr, motif_start, motif_end, motif_center,
         strand, score, nt, ll1, ll2, ll3, max_loop, min_loop, total_loop,
         topology_heuristic, has_pqs, pqs_seq)

readr::write_csv(motif_df, "results/tables/motif_catalog.csv")

# Table for the (reused, unchanged) G4ShapePredictor script: it keys on `peak_id`
# and passes it through, so we feed the MOTIF id under that column name.
readr::write_csv(
  data.frame(peak_id = motif_df$motif_id, has_pqs = TRUE, pqs_seq = motif_df$pqs_seq,
             stringsAsFactors = FALSE),
  "results/tables/motif_pqs_for_g4sp.csv")

# Strand-aware GRanges of motif spans (topology added by rule 28).
motifs_gr <- GenomicRanges::GRanges(
  seqnames = motif_df$chr,
  ranges   = IRanges::IRanges(motif_df$motif_start, motif_df$motif_end),
  strand   = motif_df$strand)
names(motifs_gr) <- motif_df$motif_id
S4Vectors::mcols(motifs_gr) <- motif_df[, c("motif_id", "peak_id", "motif_center",
                                            "score", "nt", "ll1", "ll2", "ll3",
                                            "topology_heuristic", "pqs_seq")]
saveRDS(motifs_gr, "cache/motifs_all_base.rds")   # rule 28 adds $topology -> motifs_all.rds

n_peaks_with <- length(unique(motif_df$peak_id))
message(sprintf("Done. %d motifs across %d / %d union peaks (%.1f%% with >=1 PQS); mean %.2f motifs/peak.",
                nrow(motif_df), n_peaks_with, length(union),
                100 * n_peaks_with / length(union), nrow(motif_df) / n_peaks_with))
