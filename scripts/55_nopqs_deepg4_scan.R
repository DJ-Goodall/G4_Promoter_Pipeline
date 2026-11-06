#!/usr/bin/env Rscript
# ============================================================================
# 55_nopqs_deepg4_scan.R   (Stage O, rule nopqs_deepg4_scan; TensorFlow)
#
# Control: are there G4s in the no-PQS peaks that DeepG4 would call if it were
# shown the WHOLE peak?
#
# Why this is needed. Rule 18 cuts ONE 201 bp window per peak -- centred on the
# best pqsfinder PQS, or on the peak MIDPOINT when there is none -- and rule 19
# scores that one window. For a ~1 kb no-PQS peak that is ~19% of its sequence,
# at an arbitrary position, while every PQS peak is scored at its best G4
# candidate. So the existing no-PQS vs PQS contrast (median 0.022 vs 0.125,
# rule 35) mixes window placement with sequence. Here every set is tiled the
# same way, so the comparison is like-for-like.
#
# Method. Each sequence from rule 54 is tiled into 201 bp windows every
# deepg4_step bp (plus one window flush with the 3' end, so no base is unseen).
# DeepG4::DeepG4() is called on the windows directly. DeepG4Scan() is NOT used:
# it discards every window below its threshold (default 0.5), and the
# distribution below 0.5 is the point of a control. Windows with N fraction
# > 0.1 are dropped BEFORE the call, replicating the filter inside DeepG4()
# that would otherwise shorten its output and misalign rows (see rule 19).
# Rule 54 already removed the assembly-gap peaks, so this should drop ~none.
#
# Isolated in its own rule, like rule 19, so a TensorFlow failure costs this
# rule only: rule 56 does not depend on it.
#
# Inputs:   cache/nopqs_sets.rds (rule 54), results/tables/deepg4_scores.csv
#           (rule 19: the old midpoint score, for the placement comparison)
# Outputs:  cache/nopqs_deepg4_windows.rds         per window
#           results/tables/nopqs_deepg4_scan.csv   per peak
#           results/tables/nopqs_deepg4_summary.csv
#           results/figures/26_nopqs_half_g4/nopqs_deepg4_maxprob.{pdf,png}
# Dev knob: NOPQS_ROOT (see script 54).
# ============================================================================

Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", CUDA_VISIBLE_DEVICES = "-1")  # quiet, CPU-only

suppressPackageStartupMessages({
  library(Biostrings); library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
p   <- cfg$nopqs_half_g4

root     <- Sys.getenv("NOPQS_ROOT", ".")
out      <- function(x) file.path(root, x)
fig_root <- out("results/figures")
ensure_dirs(out("cache"), out("results/tables"), file.path(fig_root, "26_nopqs_half_g4"))

win_w <- cfg$sequences$deepg4_width        # 201, DeepG4's fixed input
step  <- p$deepg4_step
cuts  <- as.numeric(unlist(p$deepg4_prob_cuts))

# --- TensorFlow via reticulate (same pattern as rule 19) --------------------
suppressPackageStartupMessages(library(reticulate))
tf_env <- cfg$deepg4$reticulate_condaenv %||% "tf-new"
tryCatch(reticulate::use_condaenv(tf_env, required = TRUE),
         error = function(e) message("use_condaenv('", tf_env, "') failed: ",
                                     conditionMessage(e)))
if (!requireNamespace("DeepG4", quietly = TRUE))
  stop("DeepG4 R package not found. See scripts/19_deepg4_predict.R for install notes.")
suppressPackageStartupMessages(library(DeepG4))

sets <- readRDS(out("cache/nopqs_sets.rds"))
meta <- sets$meta; seqs <- sets$seqs
set_levels <- levels(meta$set)

# --- tile -------------------------------------------------------------------
L <- width(seqs)
starts <- lapply(L, function(l) {
  if (l < win_w) return(integer(0))
  s <- seq.int(1L, l - win_w + 1L, by = step)
  if (tail(s, 1) != l - win_w + 1L) s <- c(s, l - win_w + 1L)   # cover the 3' end
  s
})
idx <- rep(seq_along(seqs), lengths(starts))
win <- data.frame(seq_id = names(seqs)[idx], win_start = unlist(starts),
                  stringsAsFactors = FALSE)
message(sprintf("Tiled %d sequences into %d windows (%d bp, step %d)",
                length(seqs), nrow(win), win_w, step))

# --- score (cached: this is the expensive part) -----------------------------
# The cache is a DECLARED output, so --forcerun deletes it and rescans; that is
# intended here -- the scan is this rule's whole job.
win_cache <- out("cache/nopqs_deepg4_windows.rds")
win <- cache_or_build(win_cache, {
  wseq  <- Biostrings::subseq(seqs[idx], start = win$win_start, width = win_w)
  nfrac <- as.vector(Biostrings::letterFrequency(wseq, "N", as.prob = TRUE))
  keep  <- nfrac <= 0.1
  message(sprintf("Dropping %d windows with N fraction > 0.1", sum(!keep)))
  win$prob <- NA_real_
  ki <- which(keep)
  chunks <- split(ki, ceiling(seq_along(ki) / p$deepg4_batch))
  t0 <- Sys.time()
  for (k in seq_along(chunks)) {
    ci <- chunks[[k]]
    pr <- as.numeric(DeepG4::DeepG4(wseq[ci]))
    if (length(pr) != length(ci))
      stop(sprintf("DeepG4 returned %d predictions for %d windows (chunk %d); ",
                   length(pr), length(ci), k),
           "the N filter is out of sync with DeepG4 internals.")
    win$prob[ci] <- pr
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  chunk %d/%d scored (%d windows, %.1f min elapsed)",
                    k, length(chunks), length(ci), el))
  }
  win
})

# --- per peak ---------------------------------------------------------------
per_peak <- win %>% filter(!is.na(prob)) %>% group_by(seq_id) %>%
  summarise(n_windows = n(), max_prob = max(prob),
            best_offset = win_start[which.max(prob)] + win_w %/% 2L,
            !!!setNames(lapply(cuts, function(ct) rlang::expr(sum(prob >= !!ct))),
                        sprintf("n_win_ge_%s", cuts)),
            .groups = "drop")
per_peak <- meta %>% select(seq_id, set, parent_id, width, gc, max_z) %>%
  left_join(per_peak, by = "seq_id") %>%
  mutate(n_windows = dplyr::coalesce(n_windows, 0L))

# the OLD rule-19 score (one window: best PQS, or peak midpoint for no-PQS)
old <- readr::read_csv("results/tables/deepg4_scores.csv", show_col_types = FALSE) %>%
  filter(set == "peak") %>% select(parent_id = peak_id, rule19_prob = deepg4_prob)
per_peak <- per_peak %>%
  left_join(old, by = "parent_id") %>%
  mutate(rule19_prob = ifelse(set %in% c("noPQS", "PQS_sample"), rule19_prob, NA_real_))
readr::write_csv(per_peak, out("results/tables/nopqs_deepg4_scan.csv"))

# --- summary ----------------------------------------------------------------
summ <- per_peak %>% group_by(set) %>% summarise(
  n_peaks = n(), n_scored = sum(n_windows > 0),
  windows_total = sum(n_windows),
  max_prob_median = median(max_prob, na.rm = TRUE),
  max_prob_q90 = quantile(max_prob, .9, na.rm = TRUE),
  rule19_prob_median = median(rule19_prob, na.rm = TRUE),
  !!!setNames(lapply(cuts, function(ct) rlang::expr(mean(max_prob >= !!ct, na.rm = TRUE))),
              sprintf("frac_peaks_max_ge_%s", cuts)),
  .groups = "drop")
readr::write_csv(summ, out("results/tables/nopqs_deepg4_summary.csv"))
message("DeepG4 full-peak scan summary:"); print(as.data.frame(summ))

# --- figure: ECDF of per-peak max prob --------------------------------------
set_cols <- c(noPQS = "#B2182B", noPQS_shuffled = "grey60",
              PQS_sample = "#2166AC", random_gc = "grey25")
set_lty  <- c(noPQS = "solid", noPQS_shuffled = "22", PQS_sample = "solid", random_gc = "22")
pp <- per_peak %>% filter(!is.na(max_prob)) %>%
  mutate(set = factor(set, levels = set_levels))
lab <- summ %>% mutate(set = factor(set, levels = set_levels),
                       lab = sprintf("%s: n = %d, median max %.2f, %.0f%% >= 0.5",
                                     set, n_scored, max_prob_median,
                                     100 * .data[["frac_peaks_max_ge_0.5"]]))
p_ecdf <- ggplot(pp, aes(max_prob, colour = set, linetype = set)) +
  stat_ecdf(linewidth = 0.8) +
  geom_vline(xintercept = 0.5, colour = "grey40", linetype = "dotted") +
  scale_colour_manual(values = set_cols, labels = setNames(lab$lab, as.character(lab$set))) +
  scale_linetype_manual(values = set_lty, labels = setNames(lab$lab, as.character(lab$set))) +
  labs(x = "Best DeepG4 probability anywhere in the peak (max over 201 bp windows)",
       y = "Cumulative fraction of peaks", colour = NULL, linetype = NULL,
       title = "Full-peak DeepG4 scan: do no-PQS peaks hide G4s that scored low?",
       subtitle = sprintf(paste0("Every set tiled identically (201 bp, step %d bp). ",
                                 "Dotted line = DeepG4Scan's default call threshold 0.5"),
                          step)) +
  theme_pub() + theme(legend.position = "bottom", legend.direction = "vertical")
save_plot(p_ecdf, "nopqs_deepg4_maxprob", "26_nopqs_half_g4", fig_root, width = 8, height = 6.5)
message("Done.")
