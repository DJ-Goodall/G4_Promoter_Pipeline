#!/usr/bin/env Rscript
# ============================================================================
# 03_deepg4_predict.R   (env: deepg4)
#
# Score G4-formation probability with the DeepG4 CNN for peak 201 bp windows and
# the GC-matched background. This is the ONLY step that needs Python/TensorFlow;
# isolating it means a DeepG4 setup failure does not block the topology branch.
#
# Inputs:   cache/peak_seqs_201bp.fa, cache/background_seqs.fa
# Output:   results/tables/deepg4_scores.csv  (peak_id, deepg4_prob, set)
#
# Requires the DeepG4 R package installed into the active conda env:
#   conda run -n deepg4 R -e 'remotes::install_github("raphaelmourad/DeepG4")'
# ============================================================================

Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", CUDA_VISIBLE_DEVICES = "-1")  # quiet, CPU-only

suppressPackageStartupMessages({
  library(Biostrings); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables")

# --- Point reticulate at the TensorFlow conda env --------------------------
# This script runs in SYSTEM R (which holds the DeepG4 package); only TF/keras
# come from a conda env (default "tf-new"), reached via reticulate.
suppressPackageStartupMessages(library(reticulate))
tf_env <- cfg$deepg4$reticulate_condaenv %||% "tf-new"
ok_env <- tryCatch({ reticulate::use_condaenv(tf_env, required = TRUE); TRUE },
                   error = function(e) { message("use_condaenv('", tf_env, "') failed: ",
                                                 conditionMessage(e)); FALSE })

ok_deepg4 <- requireNamespace("DeepG4", quietly = TRUE)
if (!ok_deepg4) {
  stop("DeepG4 R package not found in your R library. Install with:\n",
       "  install.packages('remotes'); remotes::install_github('raphaelmourad/DeepG4')\n",
       "and ensure conda env '", tf_env, "' has tensorflow/keras importable.")
}
suppressPackageStartupMessages(library(DeepG4))

run_deepg4 <- function(fasta_path, set_label) {
  seqs <- Biostrings::readDNAStringSet(fasta_path)
  # DeepG4 requires fixed 201 bp; keep only conforming sequences.
  w <- BiocGenerics::width(seqs)
  if (any(w != cfg$sequences$deepg4_width)) {
    keep <- w == cfg$sequences$deepg4_width
    message("  ", set_label, ": dropping ", sum(!keep),
            " sequences not exactly ", cfg$sequences$deepg4_width, " bp")
    seqs <- seqs[keep]
  }
  # CRITICAL: DeepG4::DeepG4() silently drops sequences whose N frequency > 0.1
  # *before* predicting (see DeepG4.R: `X <- X[!testNFreq]`), so its returned
  # matrix is shorter than the input and no longer aligned to names(seqs). We
  # replicate that exact filter here so peak_id stays row-aligned with the
  # predictions (and the dropped peaks are simply absent downstream).
  n_freq <- as.vector(Biostrings::letterFrequency(seqs, "N", as.prob = TRUE))
  hi_n <- n_freq > 0.1
  if (any(hi_n)) {
    message("  ", set_label, ": dropping ", sum(hi_n),
            " sequences with N frequency > 0.1 (DeepG4 cannot score these)")
    seqs <- seqs[!hi_n]
  }
  if (length(seqs) < 1) {
    warning(set_label, ": no sequences left to score after N filtering")
    return(data.frame(peak_id = character(0), deepg4_prob = numeric(0),
                      set = character(0), stringsAsFactors = FALSE))
  }
  message("Running DeepG4 on ", length(seqs), " ", set_label, " sequences...")
  pred <- as.numeric(DeepG4::DeepG4(seqs))   # N x 1 probability matrix -> vector
  if (length(pred) != length(seqs)) {
    stop(sprintf("DeepG4 returned %d predictions for %d %s sequences; ",
                 length(pred), length(seqs), set_label),
         "N-content filter is out of sync with DeepG4 internals.")
  }
  data.frame(peak_id = names(seqs),
             deepg4_prob = pred,
             set = set_label, stringsAsFactors = FALSE)
}

scores <- rbind(
  run_deepg4("cache/peak_seqs_201bp.fa", "peak"),
  run_deepg4("cache/background_seqs.fa", "background")
)
readr::write_csv(scores, "results/tables/deepg4_scores.csv")

message(sprintf("Done. DeepG4 scored %d peaks + %d background (median prob peak=%.3f, bg=%.3f)",
                sum(scores$set == "peak"), sum(scores$set == "background"),
                median(scores$deepg4_prob[scores$set == "peak"], na.rm = TRUE),
                median(scores$deepg4_prob[scores$set == "background"], na.rm = TRUE)))
