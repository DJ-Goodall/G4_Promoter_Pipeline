# ============================================================================
# scripts/_shared/half_g4.R
#
# The half-G4 site scanner for Stage O (rules 54-56). Kept OUT of helpers.R on
# purpose: helpers.R is the contract six committed tables depend on, and this
# is a new algorithm with no users yet. Sourced by scripts/56_nopqs_half_g4.R
# and verify/test_half_g4.R. Needs g4hunter_one() from helpers.R, so source
# helpers.R first.
#
# ---------------------------------------------------------------------------
# WHAT A "HALF SITE" IS HERE
#
# There is no published half-G4 caller, so this is our own definition, written
# to be transparent rather than clever:
#
#   1. Find MAXIMAL G-runs of length >= min_run on the + strand (G{min_run,})
#      and maximal C-runs on the + strand, which are G-runs of the - strand.
#   2. Link consecutive runs on the same strand into a CHAIN while the loop
#      between them is within [loop_min, loop_max].
#   3. Tier the chain by its run count:
#        T4  >= 4 runs  a regex full G4. Sanity check: pqsfinder should
#                       already have caught these, so ~0 in no-PQS peaks.
#        T3     3 runs  the 3+1 near-miss -- one tract short of intramolecular;
#                       a 4th tract from another strand or an RNA completes it.
#        T2     2 runs  the DNA:RNA hybrid G4 substrate: two DNA G-tracts are
#                       enough when the transcript supplies the other two
#                       (Zheng et al. 2013 NAR); also the bimolecular G-hairpin.
#   4. A sensitivity tier at gg_min_run (= 2, i.e. GG runs): chains of 2-3 GG+
#      runs that do NOT overlap any min_run chain on the same strand -> "T2gg";
#      GG chains of >= 4 runs -> "T4gg" (2-tetrad full G4, sanity check).
#      GG-based half sites are common in any sequence, so this tier is reported
#      but never the headline.
#
# Chains cannot overlap within a strand by construction (each run belongs to
# exactly one chain), so no collapsing step is needed -- unlike pqsfinder's
# overlapping registers, which rule 41 has to collapse.
#
# Rejected alternative: regex matching (e.g. "G{3,}.{1,7}?G{3,}") as the site
# definition. A lazy regex walks left to right and pairs runs greedily, so
# GGGAGGGAGGG is reported as ONE 2-run hit plus a leftover run instead of one
# 3-run site, and the tier then depends on the regex engine rather than the
# sequence. The chain formulation is order-free.
#
# Maximal runs mean a single long run (GGGGGGG) is ONE run, not two runs with a
# G "loop". pqsfinder would split it; for a half-site count that would make a
# homopolymer look like a two-tract substrate, which it is not.
#
# ---------------------------------------------------------------------------
# half_score
#
# G4Hunter (Bedrat et al. 2016) over the chain span, oriented so a site on its
# own G-rich strand scores positive: g4hunter_one(span) on +, -g4hunter_one on -.
# Range 0-4. It rewards long runs and penalises long loops (loop bases score 0,
# and any opposite-strand run inside a loop scores negative) -- the same logic
# pqsfinder's loop penalty encodes, without inventing new weights.
#   GGGAGGG  -> (3*3 + 0 + 3*3) / 7 = 2.571
#   GGGAAAAAAAGGG -> 18 / 13 = 1.385
# ============================================================================

# Maximal runs of `letter` with length >= min_run in one string.
# Returns a data.frame(start, end), 1-based, possibly zero rows.
.hg4_runs <- function(s, letter, min_run) {
  m <- gregexpr(sprintf("%s{%d,}", letter, min_run), s, perl = TRUE)[[1]]
  if (m[1] == -1L) return(data.frame(start = integer(0), end = integer(0)))
  st <- as.integer(m)
  data.frame(start = st, end = st + attr(m, "match.length") - 1L)
}

# Chain consecutive runs whose loop is within [loop_min, loop_max].
# Returns a list of integer vectors (row indices into `runs`), one per chain,
# including singleton "chains" (length 1) which callers drop.
.hg4_chains <- function(runs, loop_min, loop_max) {
  n <- nrow(runs)
  if (n == 0L) return(list())
  if (n == 1L) return(list(1L))
  loops <- runs$start[-1] - runs$end[-n] - 1L
  link  <- loops >= loop_min & loops <= loop_max
  # a new chain starts at run 1 and after every unlinked gap
  chain_id <- cumsum(c(TRUE, !link))
  split(seq_len(n), chain_id)
}

# All chains of >= 2 runs for one letter/min_run on one string, as a data.frame.
.hg4_chain_table <- function(s, letter, min_run, loop_min, loop_max) {
  runs <- .hg4_runs(s, letter, min_run)
  ch   <- .hg4_chains(runs, loop_min, loop_max)
  ch   <- ch[lengths(ch) >= 2L]
  if (!length(ch)) return(NULL)
  do.call(rbind, lapply(ch, function(ix) {
    r  <- runs[ix, , drop = FALSE]
    rl <- r$end - r$start + 1L
    lp <- r$start[-1] - r$end[-nrow(r)] - 1L
    data.frame(start = r$start[1], end = r$end[nrow(r)],
               n_runs = nrow(r), tetrads = min(rl),
               run_lens = paste(rl, collapse = ","),
               loops = paste(lp, collapse = ","),
               stringsAsFactors = FALSE)
  }))
}

# Scan ONE sequence (character) on both strands. Returns a data.frame with one
# row per site, or NULL. Coordinates are 1-based on the + strand of `s`.
half_g4_scan_one <- function(s, min_run = 3L, gg_min_run = 2L,
                             loop_min = 1L, loop_max = 7L, include_gg = TRUE) {
  s <- toupper(as.character(s))
  if (is.na(s) || !nzchar(s)) return(NULL)
  out <- list()
  for (strd in c("+", "-")) {
    letter <- if (strd == "+") "G" else "C"
    g3 <- .hg4_chain_table(s, letter, min_run, loop_min, loop_max)
    if (!is.null(g3)) {
      g3$tier <- ifelse(g3$n_runs >= 4L, "T4", ifelse(g3$n_runs == 3L, "T3", "T2"))
      g3$strand <- strd
      out[[length(out) + 1L]] <- g3
    }
    if (isTRUE(include_gg)) {
      gg <- .hg4_chain_table(s, letter, gg_min_run, loop_min, loop_max)
      if (!is.null(gg)) {
        # a GG chain that overlaps any min_run chain on this strand is the same
        # locus seen at a laxer run length -- not an extra site
        if (!is.null(g3)) {
          hit <- vapply(seq_len(nrow(gg)), function(i)
            any(gg$start[i] <= g3$end & gg$end[i] >= g3$start), logical(1))
          gg <- gg[!hit, , drop = FALSE]
        }
        if (nrow(gg)) {
          gg$tier <- ifelse(gg$n_runs >= 4L, "T4gg", "T2gg")
          gg$strand <- strd
          out[[length(out) + 1L]] <- gg
        }
      }
    }
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  span <- substring(s, res$start, res$end)
  sc   <- vapply(span, g4hunter_one, numeric(1), USE.NAMES = FALSE)
  res$half_score <- ifelse(res$strand == "+", sc, -sc)
  res
}

# Scan a named character vector / DNAStringSet. Returns one data.frame with a
# seq_id column (the element name) -- zero rows if nothing was found.
half_g4_scan <- function(seqs, min_run = 3L, gg_min_run = 2L,
                         loop_min = 1L, loop_max = 7L, include_gg = TRUE) {
  if (methods::is(seqs, "XStringSet")) seqs <- stats::setNames(as.character(seqs), names(seqs))
  ids <- names(seqs)
  if (is.null(ids)) ids <- as.character(seq_along(seqs))
  res <- lapply(seq_along(seqs), function(i) {
    r <- half_g4_scan_one(seqs[[i]], min_run, gg_min_run, loop_min, loop_max, include_gg)
    if (!is.null(r)) r$seq_id <- ids[i]
    r
  })
  res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(res))
    return(data.frame(start = integer(0), end = integer(0), n_runs = integer(0),
                      tetrads = integer(0), run_lens = character(0),
                      loops = character(0), tier = character(0),
                      strand = character(0), half_score = numeric(0),
                      seq_id = character(0), stringsAsFactors = FALSE))
  res <- res[, c("seq_id", "start", "end", "strand", "tier", "n_runs", "tetrads",
                 "run_lens", "loops", "half_score")]
  rownames(res) <- NULL
  res
}

# Tier vocabulary, in display order. "G3" tiers are the hypothesis; the GG
# tiers are sensitivity; T4/T4gg are sanity checks (full G4 by regex).
half_g4_tiers       <- c("T4", "T3", "T2", "T4gg", "T2gg")
half_g4_tier_labels <- c(T4   = "T4: >=4 G3 runs (full G4, sanity)",
                         T3   = "T3: 3 G3 runs (3+1 near-miss)",
                         T2   = "T2: 2 G3 runs (hybrid substrate)",
                         T4gg = "T4gg: >=4 GG runs (2-tetrad, sanity)",
                         T2gg = "T2gg: 2-3 GG runs (sensitivity)")
