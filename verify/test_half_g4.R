#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# The half-G4 scanner used by rule 56 (scripts/_shared/half_g4.R).
#
# Needs no bigWigs, no cache/ and no genome: every input is a synthetic string
# whose correct answer can be written down by hand. What these checks rule out
# are the failure modes that would be invisible in a per-kb rate:
#   - a 3-run site split into a 2-run site plus a leftover run (the tier would
#     depend on scan order, which is why the scanner chains runs instead of
#     using a lazy regex),
#   - a full 4-run G4 counted as a half site,
#   - a minus-strand (C-run) site missed or given the wrong sign,
#   - a GG site double-counted on top of the G3 site at the same locus,
#   - overlapping sites within a strand, which would inflate every rate.
# ---------------------------------------------------------------------------
source("scripts/_shared/helpers.R")
source("scripts/_shared/half_g4.R")
ok <- TRUE
chk <- function(lab, cond) { cat(sprintf("%-62s %s\n", lab, if (isTRUE(cond)) "PASS" else "FAIL")); if (!isTRUE(cond)) ok <<- FALSE }

A <- function(n) strrep("A", n)
scan1 <- function(s, ...) half_g4_scan(c(x = s), ...)

# 1. the basic T2, both strands, and coordinates
r <- scan1(paste0(A(5), "GGGAGGG", A(5)))
chk("GGGAGGG -> exactly one site", nrow(r) == 1)
chk("GGGAGGG -> tier T2 on +", identical(r$tier, "T2") && r$strand == "+")
chk("GGGAGGG -> start 6, end 12", r$start == 6 && r$end == 12)
chk("GGGAGGG -> run_lens 3,3 / loops 1 / tetrads 3",
    r$run_lens == "3,3" && r$loops == "1" && r$tetrads == 3)
chk("GGGAGGG -> half_score 18/7 = 2.571", abs(r$half_score - 18/7) < 1e-9)

r <- scan1(paste0(A(5), "CCCTCCC", A(5)))
chk("CCCTCCC -> T2 on the minus strand", nrow(r) == 1 && r$tier == "T2" && r$strand == "-")
chk("CCCTCCC -> half_score positive (oriented to its own strand)",
    abs(r$half_score - 18/7) < 1e-9)

# 2. chain tiers
r <- scan1(paste0(A(5), "GGGAGGGAGGG", A(5)))
chk("GGGAGGGAGGG -> ONE T3 (not T2 + leftover)", nrow(r) == 1 && r$tier == "T3")
r <- scan1(paste0(A(5), "GGGAGGGAGGGAGGG", A(5)))
chk("4 G3 runs -> ONE T4, no T3/T2", nrow(r) == 1 && r$tier == "T4")
r <- scan1(paste0(A(5), "GGGAGGGAGGGAGGGAGGG", A(5)))
chk("5 G3 runs -> still ONE T4 with n_runs 5", nrow(r) == 1 && r$tier == "T4" && r$n_runs == 5)

# 3. loop limits
r <- scan1(paste0(A(5), "GGG", A(7), "GGG", A(5)))
chk("loop 7 at loop_max 7 -> T2", nrow(r) == 1 && r$tier == "T2")
chk("loop 7 -> half_score 18/13", abs(r$half_score - 18/13) < 1e-9)
r <- scan1(paste0(A(5), "GGG", A(8), "GGG", A(5)))
chk("loop 8 at loop_max 7 -> nothing", nrow(r) == 0)
r <- scan1(paste0(A(5), "GGG", A(8), "GGG", A(5)), loop_max = 12L)
chk("loop 8 at loop_max 12 -> T2", nrow(r) == 1 && r$tier == "T2")
r <- scan1(paste0(A(5), "GGGAGGG", A(20), "GGGAGGG", A(5)))
chk("two T2 separated by a 20 nt gap -> two T2, not a T4", nrow(r) == 2 && all(r$tier == "T2"))

# 4. maximal runs: a homopolymer is one run, not a half site
r <- scan1(paste0(A(5), "GGGGGGG", A(5)))
chk("GGGGGGG alone -> nothing (one run, not two)", nrow(r) == 0)

# 5. GG sensitivity tier
r <- scan1(paste0(A(5), "GGAGG", A(5)))
chk("GGAGG -> T2gg only", nrow(r) == 1 && r$tier == "T2gg")
r <- scan1(paste0(A(5), "GGAGG", A(5)), include_gg = FALSE)
chk("GGAGG with include_gg = FALSE -> nothing", nrow(r) == 0)
r <- scan1(paste0(A(5), "GGGAGGG", A(5)))
chk("GGGAGGG -> no extra T2gg at the same locus", !any(r$tier == "T2gg"))
r <- scan1(paste0(A(5), "GGAGGAGGAGG", A(5)))
chk("4 GG runs -> T4gg (2-tetrad full G4)", nrow(r) == 1 && r$tier == "T4gg")
r <- scan1(paste0(A(5), "GGGAGG", A(5)))
chk("GGGAGG (one G3 run) -> T2gg with tetrads 2", nrow(r) == 1 && r$tier == "T2gg" && r$tetrads == 2)

# 6. both strands in one sequence, and an opposite-strand run inside a loop
r <- scan1(paste0("GGGAGGG", A(10), "CCCACCC"))
chk("G site and C site in one sequence -> one + and one -",
    nrow(r) == 2 && setequal(r$strand, c("+", "-")))
r <- scan1(paste0(A(5), "GGGACCCAGGG", A(5)), loop_max = 7L)
chk("C-run inside the loop lowers the score (G4Hunter)",
    nrow(r[r$strand == "+", ]) == 1 &&
    abs(r$half_score[r$strand == "+"] - (18 - 9) / 11) < 1e-9)

# 7. no overlapping sites within a strand, on a realistic random sequence
set.seed(1)
rnd <- paste(sample(c("A", "C", "G", "T"), 20000, TRUE, prob = c(.2, .3, .3, .2)), collapse = "")
r <- scan1(rnd)
no_ovl <- all(vapply(split(r, list(r$strand, r$tier %in% c("T4gg", "T2gg"))), function(d) {
  d <- d[order(d$start), ]
  nrow(d) < 2 || all(d$start[-1] > d$end[-nrow(d)])
}, logical(1)))
chk("random 20 kb: no two sites overlap within a strand/run class", no_ovl)
chk("random 20 kb: every site has >= 2 runs", all(r$n_runs >= 2))
# half_score is NOT bounded below by 0: a GG chain whose loop holds a longer
# C-run (GGCCCGG, GGTGACCCCGG) is net C-rich and scores <= 0 -- measured here,
# 17 of 623 sites, all GG-tier. That is G4Hunter correctly calling it a poor
# site, which is what the null-calibrated high-confidence threshold is for.
chk("random 20 kb: every half_score in [-4, 4]", all(abs(r$half_score) <= 4))
g3 <- r[r$tier %in% c("T4", "T3", "T2"), ]
chk("random 20 kb: G3-tier half_score all > 0", all(g3$half_score > 0))
chk("random 20 kb: G3 tiers have tetrads >= 3", all(g3$tetrads >= 3))

# 8. vector input + names, and the empty case
r <- half_g4_scan(c(a = "GGGAGGG", b = "AAAAAAA", c = "CCCTCCC"))
chk("vector input: seq_id carried through, empty element skipped",
    identical(sort(unique(r$seq_id)), c("a", "c")))
r <- half_g4_scan(c(a = "AAAA"))
chk("nothing found -> zero-row frame with the full schema",
    nrow(r) == 0 && all(c("seq_id", "tier", "half_score") %in% names(r)))

cat("\n", if (ok) "ALL HALF-G4 SCANNER TESTS PASSED" else "SOME TESTS FAILED", "\n")
quit(status = if (ok) 0 else 1)
