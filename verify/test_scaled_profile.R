#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Geometry of the per-region profile grids in helpers.R: the scaled-region one
# (scaled_bin_ranges / _profile_matrix / _positions) and the COMPOSITE one that
# replaced it in rules 52/53 (composite_bin_ranges / _profile_matrix /
# _positions / _axis), which keeps both anchor blocks linear in bp.
#
# Needs no bigWigs and no cache/: coverage is a synthetic Rle whose VALUE equals
# its POSITION, so a bin's mean must come back as that bin's midpoint. That
# turns "the binning is correct" into an exact arithmetic identity rather than
# an eyeball check on a figure, and it is what rules out the two failure modes
# that would be invisible downstream:
#   - bins that gap or overlap (signal counted twice or dropped), and
#   - a minus-strand row that is not reversed, which would average TSS signal
#     into the TES end of the curve and quietly flatten both.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(GenomicRanges); library(IRanges); library(S4Vectors)})
source("scripts/_shared/helpers.R")
ok <- TRUE
chk <- function(lab, cond) { cat(sprintf("%-58s %s\n", lab, if (isTRUE(cond)) "PASS" else "FAIL")); if(!isTRUE(cond)) ok <<- FALSE }

FL <- 1000L; NF <- 20L; NB <- 60L; nb <- 2L*NF+NB
g <- GRanges(c("chr1","chr1","chr1"),
             IRanges(c(10001, 50001, 200001), c(15000, 50600, 1200000)),
             strand = c("+","-","+"))
names(g) <- c("plus5k","minus600","plus1Mb")
r <- scaled_bin_ranges(g, FL, NF, NB)
chk("n_bins == 100", r$n_bins == nb)
st <- matrix(r$start, nrow=3); en <- matrix(r$end, nrow=3)

# 1. bins tile the full locus with no gap and no overlap
for (i in 1:3) {
  lo <- start(g)[i]-FL; hi <- end(g)[i]+FL
  chk(sprintf("gene %d: first bin starts at TSS-flank", i), st[i,1]==lo)
  chk(sprintf("gene %d: last bin ends at TES+flank", i),    en[i,nb]==hi)
  chk(sprintf("gene %d: contiguous (no gap/overlap)", i),   all(st[i,-1]==en[i,-nb]+1))
  chk(sprintf("gene %d: all bins >= 1 bp", i),              all(en[i,]>=st[i,]))
}
# 2. body bins sit exactly on the gene, flanks exactly off it
for (i in 1:3) {
  chk(sprintf("gene %d: body bin 1 starts at gene start", i), st[i,NF+1]==start(g)[i])
  chk(sprintf("gene %d: body bin %d ends at gene end", i, NB), en[i,NF+NB]==end(g)[i])
  chk(sprintf("gene %d: upstream flank is exactly %d bp", i, FL),
      en[i,NF]-st[i,1]+1 == FL)
  chk(sprintf("gene %d: downstream flank is exactly %d bp", i, FL),
      en[i,nb]-st[i,NF+NB+1]+1 == FL)
}
# 3. strand flip: synthetic coverage that ramps left->right must come back
#    ramping TSS->TES for + and reversed for -
cov <- S4Vectors::Rle(seq_len(1201000))                 # value == position
covL <- IRanges::RleList(chr1 = cov)
assign(".fake", covL, envir = globalenv())
# monkeypatch read_bigwig_rle for the test
read_bigwig_rle <- function(...) .fake
m <- scaled_profile_matrix(g, "fake.bw", FL, NF, NB)
chk("plus-strand gene: signal increases TSS -> TES", all(diff(m[1,]) > 0))
chk("minus-strand gene: signal DEcreases in genomic order after flip",
    all(diff(m[2,]) < 0))
chk("plus 1Mb gene: increases", all(diff(m[3,]) > 0))
# bin means equal the midpoint of each bin (value == position)
expect1 <- (st[1,]+en[1,])/2
chk("plus gene bin means == bin midpoints", max(abs(m[1,]-expect1)) < 1e-6)
expect2 <- rev((st[2,]+en[2,])/2)
chk("minus gene bin means == REVERSED bin midpoints", max(abs(m[2,]-expect2)) < 1e-6)
# 4. positions axis
p <- scaled_positions(FL,NF,NB)
chk("positions length == n_bins", length(p)==nb)
chk("upstream all < 0", all(p[1:NF] < 0))
chk("body in (0,1)", all(p[(NF+1):(NF+NB)] > 0 & p[(NF+1):(NF+NB)] < 1))
chk("downstream all > 1", all(p[(NF+NB+1):nb] > 1))
chk("positions strictly increasing", all(diff(p) > 0))
# 5. short-gene guard
short <- GRanges("chr1", IRanges(1000, 1030), strand="+")
chk("too-short gene raises an error",
    inherits(try(scaled_bin_ranges(short, FL, NF, NB), silent=TRUE), "try-error"))
# 6. out-of-coverage bins: NA, never clamped. A bigWig's Rle can end before the
#    BSgenome chromosome does (measured: GSM8305989 chr16 stops 199 kb short), so
#    bins past the end must come back NA -- clamping them would shrink that gene's
#    bins and rescale its body against every other gene's.
# The expectation is DERIVED from the bin ranges rather than hard-coded, so the
# test states the rule ("exactly the out-of-bounds bins are NA") instead of a
# count that depends on where the test gene happens to sit.
cov_len <- length(cov)
expect_na <- function(gr) {
  b <- scaled_bin_ranges(gr, FL, NF, NB)
  oob <- b$start < 1 | b$end > cov_len
  if (as.character(strand(gr)) == "-") rev(oob) else oob
}
edge <- GRanges("chr1", IRanges(100, 5000), strand="+")           # TSS-1000 < 1
me <- scaled_profile_matrix(edge, "fake.bw", FL, NF, NB)
chk("upstream: exactly the out-of-bounds bins are NA",
    identical(is.na(me[1, ]), expect_na(edge)))
chk("upstream: at least one bin was off the chromosome", any(expect_na(edge)))
chk("upstream: every in-bounds bin still computed",
    all(is.finite(me[1, !expect_na(edge)])))
tail_g <- GRanges("chr1", IRanges(1198000, 1200500), strand="+")  # TES+1000 past end
mt <- scaled_profile_matrix(tail_g, "fake.bw", FL, NF, NB)
chk("downstream: exactly the out-of-bounds bins are NA",
    identical(is.na(mt[1, ]), expect_na(tail_g)))
chk("downstream: at least one bin was past coverage", any(expect_na(tail_g)))
chk("downstream: a partially-covered gene keeps its covered bins",
    all(is.finite(mt[1, !expect_na(tail_g)])))
minus_edge <- GRanges("chr1", IRanges(100, 5000), strand="-")     # flip applies to NAs too
mm <- scaled_profile_matrix(minus_edge, "fake.bw", FL, NF, NB)
chk("minus strand: NA pattern is reversed with the row",
    identical(is.na(mm[1, ]), expect_na(minus_edge)))

# ===========================================================================
# 7. COMPOSITE layout (composite_bin_ranges / _positions / _profile_matrix).
#
# This is the layout rules 52/53 actually use. Its whole reason for existing is
# that the two ANCHOR blocks are linear in bp, so the assertions below are about
# exactly that: every anchor bin is exactly anchor_bin wide, the TSS sits on a
# block boundary, and the interior tiles the remainder with no gap or overlap.
# If an anchor bin were ever scaled the figure would silently go back to
# reporting a 500 bp signal decay as one number.
# ===========================================================================
cat("\n--- composite layout ---\n")
W <- 2000L; AB <- 50L; CB <- 60L
NA_ <- as.integer(2*W/AB); cnb <- 2L*NA_ + CB
cg <- GRanges(c("chr1","chr1"), IRanges(c(100001, 300001), c(120000, 305000)),
              strand = c("+","-"))
names(cg) <- c("plus20k","minus5k")
cr <- composite_bin_ranges(cg, W, AB, CB, min_interior = 1000L)
chk("composite: n_bins == 2*80 + 60 == 220", cr$n_bins == cnb && cnb == 220L)
cst <- matrix(cr$start, nrow=2); cen <- matrix(cr$end, nrow=2)
for (i in 1:2) {
  s <- start(cg)[i]; e <- end(cg)[i]
  chk(sprintf("composite gene %d: first bin starts at left-1 anchor", i), cst[i,1]==s-W)
  chk(sprintf("composite gene %d: last bin ends at right+1 anchor", i),  cen[i,cnb]==e+W)
  chk(sprintf("composite gene %d: contiguous (no gap/overlap)", i),
      all(cst[i,-1]==cen[i,-cnb]+1))
  # THE point of the layout: anchor bins are fixed bp, not scaled
  chk(sprintf("composite gene %d: all LEFT anchor bins exactly %d bp", i, AB),
      all(cen[i,1:NA_]-cst[i,1:NA_]+1 == AB))
  chk(sprintf("composite gene %d: all RIGHT anchor bins exactly %d bp", i, AB),
      all(cen[i,(NA_+CB+1):cnb]-cst[i,(NA_+CB+1):cnb]+1 == AB))
  chk(sprintf("composite gene %d: left anchor block spans exactly 2W", i),
      cen[i,NA_]-cst[i,1]+1 == 2*W)
  chk(sprintf("composite gene %d: gene start is a left-block bin boundary", i),
      s %in% cst[i,1:NA_])
  chk(sprintf("composite gene %d: gene end is a right-block bin end", i),
      e %in% cen[i,(NA_+CB+1):cnb])
  chk(sprintf("composite gene %d: interior tiles [start+W, end-W]", i),
      cst[i,NA_+1]==s+W && cen[i,NA_+CB]==e-W)
}
# value == position again: anchor bins must come back as exact midpoints
ccov <- RleList(chr1 = Rle(seq_len(400000)), compress = FALSE)
assign(".fake", ccov, envir = globalenv())
cm <- composite_profile_matrix(cg, "fake.bw", W, AB, CB, min_interior = 1000L)
chk("composite plus gene: bin means == bin midpoints",
    max(abs(cm[1,] - (cst[1,]+cen[1,])/2)) < 1e-6)
chk("composite minus gene: bin means == REVERSED bin midpoints",
    max(abs(cm[2,] - rev((cst[2,]+cen[2,])/2))) < 1e-6)
# a minus gene's TSS block must end up on the LEFT of the row after the flip
chk("composite minus gene: row is reversed (signal decreases)", all(diff(cm[2,]) < 0))
chk("composite plus gene: row is not reversed (signal increases)", all(diff(cm[1,]) > 0))

# positions / axis
cp <- composite_positions(W, AB, CB)
chk("composite positions: nrow == n_bins", nrow(cp) == cnb)
chk("composite positions: x strictly increasing", all(diff(cp$x) > 0))
chk("composite positions: block sizes 80/60/80",
    identical(as.integer(table(cp$block)), c(NA_, CB, NA_)))
chk("composite positions: bp is NA across the interior and only there",
    identical(is.na(cp$bp), cp$block == "interior"))
chk("composite positions: TSS block is real bp, first centre at -W+AB/2",
    abs(cp$x[1] - (-W + AB/2)) < 1e-9)
chk("composite positions: TES block centred on x == 4W",
    abs(mean(range(cp$x[cp$block=="tes"])) - 4*W) < 1e-9)
ax <- composite_axis(W)
chk("composite axis: TSS break at 0 and TES break at 4W",
    identical(ax$anchors, c(0, 4*as.numeric(W))))
chk("composite axis: 7 breaks, 7 labels", length(ax$breaks)==7 && length(ax$labels)==7)

# the width guard -- a gene too short for two non-overlapping anchor windows
tooshort <- GRanges("chr1", IRanges(200000, 204000), strand="+")   # 4001 bp < 2W+1000
chk("composite: gene too short for both anchor windows raises an error",
    inherits(try(composite_bin_ranges(tooshort, W, AB, CB, 1000L), silent=TRUE),
             "try-error"))
cat("\n", if (ok) "ALL GEOMETRY TESTS PASSED" else "SOME TESTS FAILED", "\n")
quit(status = if (ok) 0 else 1)
