# v1.1 — response to the 2026-08-07 audit

Internal review. Not committed to the public repo.

Every number below is re-derived from the committed v1.0 tables by
`verify/rederive_audit_numbers.R`, which `stopifnot()`s each claim. Two further
scripts measure things the audit asserted without testing:
`verify/run_heuristic_benchmark.R` and `verify/calibrate_g4_threshold.R`.

**Summary.** The audit is substantially correct. One real bug, two genuine
mislabels, and one methodological problem that changes a headline conclusion.
Against that: **two of its findings understate the problem**, **one is
factually wrong**, and **two of its prescriptions do not work as written** —
one of them would have broken a table in the opposite direction.

| § | Audit finding | Verdict |
|---|---|---|
| 2.1 | Topology strata are a chimera | **Confirmed, worse than stated** — 87%/66%, not 83%/55% |
| 2.2 | Peak calling has no null | **Confirmed** — and the repo already proves it |
| 2.3 | DeepG4 AUROC 0.563 not flagged | **Confirmed** for the figure title; **overstated** for ρ=0.13 and `g4hunter_win` |
| 2.4 | `ignore.strand` bug | **Confirmed** — and the damage is wider than described |
| 2.4b | Set `g4_threshold` to 1.2 | **Prescription fails** — would return zero rows |
| 2.5 | `p_perm` is not a permutation | **Confirmed**, but the diagnosis is half the story |
| 2.6 | Spike-in dependency unstated | **Confirmed** — and the README's control claim is unsupported |
| 2.7 | z-scoring guarantees the stability result | **Wrong** — defending the code |
| 2.7b | No accessibility control | **Confirmed** |
| 2.8 | The expression endpoint is null | **Confirmed** — promote it |

---

## 2.1 Topology strata — confirmed, and worse than stated

### The share of each stratum that was never a model call

The audit says ~83% of hybrid and ~55% of antiparallel are heuristic
back-fills. Those figures mix denominators (the per-genotype composition table
against the union-level sweep). On one consistent denominator —
`topology_composition_compare.csv` "Per peak (top PQS)" against the
`precision = 0.80` row of `topology_precision_sweep.csv`:

| stratum | final n | G4SP n | heuristic n | **% heuristic** |
|---|---|---|---|---|
| parallel | 104,336 | 95,660 | 8,676 | 8.3% |
| antiparallel | 12,098 | 4,099 | 7,999 | **66.1%** |
| hybrid | 14,429 | 1,879 | 12,550 | **87.0%** |

Heuristic total = 29,225 = `n_mixed` at precision 0.80, exactly. Use
**87.0 / 66.1**.

### Push-back: κ is the weaker statistic here

κ = −0.0098 is arithmetically correct (p_o = 0.3351, p_e = 0.3415,
n = 101,638 — reproduced). But κ is depressed by the marginal mismatch: G4SP
calls 94.1% parallel, the heuristic 34.4%. With marginals that divergent κ is
a known-unreliable summary, and "the two methods agree worse than chance"
invites a methodological rebuttal that would distract from a finding which
does not need one.

The marginal-free statement is stronger and unarguable:

- of the heuristic's **41,122 hybrid calls, 741 (1.8%)** are G4SP-hybrid
- of its **25,539 antiparallel calls, 191 (0.7%)** are G4SP-antiparallel

Lead with those. v1.1 writes both to `topology_concordance_stats.csv`.

### The audit missed the evidence that settles it

The audit writes: *"That would be defensible if the heuristic were a validated
stand-in. It is not."* That is an assertion — and the repository already
contains everything needed to test it. `data/known_topology_controls.tsv` plus
`26b_control_benchmark.py` scores **G4SP at 8/10** on literature-validated
control quadruplexes spanning all three topology classes. The heuristic had
never been run through the same test.

v1.1 adds `26d_heuristic_control_benchmark.R` and does exactly that, using the
pipeline's own pqsfinder settings. Result:

| control | reported | G4SP | loop heuristic | loops |
|---|---|---|---|---|
| cMYC_Pu27 | parallel | parallel ✓ | parallel ✓ | 1/0/1 |
| cKIT1 | parallel | parallel ✓ | antiparallel ✗ | 1/4/4 |
| cKIT2 | parallel | parallel ✓ | hybrid ✗ | 1/5/1 |
| VEGF | parallel | parallel ✓ | parallel ✓ | 1/0/1 |
| KRAS | parallel | parallel ✓ | hybrid ✗ | 1/2/4 |
| HRAS | parallel | antiparallel ✗ | antiparallel ✗ | 4/4/3 |
| BCL2_Pu39 | parallel | hybrid ✗ | hybrid ✗ | 3/5/3 |
| Telomere_human | hybrid | hybrid ✓ | parallel ✗ | 3/3/3 |
| TBA | antiparallel | antiparallel ✓ | parallel ✗ | 2/3/2 |
| Oxytricha | antiparallel | antiparallel ✓ | no PQS ✗ | — |

**G4ShapePredictor 8/10. Loop heuristic 2/10** — below the 33% expected from
guessing at a three-class problem, and the two it gets right are the trivial
1/0/1 cases. It fails every non-parallel control, calling both the human
telomere (hybrid) and the thrombin-binding aptamer (antiparallel) "parallel".

This is a far stronger justification for dropping the back-fill than κ, and it
is the single highest-value thing the audit left on the table. In fairness to
the heuristic: loop length genuinely is a determinant of topology (Guédin et
al. 2010), so it is not arbitrary — but applied to the top PQS of an isolated
sequence, with no cation, crowding or chromatin context, it does not recover
known structures.

### Push-back: `precision: 0.75` is not the cleanest fix

The audit suggests 0.75 as the working default. But look at the whole sweep:

| precision | n_mixed | n_hybrid |
|---|---|---|
| **default** | **0** | 11,062 |
| 0.70 | 59 | 11,003 |
| 0.75 | 6,090 | 7,912 |
| 0.80 | 29,225 | 1,879 |
| 0.85 | 64,960 | 37 |

At `default` (argmax, no abstention) there is **no mixed bucket at all** — so
no back-fill, no chimera — while retaining 11,062 hybrids. The mixed bucket is
created *entirely* by precision thresholding. Argmax-primary with confidence
reported as a covariate removes the problem without shrinking n. Note also
that hybrid is by far G4SP's least confident class and collapses fastest
(11,062 → 1,879 → 37), so choosing a precision *is* choosing how big the
hybrid stratum is.

### What v1.1 does

- `21_assign_topology.R`: the `!is.na(topology_heuristic)` fallback is gone.
  Abstentions terminate at `"mixed"`.
- `"mixed"` added to `topology_levels` / `labels` / `palette` in `helpers.R`,
  so every downstream figure picks it up automatically.
- κ **and** the marginal-free rates persisted to
  `topology_concordance_stats.csv`. In v1.0 κ existed only in a `message()` to
  a git-ignored `logs/` and a figure subtitle under a git-ignored
  `results/figures/` — i.e. in no committed artefact at all.
- `topology_source` aggregated into a table and a figure. **The column already
  existed** in v1.0 and was even written to `peak_topology.csv`; it was simply
  never counted or plotted, so a reader had no way to see the provenance. This
  fix was nearly free.
- `26d_heuristic_control_benchmark.R` added.

**Consequence to accept:** the "hybrid > parallel > antiparallel" TSS ordering
must be re-run on G4SP-definite peaks only. On 1,879 hybrids rather than
14,429 it may not survive. That is the honest test.

---

## 2.2 Peak calling — confirmed, and the proof is already in the repo

The audit argues from first principles that a global z ≥ 2 has no null in
which a peak can fail. `peak_summary.csv` demonstrates it outright:

| set | n_peaks |
|---|---|
| G4_BG4_WT | 107,276 |
| **G4_BG4_ERCCKO** | **114,980** |
| G4_BG4_union | 147,843 |
| **G4_BG4_ercc_union** | **150,977** |

**The ERCC control tracks yield more peaks than the real samples.** A caller
that returns ~110k peaks from a control is not thresholding on enrichment. The
audit didn't cite this; it is the sharpest available evidence.

Confirmed mechanics: `support = rowSums(sig_mat > 0)` is non-zero *coverage*
in ≥ 2 replicates, not an independent per-replicate test, and the z-test is
applied once to the replicate mean. Two things the audit missed:
`min_width: 500 == bin_size: 500` makes the width filter inert, and
`min.gapwidth = 500` with 500 bp bins merges only *directly adjacent* bins
(one skipped bin gives a gap of exactly 500, which is not `< 500`).

**Mild push-back on "genuine reproduction".** Read in context (README:64–66)
the sentence contrasts *recomputing* peaks with *downloading* a deposited peak
file. It does not claim to reproduce Sato et al.'s MACS3 calls. It is genuinely
ambiguous and is reworded in `README-deltas.md`, but it is not the overclaim
the audit reads it as.

**Also worth stating:** the audit's caution that cross-genotype peak counts
carry no KO information is correct, but the pipeline never reads them that
way. `23_gain_loss.R` already does limma-voom differential signal on the
**union** peak set, which is the right analysis and is immune to per-genotype
z-scoring. The README does not present peak counts as a KO effect.

**v1.1:** `z_per_replicate` option in `call_peaks_from_bigwigs()` (default
`FALSE`, so v1.0 stays reproducible), `min_width` → 1000. No IgG null is
possible — there is no IgG track in the GEO submission — so the ERCC tracks
should be used as the empirical control arm and the overlap reported.

---

## 2.3 DeepG4 AUROC — confirmed on the framing, overstated on the rest

**Confirmed.** AUROC = 0.563 against GC-matched background, and the v1.0
figure title asserted *"DeepG4 discriminates G4 peaks from GC-matched
background"*. That title cannot stand at 0.563. v1.1 puts the number in the
title, adds a subtitle when AUROC < 0.65, and raises a `warning()`.

**Push-back 1 — ρ = 0.13 was not buried.** The audit says it "belongs in the
README, not a table". But `22_sanity_check.R:3-13` documents *why* DeepG4
cannot correlate with continuous occupancy ("DeepG4 is a *classifier* of
G4-formation propensity, not an occupancy predictor... Phase-1 finding:
Spearman rho ~ 0.13"), and the metrics are suffixed `_SUPPL` deliberately.
That is a documented demotion with a stated reason, not an oversight.

**Push-back 2 — `g4hunter_win` should not be replaced or dropped.** The audit
says "that metric is dead... It should either be replaced with a
max-over-sliding-25-nt G4Hunter score or dropped." The replacement **already
exists in the same table**:

| pair | ρ |
|---|---|
| `g4hunter_pqs` vs `pqsfinder_score` | **0.635** |
| `g4hunter_pqs` vs `deepg4_prob` | **0.330** |
| `g4hunter_win` vs `pqsfinder_score` | 0.00035 |
| `g4hunter_win` vs `deepg4_prob` | 0.021 |

`g4hunter_pqs` is G4Hunter on the PQS subsequence and behaves exactly as it
should. `g4hunter_win` is the same score averaged over 201 bp, and its ρ ≈ 0
next to it *is* the demonstration that a window mean dilutes a ~25 nt motif
into noise — the identical failure mode fixed in `48_strand_asymmetry.R`.
`25_propensity_metrics.R:84-85` already says so.

What *is* wrong is the script header's claim that "High agreement validates all
three", which the table does not support. v1.1 corrects the header, keeps both
metrics, and labels `g4hunter_win` a dilution control on the figure.

**v1.1 also adds** an unmatched-background arm as the audit suggests, and a
realised-vs-target per-GC-bin diagnostic: `18_extract_sequences.R` does
`take <- min(length(avail), target[[b]])`, which truncates silently, so the
delivered background is under-matched exactly at the high-GC end where G4s
live.

---

## 2.4 The `ignore.strand` bug — confirmed, and the damage is wider

`48_strand_asymmetry.R:135` passes a peak GRanges stamped with
`g4_detected_strand` into `findOverlaps()`, which defaults to
`ignore.strand = FALSE`. A `+` peak can only match a `+` gene, so
`gene_strand == g4_detected_strand` is true by construction. Confirmed: all 16
rows of `strand_bias_relative.tsv` read `count_against = 0`, `pct_with = 100`,
p ~ 1e-243. Confirmed that this is the only missing `ignore.strand` in the
repo; all 13 other overlap calls pass it, and lines 67 and 82 are safe because
those peaks carry `*`.

**Three consequences the audit did not list:**

1. `metaprofile_template_vs_coding` is **structurally a single-curve figure**.
   `g4_strand_class` can never take the value `"template"`, so `collect()`
   returns an all-NA profile that is filtered out at line 295.
2. **Peaks were silently discarded, not counted as "against".** Across all 16
   cells only 7,636 of 16,000 peaks (47.7%) survive into the gene-relative
   test; WT promoter 250 keeps 806 of 1,767 (45.6%). The rest matched no
   same-strand gene, returned NA, and vanished.
3. The script header advertises a **Wilcoxon and a chi-square test**. Neither
   exists anywhere in the script. v1.1 implements both.

### Push-back: the prescribed `g4_threshold` fix would empty the table

The audit's diagnosis is right — `g4_threshold: 0.0` makes `has_g4_motif` a
tautology and all 32 rows read 100%. But its fix ("set it to a real G4Hunter
cutoff; 1.2 is the conventional value, 1.5 for stringent") does not work,
because `g4_hunter_score` is a **mean over the whole 501/1001 bp window**, not
a sliding-window max. Measured:

- observed v1.0 scores span **0.088–0.153** (`g4_score_stats.tsv`)
- at 1.2, **0 of 10** literature control G4s pass (`verify/test_g4hunter.R`)
- on 5,000 random mm10 windows, the v1.0 statistic has FPR **0.000** at 1.2

Applying the audit's fix as written would have broken the table in the
opposite direction — from "everything passes" to "nothing passes".

### And the obvious repair is also not enough

Switching to a sliding-25 nt max (what the literature cutoffs are actually
defined against) fixes sensitivity — 9/10 controls pass at 1.2 — but 1.2 is
then far *too permissive* as a per-window classifier, because a max over ~477
overlapping windows is an extreme-value statistic. Measured against real mm10
(`verify/calibrate_g4_threshold.R`):

| extracted width | FPR at 1.2 | FPR at 1.5 | calibrated @ 5% FPR |
|---|---|---|---|
| 501 bp | 0.441 | 0.230 | **2.28** |
| 1001 bp | 0.635 | 0.396 | **2.60** |

The published 1.2/1.5 values are genome-*scanning* thresholds for locating
candidate sites, not per-locus decision thresholds — and the FPR rises with
window width, so no single fixed number can serve both sizes.

**v1.1 therefore calibrates the threshold per window size to a stated false
positive rate**, writes the provenance to
`strand_g4_threshold_calibration.tsv`, and leaves a manual override that
reports the FPR as unknown when used. Reported, not asserted.

**Confirmed and acted on:** `g4_hunter_score(revcomp(s)) == −g4_hunter_score(s)`
exactly (reverse complement maps every G-run of length L to a C-run of length
L; `rle` lengths are preserved, values flip sign, `mean()` is order-invariant).
So `fwd >= rev` was just `fwd >= 0` and `pmax(abs(fwd), abs(rev))` was
`abs(fwd)`. The reverse pass is dropped.

---

## 2.5 `p_perm` — confirmed, but the diagnosis is half the story

**Confirmed:** no shuffling anywhere in the repo. `p_perm` is a `binom.test`,
and the name leaks into four README lines ("Permutation testing confirms this
quantitatively", "permutation p < 0.001"). Renamed `p_binom` in v1.1, with
README replacement text supplied.

**Context the audit omits:** `12_peak_region_enrichment.R:8` documents the
substitution ("binomial replaces regioneR permutation, which hangs on
Windows") and the plot *subtitle* names the binomial null correctly. This is a
naming leak, not a hidden swap. Worth saying, because "no shuffling occurs"
reads as concealment.

**Push-back: the mechanism is only half right.** The audit identifies the
missing peak-width term, correctly. But the region sets are **also never
globally `reduce()`d** — promoters are one TSS ± 2 kb window per gene, so
24,483 × 4,001 bp = **97.9 Mb** of un-merged footprint feeds `p_expect`
(verified: `exp_mean/n_peaks × genome_bp` reproduces 97.93 Mb exactly). That
**inflates** `p_expect` and *deflates* fold — the opposite bias to ignoring
peak width. Fixing only the width term moves the numbers the wrong way.

v1.1 fixes both, with an exact analytic null: a peak of width *w* overlaps
region [s,e] iff its start lies in [s−w+1, e], so the expectation is the
`reduce()`d coverage of the width-dilated region set. Verified against
brute-force simulation to within 1% at every width from 1 bp to 2 kb
(`verify/test_enrichment_null.R`). Runs in seconds, no Windows problem.

Effect on the committed BG4/WT numbers:

| region | fold v1.0 | width factor | fold v1.1 (approx) |
|---|---|---|---|
| promoter | 1.99 | 1.23 | 1.62 |
| enhancer | 5.56 | 2.63 | 2.12 |
| 5′ UTR | **11.90** | 3.71 | **3.21** |
| 1st intron | 1.50 | 1.06 | 1.42 |
| gene body | 1.18 | 1.02 | 1.16 |

The monotone-in-size artefact largely disappears. **This is the artefact that
caused the fold-enrichment figure to be dropped from the README in `a4acb9f`**
— the commit message diagnosed the symptom ("fold enrichment scales inversely
with region size") and treated it as a presentation problem. It was a
misspecified null. The figure is worth reinstating once corrected.

Two points the audit makes that v1.1 documents rather than fixes: peaks are
not independent Bernoulli draws, and with n ≈ 100,000 the p-values do no work
(36 of 60 rows underflow to exactly 0). `regioneR::permTest` with
`circularRandomizeRegions` remains the gold standard on Linux.

---

## 2.6 Spike-in dependency — confirmed, and the README's control claim fails

Confirmed: the WT < FANCJKO < DHX36KO ≈ dKO ordering is read off
*Drosophila*-spike-in-scaled bigWigs, with no independent check and no
internal reference.

**Going further than the audit:** README:465 claims the ERCC tracks being
weaker at promoters shows *"the enrichment is antibody-specific"*. It does
not. The ERCC tracks are **BG4 CUT&Tag with the same antibody**, differing by
batch — a batch control, not an antibody control. And the R-loop side inverts:

| track | S9.6 promoter fold |
|---|---|
| real WT | **0.987** |
| ERCCWT | **1.367** |
| ERCCKO | **1.308** |

The "control" scores *higher* than the real sample. That sentence should be
deleted; replacement text is in `README-deltas.md`. The audit's two suggested
internal controls (background-median normalisation, which
`39_g4_rloop_promoter_correlation.R:bg_level_for` already implements, and a
no-PQS negative locus set) are both good and are specified for a later pass.

---

## 2.7 G4 × R-loop — the z-scoring claim is wrong

> *"z-scoring within genotype guarantees the 'stable across genotypes'
> conclusion. Any global magnitude difference is removed before correlating,
> so 'helicase loss does not change that relationship' is partly true by
> construction."*

**This is incorrect on two independent grounds, and I would defend the code as
written.**

1. **The committed correlations are not computed on the z-scores.**
   `39_g4_rloop_promoter_correlation.R:248-254` passes `log_g4` and `log_rloop`
   to `corr_row()`. The z-scores exist only to put G4 and R-loop on a shared
   plot axis, and the script header at lines 26–28 states this explicitly:
   *"The z-score is monotone-affine, so Pearson r / Spearman rho (reported on
   the log2 fold-over-background values) are UNCHANGED by it."*

2. **Even if they were, it would not matter.** Pearson and Spearman are
   location–scale invariant. z-scoring cannot change `r` under any
   circumstances.

The residual observation — that a within-genotype correlation cannot express a
cross-genotype magnitude shift — is true of every correlation ever computed. It
is not a consequence of z-scoring, and the script already pre-empts it. This
finding should be dropped from the audit, not acted on.

**The second half of 2.7 is fair and unaddressed.** BG4 and S9.6 both require
accessible chromatin, and without an ATAC or IgG covariate the r = 0.51–0.57
cannot be separated from two assays tracking the same open regions. A partial
correlation against promoter ATAC is the right test. Flagged as a data
dependency — there is no ATAC track in the current manifest — and specified
but not implemented in v1.1.

One weak counter-argument worth noting and then discounting: BG4 shows
promoter fold enrichment ~1.99 while S9.6 shows ~0.99, which is not what you'd
expect if both were simply tracking open chromatin. But with 337,757 S9.6
peaks covering 10.5% of the genome, overlap-based enrichment saturates, so
this says more about the peak caller (§2.2) than about specificity. Not a
defence.

---

## 2.8 The expression endpoint is null — confirmed, promote it

Confirmed in full. The headline test (parallel-G4 promoters vs no-G4
promoters) is non-significant in all three KOs: Wilcoxon p = **0.169 / 0.179 /
0.779**, with median |LFC| differences of 0.0019 / 0.0010 / 0.0015. G4 count
per promoter vs log baseMean gives ρ = 0.0217 — significant at p = 0.008 only
because n ≈ 20,000, explaining ~0.05% of variance. All three Kruskal tests on
|LFC| by G4-count bucket are null (p = 0.26–0.58). The one nominal Kruskal
(FANCJKO 0.044) is uncorrected across three tests and sits on a 0.001
difference.

Grep confirms the README never mentions it: `promoter_g4_expression|multiG4|
kruskal|dysregulat` returns only the Stage-K table row at line 163.

Agreed on the framing, and agreed it is the better result. Text supplied in
`README-deltas.md`. The argument is stronger with §2.1's control benchmark
attached: bulk occupancy does not predict transcriptional consequence, *and*
the best available sequence methods cannot agree on topology (κ ≈ 0; 8/10 vs
2/10 on controls). Both point the same way — toward a system where topology,
position and strand are defined on the template rather than inferred.

---

## Findings from this pass that were in neither document

1. **`peak_counts_per_region.tsv` and `region_fold_enrichment.tsv` disagree 4×
   under the same labels.** WT 5′ UTR = **878** (priority-exclusive) vs
   **3,410** (non-exclusive `countOverlaps`). Same directory, same script
   family. Documented in both scripts and in the README Limitations.
2. **README:154** "peaks present in ≥ 2 replicates" misdescribes
   `rowSums(sig_mat > 0)`.
3. **README:158** "201 bp peak-summit sequences" — there is no summit calling
   anywhere in the pipeline.
4. **README:194** advertises "11 literature-validated control G4s"; the file
   has **10**.
5. **`.bw_rle_cache` was keyed on absolute path only**, not on the requested
   seqlevels, so the first call decided the chromosome set for every later call
   on that file. Latent in v1.0 (all callers pass the same `std_chroms`) but a
   trap. Fixed, along with `evict_bigwigs()`, which had to learn the new key or
   it would silently free nothing — and those eviction calls are load-bearing
   for memory, not hygiene.
6. **`unclassified.subsample: 16980`** is a hard-coded magic number equal to
   the observed no-PQS peak count of the v1.0 run.
7. **κ and `topology_source` existed only in uncommitted artefacts** — a
   `message()` into a git-ignored `logs/` and a figure under a git-ignored
   `results/figures/`.

---

## What is deliberately not done here

- **No re-run.** `cache/` is empty, so regenerating means redoing the ~80 min
  pqsfinder sweep plus everything downstream. Scripts are corrected and
  statically verified; execution is a separate decision.
- **`regioneR` permutation** — the analytic width-aware null is exact for
  uniform placement and Windows-safe. A real permutation preserving
  inter-peak distances remains better and is noted in the script header.
- **Accessibility covariate for §2.7** — needs an ATAC track that is not in
  the manifest.
- **IgG/input null for the peak caller** — no such track exists in the GEO
  submission. The ERCC arm is the available substitute.
