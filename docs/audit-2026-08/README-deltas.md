# README changes for v1.1

Exact replacement text for every README claim the audit found, plus the ones
found while checking it. Line numbers are against `README.md` at `3c1b264`.

---

## 1. The "permutation" claims — four places, all wrong

No shuffling has ever happened in `12_peak_region_enrichment.R`. The script
header documented the substitution ("binomial replaces regioneR permutation,
which hangs on Windows") and the plot subtitle named the binomial null
correctly, so this was a naming leak rather than a hidden swap — but it leaked
into four reader-facing lines.

**Line 129** (mermaid diagram)

> `permutation fold-enrichment · signal · metaprofiles`

→

> `binomial fold-enrichment · signal · metaprofiles`

**Line 156** (stage table)

> Permutation-based fold enrichment of peaks in each region class

→

> Width-aware binomial fold enrichment of peaks in each region class

**Lines 463–469** — replace the paragraph wholesale:

> Fold enrichment is quantified against a width-aware binomial null
> (`12_peak_region_enrichment.R`, table `region_fold_enrichment.tsv`). The
> expectation for each peak is the reduce()d region footprint dilated by the
> peak's own width — a 966 bp peak overlaps a region set whenever its start
> falls within `region_width + peak_width − 1`, not merely within the region.
> This is an analytic test, **not a permutation**: `regioneR::permTest` with
> `circularRandomizeRegions` would be the gold standard but hangs on Windows,
> which is why the binomial was substituted. BG4 peaks are enriched at
> promoters, distal enhancers and 5′ UTRs in every genotype.
>
> Two caveats. First, with n ≈ 100,000 peaks any deviation returns p ≈ 0
> (36 of 60 rows underflow to exactly zero), so the p-values are not doing
> work — read the fold values, not the stars. Second, peaks are not
> independent Bernoulli draws; they cluster with mappability, GC and
> accessibility, which the binomial cannot express.

Note the v1.0 sentence about fold enrichment scaling inversely with region
size can now be **deleted**: that artefact was the missing peak-width term,
and correcting the null largely removes it (5′ UTR falls from ~11.9 toward the
promoter value rather than dwarfing it). The fold-enrichment figure dropped
from the showcase in `a4acb9f` is worth reinstating once the corrected numbers
are in.

---

## 2. The ERCC "antibody-specific" claim — delete it

**Line 465–466**

> while the ERCC1 control tracks are far weaker at promoters — the enrichment
> is antibody-specific, not an artefact of the peak caller

This does not support the conclusion. The ERCC tracks are **BG4 CUT&Tag with
the same antibody**, from a different batch — a batch control, not an antibody
control. And the R-loop side inverts: ERCCWT promoter fold is **1.367** and
ERCCKO **1.308**, both *higher* than real S9.6 WT at **0.987**. Replace with:

> The ERCC batch scores lower at promoters (BG4 fold 1.27–1.45 vs 1.92–2.12),
> but these tracks use the same BG4 antibody and differ from the main panel by
> batch, so this is a batch comparison rather than an antibody-specificity
> control. The pipeline has no IgG or input track and cannot establish
> antibody specificity. Note also that the ERCC **R-loop** tracks score
> *higher* at promoters (1.31–1.37) than the real S9.6 WT track (0.99), which
> is the opposite of what a specificity argument would predict.

---

## 3. Peak-calling claims — two places

**Line 13**

> replicate-aware peak calling

→

> z-score peak calling

(or keep "replicate-aware" only if `peak_calling.z_per_replicate: true` is set,
in which case it becomes accurate for the first time.)

**Line 154** (stage table) — this describes something the code does not do:

> 500 bp genome bins → replicate-mean signal → global z-score (z ≥ 2) → peaks
> present in ≥ 2 replicates, merged across 500 bp gaps

The v1.0 replicate criterion was `rowSums(sig_mat > 0) >= 2` — non-zero
*coverage* in ≥ 2 replicates, which at CUT&Tag depth is close to a no-op — and
the z-test was applied once to the replicate mean, so one strong replicate
could carry a bin. Replace with either the honest v1.0 description:

> 500 bp genome bins → replicate-mean signal → global z-score (z ≥ 2), with
> non-zero coverage required in ≥ 2 replicates, merged across 500 bp gaps

or, with `z_per_replicate: true`, the accurate v1.1 one:

> 500 bp genome bins → per-replicate global z-score → bins passing z ≥ 2 in
> ≥ 2 replicates, merged across 500 bp gaps

**Line 65** — reword the ambiguity:

> so peak counts are a genuine reproduction of the calling step, not a copy of it

The intent (recomputing rather than downloading a deposited peak file) is
defensible, but the sentence reads as a claim to reproduce Sato et al.'s MACS3
calls, which this is not. Replace:

> so peaks are recomputed here rather than copied from a deposited peak file.
> Note this is a *different* caller from the source paper, which used MACS3
> (q < 0.05) with consensus across all three replicates; the z-score caller
> here is more permissive (BG4 union covers 5.2 % of mm10, S9.6 union 10.5 %)
> and the counts are not comparable to theirs.

---

## 4. Topology claims

**Line 158** — two errors in one row.

> 201 bp peak-summit sequences + GC-matched background

There is no summit calling anywhere in the pipeline; peaks are merged 500 bp
z-score bins. Windows are centred on the top PQS, with the peak midpoint as
fallback. →

> 201 bp PQS-centred sequences (peak midpoint as fallback) + GC-matched background

> loop-length heuristic cross-check

It is not a cross-check in v1.0 — it *replaced* G4SP's judgment for 22.3 % of
peaks. In v1.1 it is genuinely only a comparator. →

> loop-length heuristic reported alongside for comparison (not used to assign labels)

**Line 194**

> 11 literature-validated control G4s

`data/known_topology_controls.tsv` contains **10**. →

> 10 literature-validated control G4s

---

## 5. New Limitations section

Add before "Selected results". These are the numbers a reader needs and none
of them appear in the v1.0 README.

> ## Limitations
>
> **Topology assignment is unresolved, and v1.1 says so.** At the configured
> precision G4ShapePredictor declines to classify 22.3 % of PQS-bearing peaks.
> v1.0 back-filled those with a pqsfinder loop-length heuristic, which
> supplied 87 % of the reported hybrid stratum and 66 % of the antiparallel
> stratum. Scored on the same 10 literature control quadruplexes, G4SP gets
> **8/10** and the loop heuristic **2/10** — below the 33 % expected from
> guessing at a three-class problem
> (`topology_control_benchmark_compare.csv`). v1.1 therefore keeps G4SP's
> abstentions as their own `mixed` stratum and makes no substitutions. Call
> provenance is reported in `topology_source_composition.csv`.
>
> **The two methods agree at chance level.** Cohen's κ = −0.0098 over the
> 101,638 peaks where both give a definite call. κ is depressed here by very
> different marginals (G4SP 94.1 % parallel, heuristic 34.4 %), so the
> marginal-free statement is the fairer one: of the heuristic's 41,122 hybrid
> calls only 741 (1.8 %) are also G4SP-hybrid; of its 25,539 antiparallel
> calls only 191 (0.7 %) agree.
>
> **DeepG4 AUROC is 0.563** against GC-matched background, with Spearman
> ρ = 0.13 against peak strength. DeepG4 is a classifier of formation
> propensity rather than an occupancy predictor, so the weak ρ is expected —
> but an AUROC this close to 0.5 on the positive control is not. Either the GC
> matching removes the signal (GC drives both DeepG4 score and BG4 occupancy),
> or a substantial fraction of called peaks are not G4 sites. The DeepG4
> panels are diagnostic, not validation.
>
> **The peak caller has no null.** A relative global z ≥ 2 on each sample's own
> signal cannot fail to return peaks: the ERCC control tracks yield *more*
> peaks than the real samples (114,980 vs 107,276 for WT). Peak counts are not
> comparable across genotypes either, since μ and σ are recomputed per
> genotype — a global occupancy change is normalised away before thresholding.
> Read the limma-voom differential signal on the union peak set
> (`23_gain_loss.R`) for the KO effect, not `peak_summary.csv`.
>
> **The magnitude claim depends on someone else's spike-in scaling.** The
> WT < FANCJKO < DHX36KO ≈ dKO ordering is read off *Drosophila*-spike-in-scaled
> bigWigs from Sato et al. This pipeline contains no independent check on that
> scaling and no internal reference. The ERCC batch clustering at r ≈ 0.70–0.80
> while within-batch replicates reach r ≈ 0.94–0.98 shows batch effects of the
> relevant magnitude exist in this assay.
>
> **`peak_counts_per_region.tsv` and `region_fold_enrichment.tsv` count
> differently** under the same region labels. The first is priority-exclusive
> (a peak is assigned to promoter > enhancer > 5′ UTR > first intron > gene
> body and counted once): WT 5′ UTR = 878. The second uses non-exclusive
> overlap: the same cell is 3,410. Neither is wrong; they answer different
> questions.

---

## 6. Promote the null expression result into "Selected results"

Currently absent from the README narrative entirely — grep for
`promoter_g4_expression|multiG4|kruskal` returns only the Stage-K table row at
line 163. Add:

> ### Promoter G4 content does not predict which genes respond to helicase loss
>
> Neither the number of G4s at a promoter nor their assigned topology predicts
> which genes change expression in the knockouts
> (`29_promoter_g4_expression.R`, `47_multiG4_promoter_expression.R`). The
> headline test — parallel-G4 promoters vs no-G4 promoters — is
> non-significant in all three KOs (Wilcoxon p = 0.169 / 0.179 / 0.779), with
> median |LFC| differences of ~0.001. G4 count per promoter correlates with
> baseline expression at ρ = 0.022 (p = 0.008, significant only because
> n ≈ 20,000) and with |LFC| not at all (Kruskal p = 0.26–0.58). The one
> nominally significant result (FANCJKO Kruskal p = 0.044) sits on a median
> difference of 0.001 and would not survive correction across the three tests.
>
> This is a real result, not a failed one. Steady-state G4 occupancy at a
> promoter, measured in bulk chromatin, is **not sufficient** to predict
> transcriptional consequence. Occupancy conflates formation with persistence;
> position relative to the TSS is collapsed into presence/absence; and strand
> is not resolved. Separating those requires a system where topology, position
> and strand are defined on the template rather than inferred from sequence —
> which is precisely what the κ ≈ 0 and 2/10 control results above say cannot
> currently be done computationally.
