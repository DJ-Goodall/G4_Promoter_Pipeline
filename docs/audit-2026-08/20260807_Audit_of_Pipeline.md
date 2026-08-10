
I cloned it and read the code rather than the README. Full audit above; here's the mechanics and the four things that matter.

**How it works.** Everything hangs off one primitive: `read_bigwig_rle()` imports each track as an `RleList` cached in an environment, and every downstream number — peak seeds, region means, metaprofiles — is a `Views`/`viewMeans` over that Rle. Peaks are 500 bp bins → replicate mean → `log1p` → global z → `z ≥ 2` → merge. The topology branch runs pqsfinder on each union peak, feeds the top PQS to a G4ShapePredictor RandomForest, and separately scores a 201 bp PQS-centred window with DeepG4. `cache_or_build()` exploits R's lazy argument evaluation so the 80-minute pqsfinder sweep never re-fires.

**The critical one.** At `precision: 0.80`, G4SP calls 22.3% of PQS "mixed" and only **1.4% hybrid**. But `topology_composition.csv` reports 10.0% hybrid and no mixed bucket — because `21_assign_topology.R` back-fills every mixed call with the loop-length heuristic. So ~83% of your hybrid stratum and ~55% of antiparallel are heuristic calls. And computing Cohen's κ from your own committed concordance table: **κ = −0.0098** (33.5% agreement, 34.2% expected). The two methods agree slightly _worse_ than chance. The "hybrid > parallel > antiparallel" TSS result is currently indistinguishable from a loop-length artefact.

**A real bug.** `48_strand_asymmetry.R:135` — `findOverlaps` defaults to strand-aware, and you pass peaks stamped with the G4 strand. So `gene_strand == g4_detected_strand` is true by construction: every row of `strand_bias_relative.tsv` reads 100%, p ≈ 1e-243. One argument fixes it, and it's the one strand analysis with real content. Separately `g4_threshold: 0.0` makes `g4_detection_by_strand.tsv` 100% in all 32 rows.

**Two framing issues.** `p_perm` in script 12 is an analytic `binom.test` against region_bp/genome_bp — no shuffling, and the null ignores peak width, which is exactly why fold-enrichment inverts with region size. And your `sanity_metrics.csv` has DeepG4 AUROC = **0.563** against GC-matched background, with ρ = 0.13 vs peak strength — that belongs in the README, not a table.

The thing I'd actually push on: `promoter_g4_expression_tests.csv` and `multiG4_expression_stats.csv` are null. Promoter G4 count and topology don't predict which genes move in the KOs (ρ = 0.02; Kruskal p = 0.26–0.58). You've buried it, and I think it's the better result — it says bulk occupancy is insufficient to predict transcriptional consequence, which is the strongest possible motivation for reconstituting this in NPE where you define topology and position instead of inferring them.

---

# Audit — `DJ-Goodall/G4_Promoter_Pipeline` @ `3c1b264`

Clone inspected locally: 52 Snakemake rules, 50 numbered scripts, 647-line shared helper library, 9,657 lines total. Data: Sato et al. GSE269081/82/84 (mESC, mm10; BG4 CUT&Tag ×3 reps, S9.6 ×2 reps, CEL-seq2 RNA-seq ×3, WT / _Dhx36_⁻ᐟ⁻ / _Fancj_⁻ᐟ⁻ / dKO, plus an ERCC batch).

Everything below was checked against the committed `results/tables/` from the v1.0.0 run, so the numbers are the pipeline's own output, not a re-run.

---

## 1. Architecture — how it actually works

### The spine

```
bigWigs ──► RleList (rtracklayer, memoised in .bw_rle_cache)
              │
              ├─ 500 bp genome bins ─► replicate mean ─► log1p ─► global z
              │                            └─► z ≥ 2 seeds ─► reduce(gap 500) ─► peaks
              │
              └─ IRanges::Views ─► viewMeans / binnedAverage ─► signal + metaprofiles
```

Three design decisions carry almost all the weight:

**(a) Signal extraction is `Views` over an in-memory `RleList`.** `read_bigwig_rle()` (`helpers.R:117`) imports each whole bigWig as an RleList and caches it in an environment keyed by absolute path; `evict_bigwigs()` frees it. Every downstream quantity — peak seeds, per-region means, profile matrices — is a `viewMeans`/`binnedAverage` over that Rle. This is why the pipeline is fast and also why memory is the binding constraint. The eviction calls after each genotype loop are load-bearing, not hygiene.

**(b) `compute_profile_matrix()` (`helpers.R:180`) is the metaprofile engine.** Per chromosome it builds a Views object over all windows, `rbind`s them into a regions × bp matrix, bins with `cut(seq_len(w), breaks = n_bins)`, then reverses rows for `-`-strand features. The `rbind` is only safe because every caller passes fixed-width windows via `center_window()` — I checked all 13 call sites and they all do. If anyone ever passes variable-width regions (gene bodies) directly, `rbind` will silently recycle short rows and produce garbage rather than error.

**(c) `cache_or_build()` relies on R's lazy argument evaluation.** `cache_or_build(file, expr)` returns `readRDS(file)` before `expr` is ever forced. That is correct but fragile — it is the reason an 80-minute `pqsfinder` sweep doesn't run on every invocation. The PQS cache is additionally keyed on peak count (`peak_pqs_n%d.rds`), which is a nice touch: the smoke-test subsample can't contaminate a full run.

### The topology branch (scripts 17–21, the intellectual core)

```
union peaks ──► pqsfinder (min_score 20, top PQS per peak)
                    │                    │
                    │                    └─► loop lengths ll1/ll2/ll3 ─► heuristic
                    │
                    ├─► PQS subsequence ─► G4ShapePredictor RF (precision 0.80)
                    │                          → parallel / antiparallel / hybrid / mixed
                    │
                    └─► 201 bp window centred on the PQS ─► DeepG4 → P(G4 forms)
                                                              (sanity/ROC only)
```

`21_assign_topology.R` merges them with this precedence:

```r
topology_final = case_when(
  no_pqs        ~ "no_canonical_PQS",
  g4sp_definite ~ topology_g4sp,          # G4SP won
  !is.na(topology_heuristic) ~ topology_heuristic,   # <- the "mixed" bucket lands here
  TRUE ~ "no_canonical_PQS")
```

That last-but-one line is where the biggest problem in the repository lives (§2.1).

---

## 2. Findings, ranked by how much they change the conclusions

### 2.1 The topology strata are a chimera, and the two methods agree at chance level — **critical**

`topology_precision_sweep.csv` shows what G4ShapePredictor alone calls at the configured `precision: 0.80`, over the 130,863 union peaks with a PQS:

||parallel|antiparallel|hybrid|**mixed**|
|---|---|---|---|---|
|G4SP @ 0.80|95,660 (73.1 %)|4,099 (3.1 %)|1,879 (**1.4 %**)|29,225 (**22.3 %**)|

`topology_composition.csv` — what the figures actually use — has no `mixed` category:

||parallel|antiparallel|hybrid|no PQS|
|---|---|---|---|---|
|final label (mean over genotypes)|79,338 (72.9 %)|8,643 (7.9 %)|10,934 (10.0 %)|10,077 (9.1 %)|

The 29,225 peaks G4SP explicitly declined to classify were reassigned by the loop-length heuristic. Arithmetically that means:

- **~83 % of the "hybrid (3+1)" stratum** and **~55 % of the "antiparallel (2+2)" stratum** are loop-length calls, not model calls.
- The parallel class is the only one that is predominantly G4SP.

That would be defensible if the heuristic were a validated stand-in. It is not. Computing Cohen's κ from the committed `topology_concordance.csv` (n = 101,638 peaks where both methods gave a definite call):

```
raw agreement = 33.5 %   expected by chance = 34.2 %   κ = −0.0098
```

**The two methods agree slightly worse than chance.** Their marginals are also wildly different (G4SP: 94 % parallel; heuristic: 34 % parallel / 40 % hybrid / 25 % antiparallel). The script computes κ and prints it into the concordance figure subtitle, but the value never surfaces in the README, which describes the heuristic as a "cross-check."

**Consequence.** The headline topology result — _hybrid > parallel > antiparallel TSS peak height_ (`10_metaprofiles_tss_signal_by_topology.png`) — is largely a statement about "PQS with exactly one loop longer than 3 nt," which is not the same object as a hybrid G4. The ordering may still be real; the current evidence can't distinguish it from a loop-length artefact.

**Suggested fixes**, in order of effort:

1. Keep `mixed` as its own stratum. Set the fallback to `"mixed"` rather than the heuristic, and plot it. If the hybrid TSS peak survives on the 1,879 G4SP-definite hybrids, that's a real result.
2. Report κ in the README, and report `topology_source` composition (`g4sp` vs `heuristic`) alongside every topology figure.
3. Consider `precision: 0.75` (6,090 mixed, 4.7 %) as the working default, with 0.80 as a sensitivity analysis. The sweep table already gives you the data to justify the choice.

### 2.2 Peak calling is a signal quantile, not an enrichment test — **critical**

`call_peaks_from_bigwigs()` (`helpers.R:523`): bin the genome at 500 bp, take the replicate mean, `log1p`, z-score **globally within that genotype**, keep `z ≥ 2`, merge across 500 bp gaps.

The output (`peak_summary.csv`):

|set|peaks|mean bp|total Mb|% of mm10|
|---|---|---|---|---|
|BG4 union|147,843|966|142.8|**5.2 %**|
|S9.6 union|337,757|848|286.4|**10.5 %**|

For comparison, the source paper called peaks with MACS3 (q < 0.05) requiring consensus across all three replicates. This pipeline reimplements a different, far more permissive caller and the README frames it as "a genuine reproduction of the calling step" — it is a reproduction of _a_ calling step, not of Sato et al.'s. Three specific problems:

- **No input/IgG normalisation.** z ≥ 2 on a log1p-transformed distribution selects roughly the top 2.3 % of bins whatever the data are. There is no null model in which a peak can fail to be called.
- **`min_replicates` is close to a no-op.** The replicate criterion is `rowSums(sig_mat > 0) >= min_replicates` — signal strictly greater than zero in ≥2 replicates. In a bin already at z ≥ 2 this is almost always satisfied. The z-score itself is computed on the replicate _mean_, so a single strong replicate can carry a bin. The README's "replicate-aware peak calling" and "peaks present in ≥2 replicates" overstate this. A real version would z-score each replicate separately and require `z ≥ 2` in ≥2 of them.
- **Cross-genotype peak counts are not comparable by construction.** μ and σ are recomputed per genotype, so a global increase in G4 occupancy in a KO is normalised away before thresholding. This is fine — arguably desirable — but it means `peak_summary.csv` counts (107k WT vs 114k dKO) carry no information about the KO effect, and shouldn't be read as if they did.

The independent evidence that the peak set is noisy is already in the repo: see §2.3.

### 2.3 DeepG4 AUROC = 0.563 — the sanity check fails, and it isn't flagged

`sanity_metrics.csv`:

```
auroc_peak_vs_background   0.563
median_prob_peak           0.088
median_prob_background     0.061
spearman_prob_vs_maxz      0.131
```

DeepG4 barely separates called BG4 peaks from GC-matched random genomic windows, and its predicted G4-formation probability is essentially uncorrelated with peak strength. Two non-exclusive readings:

1. The GC-matching (`18_extract_sequences.R`, quantile-binned to the peak GC distribution) is so aggressive that it removes the signal — GC content is the dominant determinant of both DeepG4 score and BG4 occupancy, so matching on it matches away the biology. This is a real possibility and worth testing with an unmatched background as a second arm.
2. A large fraction of the 147,843 peaks are not G4 sites. Given §2.2, this is likely a substantial contributor.

Either way, an AUROC of 0.56 on your own positive control should be prominent, not a line in `sanity_metrics.csv`. Right now the README calls DeepG4 "the one fragile dependency" whose failure "costs a single figure" — the framing is about whether TensorFlow loads, not about what the result says.

Relatedly, `propensity_crossvalidation.csv` shows `g4hunter_win` vs `pqsfinder_score` ρ = 0.00035 and `g4hunter_win` vs `deepg4_prob` ρ = 0.021. That metric is dead: G4Hunter averaged over a 201 bp window sums signed G- and C-run contributions, which cancel to near zero over any window much larger than the motif. It should either be replaced with a max-over-sliding-25-nt G4Hunter score or dropped.

### 2.4 `strand_bias_relative` is an artefact of a missing `ignore.strand` — **bug**

`48_strand_asymmetry.R:130–142`:

```r
gr_peaks <- GRanges(..., strand = ifelse(g4_detected_strand %in% c("+","-"),
                                          g4_detected_strand, "*"))
gs_ov <- findOverlaps(gr_peaks, genes_gr, select = "first")   # <- strand-aware by default
sequences_df$gene_strand <- as.character(strand(genes_gr))[gs_ov]
sequences_df$strand_class <- ifelse(g4_detected_strand == gene_strand,
                                    "with_strand", "against_strand")
```

`findOverlaps` on two stranded GRanges defaults to `ignore.strand = FALSE`. A peak stamped `+` can therefore only match a `+` gene, so `gene_strand == g4_detected_strand` is **true by construction**. The committed table shows exactly that:

```
WT  promoter 250   count_with = 806   count_against = 0   pct_with = 100   p = 4.7e-243
```

Every row: 100 %, p ~ 1e-240. Fix is one argument: `findOverlaps(gr_peaks, genes_gr, select = "first", ignore.strand = TRUE)`. This is the only place in the codebase where `ignore.strand` is missing — lines 67 and 82 of the same script are fine because those peaks still carry `*`. Worth re-running: the coding-vs-template asymmetry is the one strand result with real biological content, and right now you don't have it.

Two related issues in the same script:

- `g4_threshold: 0.0` makes `has_g4_motif = pmax(abs(fwd), abs(rev)) >= 0` always true. Hence `g4_detection_by_strand.tsv` is 100 % in all 32 rows — a table that cannot report anything. Set it to a real G4Hunter cutoff (1.2 is the conventional value; 1.5 for stringent).
- `strand_bias_summary` (reference `+` vs `-`) is a null test by design: promoters are ~evenly split between `+` and `-` genes, so 50/50 is the expected answer. The result (padj 0.79–0.96) is correct and uninformative. Only the gene-relative version is meaningful — which is the one broken by the bug above.

Note that `g4_hunter_score(revcomp(s)) == -g4_hunter_score(s)` identically under this implementation, so `fwd >= rev` reduces to `fwd >= 0`. That's a fine definition of "which strand is G-rich," but it means `rev` is a redundant computation and `pmax(abs(fwd), abs(rev))` is just `abs(fwd)`.

### 2.5 "Permutation-based fold enrichment" is an analytic binomial test — **naming + null model**

`12_peak_region_enrichment.R:80–90`. The column is called `p_perm`, the README says "permutation p < 0.001" and "Permutation testing confirms this quantitatively," but the code is:

```r
p_expect <- sum(width(rg)) / genome_bp
obs      <- sum(countOverlaps(pk, rg) > 0)
p_perm   <- binom.test(obs, n_peaks, p_expect, alternative = "greater")$p.value
```

No shuffling occurs. Beyond the naming, the null is mis-specified for interval overlap in three ways:

1. **Peak width is ignored.** A 966 bp peak overlaps a region set with probability ≈ (region_bp + n_regions × peak_width) / genome_bp, not region_bp / genome_bp. Expected overlap is underestimated, fold enrichment inflated. This is also the mechanism behind the effect the README already noticed — "fold-enrichment values scale inversely with region size" — small features (5′ UTRs) are hit hardest.
2. **Peaks are not independent Bernoulli draws.** They cluster with mappability, GC and accessibility.
3. **n ≈ 100,000** means any deviation returns p ≈ 0. The p-values are not doing work.

`regioneR::permTest` with `circularRandomizeRegions` (preserves inter-peak distances) or shuffling within an accessibility-matched universe would give an honest answer, and would probably shrink the enrichment considerably. The v1.0.0 commit history shows the fold-enrichment plot was dropped from the README showcase (`a4acb9f`), which suggests you already had doubts — the underlying table is still in `rule all` and still labelled as a permutation.

### 2.6 The magnitude claim rests entirely on someone else's spike-in scaling

The headline occupancy result — WT < FANCJKO < DHX36KO ≈ dKO at promoters — is read off raw bigWig values. Those bigWigs are the _Drosophila_-spike-in-scaled tracks from Sato et al.'s Snakemake workflow. The pipeline contains no independent check on that scaling and no internal reference (no IgG, no unchanged-region normaliser). If the spike-in factors are off, the entire genotype ordering moves.

This is not a criticism of the analysis so much as an unstated dependency: the pipeline can reproduce the ordering, but it cannot validate it. Two cheap internal controls are available:

- Normalise each track to its own background median (you already implement exactly this in `39_g4_rloop_promoter_correlation.R:bg_level_for`) and show the promoter metaprofile ordering survives. If it does, the claim is much stronger.
- Show the ordering at a set of loci with no PQS as a negative control. If WT < KO holds there too, it's a scaling/depth effect.

The ERCC batch clustering separately at r ≈ 0.70–0.80 while within-batch replicates hit r ≈ 0.94–0.98 demonstrates that batch effects of the relevant magnitude do exist in this assay.

### 2.7 The G4 × R-loop correlation may be an open-chromatin correlation

`39_g4_rloop_promoter_correlation.R` is the most carefully built module in the repo — shared promoter subsample across genotypes, per-replicate background-median depth normalisation, z-scoring before the loess. Result: Pearson r = 0.51–0.57, with a threshold-like loess that only rises above ~2.5 z.

Two caveats:

- **z-scoring within genotype guarantees the "stable across genotypes" conclusion.** Any global magnitude difference is removed before correlating, so "helicase loss does not change that relationship" is partly true by construction. Worth stating.
- **No accessibility control.** BG4 and S9.6 CUT&Tag both require accessible chromatin and both concentrate at active TSSs. Without an ATAC/IgG covariate you can't separate genuine G4–R-loop coupling from two assays independently tracking the same open regions. A partial correlation against promoter ATAC signal, or against the ERCC-batch tracks as an antibody-nonspecific reference, would address this. The non-linear threshold shape is the most interesting part of the result and also the part most vulnerable to this confound.

### 2.8 The biological endpoint is null, and the README doesn't say so

This is the finding I'd most want you to sit with, because it bears directly on the Veni framing.

`promoter_g4_expression_tests.csv`:

|KO|Kruskal p (topology vs \|LFC\|)|median \|LFC\| parallel|median \|LFC\| no G4|
|---|---|---|---|
|DHX36KO|0.112|0.1277|0.1259|
|FANCJKO|0.044|0.0256|0.0246|
|dKO|0.308|0.1763|0.1748|

`multiG4_expression_stats.csv`:

```
spearman nG4 vs log baseMean   rho = 0.0217   p = 0.0078      # significant, effect ≈ 0
kruskal |LFC| by G4-count bucket   p = 0.58 / 0.31 / 0.26     # null
```

So: **promoter G4 content — neither the number of G4s nor their assigned topology — predicts which genes change expression when the resolving helicases are removed.** The one nominally significant Kruskal (FANCJKO, p = 0.044) has a median |LFC| difference of 0.001 and would not survive correction across the three tests.

The README's "Selected results" section presents the occupancy gains and the topology-resolved TSS peaks, and stops there. It never states the null. I'd argue the null is the more publishable finding, and it is a _better_ motivation for your NPE work than a weak positive would be: it says that steady-state G4 occupancy at a promoter, measured in bulk chromatin, is not sufficient to predict transcriptional consequence — which is exactly the argument for needing a reconstituted system where you control structure, topology and position independently. Consider promoting it.

---

## 3. What the pipeline does well

Not everything here is a caveat, and several things are better than typical:

- **Provenance discipline.** One rule per numbered script, one log per rule, script headers document inputs/outputs/config keys consumed, every seed in `config/config.yaml`, `sessionInfo()` in logs. The `peak_pqs_n%d.rds` cache-keying is a genuinely good idea that prevents a whole class of silent staleness bugs.
- **The G4SP reimplementation is careful.** `20_g4sp_topology.py` reproduces the upstream GUI's exact inference path (`seqs_converter` → `predict_proba` → `get_threshold`) rather than guessing, and the scikit-learn 1.0.2 pin with the documented reason (tree node dtype changed in 1.3, newer versions can't unpickle and would silently fall back) is exactly the kind of thing that usually goes unrecorded.
- **PQS-centred DeepG4 windows.** Centring the 201 bp window on the top PQS rather than the peak midpoint (`18_extract_sequences.R`) is correct and non-obvious given 500 bp – 5.5 kb peaks.
- **Graceful degradation is real,** and the failure modes are documented with the diagnostic to check.
- **Honest negative results are committed** — the strand-bias nulls, the near-zero propensity correlations, the flat expression tests are all in `results/tables/`. Nothing is hidden; the issue is emphasis, not integrity.

---

## 4. Priority list

|#|Action|Effort|Changes conclusions?|
|---|---|---|---|
|1|Stop back-filling `mixed` with the heuristic; keep it as a stratum; report κ|~1 h|**Yes** — topology figures|
|2|Fix `ignore.strand = TRUE` at `48_strand_asymmetry.R:135`; set `g4_threshold` > 0|~10 min|**Yes** — gives you a real result|
|3|Rename `p_perm` → `p_binom`, or implement a real `regioneR` permutation|10 min / 1 day|Yes — enrichment magnitudes|
|4|Add a per-replicate z threshold so `min_replicates` bites; report peak counts at MACS3-comparable stringency as a sensitivity arm|~half day|Possibly — everything downstream|
|5|Background-median-normalise the promoter metaprofiles as a check on spike-in scaling|~2 h|Confirms or breaks the headline|
|6|Add an unmatched-background DeepG4 AUROC arm; drop or fix `g4hunter_win`|~1 h|Interprets the 0.563|
|7|Promote the expression null into "Selected results"|~20 min|Framing|

---

## 5. Relevance to the transcription-competent NPE work

Three things transfer directly:

- The **topology-assignment problem is unsolved computationally**, and this audit quantifies how unsolved (κ ≈ 0 between the two best available sequence-based methods). That is a strong argument for the NPE approach: you can _define_ topology on the template by construction — parallel, antiparallel and hybrid inserts with 7-deaza-dG ablation controls — instead of inferring it. The pipeline's failure is your proposal's premise.
- The **null expression result** says bulk promoter G4 occupancy doesn't predict transcriptional outcome in cells. The obvious explanations are that occupancy conflates formation with persistence, that position relative to the TSS matters more than presence, and that the relevant variable is strand. All three are exactly what a defined plasmid template lets you separate, and the second is already your positional-rules aim.
- The **G4 × R-loop threshold behaviour** (flat below ~2.5 z, steep above) is the most mechanistically suggestive result in the repo and survives most of the criticisms above. If real, it implies coupling is cooperative rather than proportional — which is what a G-loop model predicts, and which you could test in NPE by titrating transcript and asking whether R-loop-dependent G4 stabilisation is switch-like.
