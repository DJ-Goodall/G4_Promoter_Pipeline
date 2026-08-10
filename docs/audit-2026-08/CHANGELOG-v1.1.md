# Changelog — v1.1 (internal review build)

Semantic versioning in this repo says the **major** version changes when the
scientific output changes. By that rule this is not a 1.1 — the topology
strata change, so a public release of these fixes would be **2.0.0**. It is
labelled v1.1 here because it is an internal review overlay, not a release.

Not committed, not pushed. Depends on the v1.0 tree for every unchanged file.

---

## Changed — results move

- **Topology labels no longer mix two methods.** `21_assign_topology.R` no
  longer back-fills G4ShapePredictor's abstentions with the pqsfinder
  loop-length heuristic. Peaks G4SP declines to classify are labelled `mixed`
  and stay there. At `precision: 0.80` this affects 29,225 of 130,863
  PQS-bearing union peaks, which had been supplying **87.0%** of the reported
  hybrid stratum and **66.1%** of the antiparallel stratum.
  → The "hybrid > parallel > antiparallel" TSS ordering must be re-derived.

- **Gene-relative strand analysis now returns a result.**
  `48_strand_asymmetry.R:135` was missing `ignore.strand = TRUE`, making
  `gene_strand == g4_detected_strand` true by construction: all 16 rows of
  `strand_bias_relative.tsv` read 100% / p ~ 1e-243, and ~52% of peaks were
  silently discarded rather than counted as against-strand. Also unblocks
  `metaprofile_template_vs_coding`, which could previously only ever plot one
  of its two curves.

- **G4Hunter switched from window mean to sliding-25 nt max**, and the motif
  threshold is now calibrated per window size against random mm10 at a stated
  false-positive rate (~2.28 at 501 bp, ~2.60 at 1001 bp for 5% FPR) instead of
  the tautological `g4_threshold: 0.0`.

- **Region fold-enrichment null corrected for peak width and un-merged
  regions.** Expected overlap is now the `reduce()`d region footprint dilated
  by the peak's own width. Largely removes the inverse-region-size artefact
  that caused the fold-enrichment figure to be dropped in `a4acb9f`
  (5′ UTR ~11.9 → ~3.2; promoter 1.99 → 1.62).

- **`peak_calling.z_per_replicate`** (default `FALSE`, so v1.0 reproduces
  exactly). When `TRUE`, each replicate is z-scored independently and the
  threshold must be met in ≥ `min_replicates` of them — what the README always
  claimed. `min_width` 500 → 1000, since 500 == `bin_size` made it inert.

## Changed — naming and framing only

- `region_fold_enrichment.tsv`: column `p_perm` → **`p_binom`**. No shuffling
  has ever happened in this script. Plot title no longer says "relative to
  random placement".
- `22_sanity_check.R`: the ROC figure title states the AUROC (0.563) instead of
  asserting that DeepG4 "discriminates"; warns when AUROC < 0.65.
- `25_propensity_metrics.R`: header no longer claims "high agreement validates
  all three"; `g4hunter_win` is labelled a dilution control.

## Added

- **`26d_heuristic_control_benchmark.R`** — scores the loop-length heuristic on
  the same 10 literature control G4s as G4SP. Result: **G4SP 8/10, heuristic
  2/10**. This is the evidence that justifies removing the back-fill.
- `topology_concordance_stats.csv` — κ (−0.0098) plus the marginal-free
  cross-classification rates (1.8% / 0.7%), which are the fairer summary given
  the divergent marginals. Previously κ reached no committed artefact at all.
- `topology_source_composition.csv` + figure — the `topology_source` column
  already existed in v1.0 but was never aggregated anywhere.
- `strand_g4_threshold_calibration.tsv`, `strand_score_tests.tsv` (Wilcoxon),
  `strand_detection_tests.tsv` (chi-square). The latter two were advertised in
  the v1.0 script header but never implemented.
- `background_gc_matching.csv` + `cache/background_seqs_unmatched.fa` — the
  unmatched-background arm, and a realised-vs-target GC diagnostic for a
  matching step that truncated silently.
- `"mixed"` added to `topology_levels` / `topology_labels` /
  `topology_palette`; `topology_definite` exported so scripts stop redefining
  it locally.

## Fixed — latent

- `.bw_rle_cache` was keyed on absolute path only, ignoring the requested
  seqlevels, so the first call fixed the chromosome set for every later call on
  that file. Now keyed on path + seqlevels, with `evict_bigwigs()` updated to
  match — without which eviction would silently free nothing, and those calls
  are load-bearing for memory.

## Documentation

- `README-deltas.md` — replacement text for 9 claims: the four "permutation"
  lines, "replicate-aware", "peaks present in ≥ 2 replicates", "peak-summit
  sequences", "11 literature-validated controls" (there are 10), "genuine
  reproduction of the calling step", and the unsupported ERCC
  "antibody-specific" sentence. Plus a Limitations section and the promoted
  null expression result.
- `REVIEW.md` — verdict on each audit point, including the three places I
  disagree.

## Not done

No re-run. `cache/` is empty, so regenerating means redoing the ~80 min
pqsfinder sweep and everything downstream. Also deferred: a real `regioneR`
permutation (Linux only), an accessibility covariate for the G4 × R-loop
correlation (no ATAC track in the manifest), and an IgG null for the peak
caller (no IgG track exists in the GEO submission).
