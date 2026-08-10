# Changelog

All notable changes to this workflow are recorded here. Versions follow
[semantic versioning](https://semver.org/): the major version changes when the
scientific output changes, the minor version when rules or parameters are added,
and the patch version for documentation and fixes that leave results unchanged.

## [Unreleased] — targeting 2.0.0

Response to an external audit of `3c1b264` dated 2026-08-07. The point-by-point
response, including the three findings I disagree with, is in
[`docs/audit-2026-08/REVIEW.md`](docs/audit-2026-08/REVIEW.md); every quantitative
claim in it is re-derived by [`verify/run_all.sh`](verify/).

**This is a MAJOR change by the convention stated above** — the topology strata
change, so the scientific output changes. The branch is dated rather than
numbered because these corrections have not yet been run end to end.

### Changed — results move

- **Topology labels no longer mix two methods.** `21_assign_topology.R` and
  `43_assign_motif_topology.R` no longer back-fill G4ShapePredictor's
  abstentions with the pqsfinder loop-length heuristic. Peaks G4SP declines to
  classify are labelled `mixed` and stay there. At `precision: 0.80` this
  affects 29,225 of 130,863 PQS-bearing union peaks, which had been supplying
  **87.0 %** of the reported hybrid stratum and **66.1 %** of the antiparallel
  stratum. The "hybrid > parallel > antiparallel" TSS ordering must be
  re-derived. Both scripts had the same back-fill; correcting only the peak
  level would have left `topology_composition_motif_vs_peak` comparing a
  de-chimerised arm against a still-chimeric one.
- **Gene-relative strand analysis now returns a result.**
  `48_strand_asymmetry.R` was missing `ignore.strand = TRUE` on one
  `findOverlaps`, making `gene_strand == g4_detected_strand` true by
  construction: all 16 rows of `strand_bias_relative.tsv` read 100 % /
  p ~ 1e-243, and ~52 % of peaks were silently discarded rather than counted as
  against-strand. Also unblocks `metaprofile_template_vs_coding`, which could
  previously only ever plot one of its two curves.
- **G4Hunter switched from window mean to sliding-25 nt max**, and the motif
  threshold is now calibrated per window size against random mm10 at a stated
  false-positive rate (~2.28 at 501 bp, ~2.60 at 1001 bp for 5 % FPR) instead of
  the tautological `g4_threshold: 0.0`, which made `has_g4_motif` true for every
  row.
- **Region fold-enrichment null corrected for peak width and un-merged
  regions.** Expected overlap is now the `reduce()`d region footprint dilated by
  the peak's own width. Largely removes the inverse-region-size artefact that
  caused the fold-enrichment figure to be dropped in `a4acb9f`
  (5′ UTR ~11.9 → ~3.2; promoter 1.99 → 1.62).
- **`peak_calling.z_per_replicate: true`** — each replicate is z-scored
  independently and the threshold must be met in ≥ `min_replicates` of them,
  which is what the README always claimed. Set `false` to reproduce v1.0
  exactly. `min_width` 500 → 1000, since 500 == `bin_size` made it inert.

### Added

- **`26d_heuristic_control_benchmark.R`** — scores the loop-length heuristic on
  the same 10 literature control G4s as G4SP. Result: **G4SP 8/10, heuristic
  2/10**, below the 33 % expected from guessing at three classes. This is the
  evidence that justifies removing the back-fill.
- **`21b_topology_sensitivity.R`** and rule `tss_metaprofiles_argmax` — the
  confidence-stratified sensitivity arm. Recomputes topology at plain argmax
  (no abstention, 11,062 hybrids retained) from the per-class probabilities
  already in `g4sp_topology.csv`, and renders the TSS metaprofile under both
  arms, so the ordering is reported against the precision threshold rather than
  at one setting of it.
- **`verify/`** — five standalone checks, no bigWigs or `cache/` needed, ~1 min.
- `topology_concordance_stats.csv` — κ (−0.0098) plus the marginal-free
  cross-classification rates (1.8 % / 0.7 %), which are the fairer summary given
  the divergent marginals. Previously κ reached no committed artefact at all.
- `topology_source_composition.csv` + figure — the `topology_source` column
  existed in v1.0 but was never aggregated anywhere.
- `strand_g4_threshold_calibration.tsv`, `strand_score_tests.tsv` (Wilcoxon) and
  `strand_detection_tests.tsv` (chi-square). The latter two were advertised in
  the v1.0 script header but never implemented.
- `background_gc_matching.csv` + `cache/background_seqs_unmatched.fa` — the
  unmatched-background arm and a realised-vs-target GC diagnostic for a matching
  step that truncated silently.
- A **Limitations** section in the README, and the null promoter-G4-expression
  result promoted into "Selected results".

### Fixed

- `.bw_rle_cache` was keyed on absolute path only, ignoring the requested
  seqlevels, so the first call fixed the chromosome set for every later call on
  that file. Now keyed on path + seqlevels, with `evict_bigwigs()` updated to
  match — without which eviction would silently free nothing, and those calls
  are load-bearing for memory.
- `07_call_peaks.R` never forwarded `z_per_replicate` to
  `call_peaks_from_bigwigs()`, and `37_unclassified_peak_motifs.R` fell back to
  the literal `16980` when `unclassified.subsample` was null — both would have
  made their config keys silently inert.
- The new `mixed` stratum is now carried through every consumer of
  `peak_topology.csv`. Nine scripts would otherwise have dropped it: three
  metaprofile scripts subset to a hard-coded `c(definite, "none")`, and three
  more built factor levels without it, which turns those rows into `NA` and
  renormalises the remaining percentages over a smaller denominator.

### Documentation

- Corrected nine README claims: the four "permutation" lines, "replicate-aware",
  "peaks present in ≥ 2 replicates", "peak-summit sequences", "11
  literature-validated controls" (there are 10), "genuine reproduction of the
  calling step", and the unsupported ERCC "antibody-specific" sentence — the
  ERCC tracks use the same BG4 antibody and differ by batch, and the R-loop side
  inverts (ERCC 1.31–1.37 vs real S9.6 WT 0.99).

### Not done

No end-to-end re-run yet. Also deferred: a real `regioneR::permTest` with
`circularRandomizeRegions` (Linux only; hangs on Windows), an accessibility
covariate for the G4 × R-loop correlation (no ATAC track in the manifest), and
an IgG null for the peak caller (no IgG track exists in the GEO submission).

## [1.0.0] — 2026-07-30

Initial public release.

- One Snakemake workflow (~50 rules, 50 numbered R/Python scripts) covering the
  full analysis from public CUT&Tag bigWigs plus an RNA-seq count table through
  to every figure and table, in twelve stages (A–L).
- **Stage A** — RNA-seq QC, DESeq2 across five contrasts with ashr shrinkage,
  GO/GSEA enrichment, gene-biotype resolved volcano plots.
- **Stage B** — replicate-aware z-score peak calling for G4 (BG4) and R-loop
  (S9.6) CUT&Tag, and construction of promoter, 5′ UTR, first-intron, gene-body
  and ENCODE cCRE enhancer region sets.
- **Stages C–E** — genome-wide CUT&Tag QC, regional enrichment and metaprofiles,
  and integration of differential expression with promoter/TSS occupancy.
- **Stage F** — sequence-level G-quadruplex topology from pqsfinder, DeepG4 and
  G4ShapePredictor, calibrated against literature-validated control quadruplexes,
  plus gain/loss, propensity and feature-partition analyses.
- **Stages G–L** — topology-resolved metaprofiles, characterisation of peaks
  lacking a canonical quadruplex motif, lncRNA G4 profiles, G4 × R-loop
  correlation and co-occurrence, all-PQS/multi-G4 promoter analysis at motif
  resolution, strand asymmetry, and Gviz locus tracks.
- Repository-local data layout with `scripts/download_data.ps1` / `.sh` to fetch
  the 28 CUT&Tag bigWigs and the RNA-seq count table from NCBI GEO.
- `envs/install_r_packages.R` to install the required CRAN, Bioconductor and
  GitHub R packages; conda environment files for the Python and Linux/WSL routes.
- MIT licence for the code, CC-BY-4.0 for the derived tables and figures,
  `CITATION.cff` and `.zenodo.json` for citation and archiving.
