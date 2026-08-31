# `verify/` — standalone checks for the v1.1 audit corrections

These are not pipeline rules. They are independent checks that the corrections
made in response to the [2026-08-07 audit](../docs/audit-2026-08/) are
justified by data rather than asserted — every quantitative claim in
[`REVIEW.md`](../docs/audit-2026-08/REVIEW.md) is re-derived here.

They need **no bigWigs and no `cache/`**. They read the committed summary
tables under `results/tables/`, the control TSV under `data/`, and BSgenome
mm10. The whole suite takes about a minute.

```bash
# from the repository root
bash verify/run_all.sh

# or point at a specific R, e.g. on Windows
RSCRIPT="/c/Program Files/R/R-4.5.1/bin/Rscript.exe" bash verify/run_all.sh
```

Run this with the **system R**, not inside `conda activate G4` — the conda
environment shadows R 4.5.1 with a build that has no Bioconductor packages, and
checks 3–5 die with `no package called 'GenomicRanges'`.

Check 1 asserts the **v1.0** numbers `REVIEW.md` was written against, so it must
read v1.0 tables. Once the v1.1 pipeline has been run, `results/tables/` no
longer holds them, so `run_all.sh` materialises them from `$V10_REF` (default
`3c1b264`) with `git archive` into a temp dir and passes the path in
`$V11_TABLES`. Nothing is written to the index or the working tree.

| Check | What it establishes |
|---|---|
| `rederive_audit_numbers.R` | Every number quoted in `REVIEW.md`, `stopifnot()`-ed against the committed v1.0 tables (read from `$V11_TABLES`) |
| `test_g4hunter.R` | The v1.0 window-mean statistic calls **0/10** literature controls at the conventional 1.2 cutoff; the v1.1 sliding-25 nt max calls 9/10 but at a 36% false-positive rate — which is why the threshold is calibrated rather than hard-coded |
| `test_enrichment_null.R` | The analytic width-aware expectation used by `12_peak_region_enrichment.R` matches brute-force simulation within 1% at every peak width from 1 to 2000 bp |
| `run_heuristic_benchmark.R` | **G4ShapePredictor 8/10 vs the loop heuristic 2/10** on the same 10 literature controls — the evidence for removing the back-fill in rules 21 and 43 |
| `calibrate_g4_threshold.R` | The G4Hunter threshold calibration against real mm10 (~2.28 at 501 bp, ~2.60 at 1001 bp for a 5% FPR) |

`run_all.sh` also parse-checks every script under `scripts/`, so a syntax error
anywhere in the workflow fails the suite.

Two CSVs are committed alongside — `topology_control_benchmark_compare.csv` and
`g4_threshold_calibration.csv` — so the results can be read without re-running.
Both are regenerated in place by the suite.

## A note on what these do and do not test

They test the *reasoning*, not the pipeline. They confirm that the corrections
are warranted and that the replacement statistics behave as claimed on external
standards. They do not run the workflow — for that, see the smoke-test route in
the main [README](../README.md#reproducing).
