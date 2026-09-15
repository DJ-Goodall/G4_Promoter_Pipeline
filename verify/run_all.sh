#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Run every v1.1 verification check. None of these needs cache/ or the bigWigs;
# they read the committed v1.0 tables, the control TSV, and BSgenome mm10.
#
# From the repository root:  bash verify/run_all.sh
# ---------------------------------------------------------------------------
set -u
RSCRIPT="${RSCRIPT:-Rscript}"
fail=0

run () {
  echo ""
  echo "=============================================================="
  echo "  $1"
  echo "=============================================================="
  # Do NOT wrap this in `if ...; then :; fi` -- the `:` runs after the pipeline
  # and resets PIPESTATUS, which silently made every check report success.
  "$RSCRIPT" "$2" 2>&1 | grep -v "built under R version"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]; then echo "!! FAILED: $2 (exit $rc)"; fail=1; fi
}

echo "### Static syntax check ###"
for f in scripts/*.R scripts/_shared/*.R verify/*.R; do
  if "$RSCRIPT" -e "invisible(parse('$f'))" >/dev/null 2>&1; then
    echo "OK   $f"
  else
    echo "FAIL $f"; fail=1
  fi
done

# Check 1 asserts the v1.0 numbers REVIEW.md was written against, so it must read
# v1.0 tables. Once v1.1 has been run, results/tables/ holds the new ones and every
# assertion mismatches by design -- so materialise the v1.0 tables from git instead.
V10_REF="${V10_REF:-3c1b264}"
V10_TREE="$(mktemp -d)"
# git archive is used rather than `git checkout -- <path>`: it writes only to the
# temp dir and never touches the index or the working tree.
if git rev-parse --verify --quiet "$V10_REF" >/dev/null &&
   git archive "$V10_REF" results/tables | tar -x -C "$V10_TREE"; then
  export V11_TABLES="$V10_TREE/results/tables"
  echo ""
  echo "Check 1 reads the v1.0 tables from $V10_REF (V11_TABLES=$V11_TABLES)"
else
  echo ""
  echo "!! v1.0 ref '$V10_REF' not found -- check 1 will read the working tree and"
  echo "   will fail if v1.1 has been run. Set V10_REF to the v1.0 commit."
fi

run "1. Re-derive every REVIEW.md number from the committed v1.0 tables" \
    verify/rederive_audit_numbers.R

unset V11_TABLES
rm -rf "$V10_TREE"

run "2. G4Hunter: window mean vs sliding max, and why 1.2 works for neither" \
    verify/test_g4hunter.R

run "3. Width-aware enrichment null vs brute-force simulation" \
    verify/test_enrichment_null.R

run "4. Loop heuristic vs G4ShapePredictor on the literature controls" \
    verify/run_heuristic_benchmark.R

run "5. Calibrate the G4Hunter motif threshold against real mm10" \
    verify/calibrate_g4_threshold.R

run "6. Scaled-region profile geometry (rules 52/53)" \
    verify/test_scaled_profile.R

echo ""
if [ "$fail" -eq 0 ]; then
  echo "=============================================================="
  echo "  ALL v1.1 VERIFICATION CHECKS PASSED"
  echo "=============================================================="
else
  echo "!! Some checks failed -- see above."
  exit 1
fi
