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
  if "$RSCRIPT" "$2" 2>&1 | grep -v "built under R version"; then :; fi
  # shellcheck disable=SC2181
  if [ "${PIPESTATUS[0]}" -ne 0 ]; then echo "!! FAILED: $2"; fail=1; fi
}

echo "### Static syntax check ###"
for f in scripts/*.R scripts/_shared/*.R verify/*.R; do
  if "$RSCRIPT" -e "invisible(parse('$f'))" >/dev/null 2>&1; then
    echo "OK   $f"
  else
    echo "FAIL $f"; fail=1
  fi
done

run "1. Re-derive every REVIEW.md number from the committed v1.0 tables" \
    verify/rederive_audit_numbers.R

run "2. G4Hunter: window mean vs sliding max, and why 1.2 works for neither" \
    verify/test_g4hunter.R

run "3. Width-aware enrichment null vs brute-force simulation" \
    verify/test_enrichment_null.R

run "4. Loop heuristic vs G4ShapePredictor on the literature controls" \
    verify/run_heuristic_benchmark.R

run "5. Calibrate the G4Hunter motif threshold against real mm10" \
    verify/calibrate_g4_threshold.R

echo ""
if [ "$fail" -eq 0 ]; then
  echo "=============================================================="
  echo "  ALL v1.1 VERIFICATION CHECKS PASSED"
  echo "=============================================================="
else
  echo "!! Some checks failed -- see above."
  exit 1
fi
