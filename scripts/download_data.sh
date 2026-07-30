#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Fetch the public input data for G4_Promoter_Pipeline from NCBI GEO.
#
# Downloads into the repository-local layout expected by config/config.yaml:
#   data/bigwig/   28 CUT&Tag coverage tracks (GSE269084; ~3.4 GB total)
#   data/rnaseq/   GSE269081_count_table.tsv  (RNA-seq counts)
# and clones the G4ShapePredictor topology model into external/.
#
# The file list, expected sizes and sample metadata come from
# data/geo_manifest.tsv. Complete files are skipped and partial transfers are
# resumed, so the script is safe to re-run after an interruption.
#
# All 28 bigWigs are needed for a complete run: the ERCC tracks feed the ERCC
# peak sets built by scripts/07_call_peaks.R, which `rule all` requires. Use
# --skip-ercc only if you also drop ERCCWT/ERCCKO from peak_calling.genotypes
# in config/config.yaml.
#
# Usage:
#   scripts/download_data.sh [--refs] [--skip-ercc] [--dry-run]
#     --refs       also fetch GENCODE vM25 GTF (~28 MB) and mm10 RepeatMasker
#                  (~142 MB) into data/ (otherwise the workflow fetches them
#                  itself on first run)
#     --skip-ercc  skip the 8 ERCC control tracks (~1.3 GB)
#     --dry-run    print the plan without downloading
#
# Already have the bigWigs elsewhere? Don't copy them — point the workflow at
# them with a local overlay instead:
#   cp config/config.local.example.yaml config/config.local.yaml   # edit paths
#   snakemake --cores 4 --configfile config/config.local.yaml
# ---------------------------------------------------------------------------
set -euo pipefail

REFS=0; SKIP_ERCC=0; DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --refs)      REFS=1 ;;
    --skip-ercc) SKIP_ERCC=1 ;;
    --dry-run|-n) DRY_RUN=1 ;;
    -h|--help)   sed -n '2,32p' "$0"; exit 0 ;;
    *) echo "Unknown option: $arg (try --help)" >&2; exit 2 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$REPO_ROOT/data/geo_manifest.tsv"
BIGWIG_DIR="$REPO_ROOT/data/bigwig"
RNASEQ_DIR="$REPO_ROOT/data/rnaseq"

GEO_SAMPLES="https://ftp.ncbi.nlm.nih.gov/geo/samples"
COUNTS_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269081/suppl/GSE269081_count_table.tsv.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
RMSK_URL="https://hgdownload.gi.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"

command -v curl >/dev/null || { echo "curl is required" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

file_size() { [[ -f "$1" ]] && stat -c%s "$1" 2>/dev/null || stat -f%z "$1" 2>/dev/null || echo 0; }

failed=()

fetch() {
  # fetch <url> <dest> [expected_bytes]
  local url="$1" dest="$2" expected="${3:-0}" name
  name="$(basename "$dest")"
  if [[ "$expected" -gt 0 && -f "$dest" && "$(file_size "$dest")" -eq "$expected" ]]; then
    printf '  [skip] %s (already complete)\n' "$name"; return 0
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '  [plan] %s\n         <- %s\n' "$name" "$url"; return 0
  fi
  printf '  [get ] %s\n' "$name"
  if ! curl --fail --location --retry 3 --retry-delay 5 --continue-at - \
            --progress-bar --output "$dest" "$url"; then
    echo "  download failed: $url" >&2; failed+=("$name"); return 1
  fi
  if [[ "$expected" -gt 0 ]]; then
    local actual; actual="$(file_size "$dest")"
    if [[ "$actual" -ne "$expected" ]]; then
      echo "  size mismatch for $name: got $actual bytes, expected $expected" >&2
      failed+=("$name"); return 1
    fi
  fi
}

# GEO buckets samples by GSM accession with the last three digits masked:
#   GSM8305989 -> .../samples/GSM8305nnn/GSM8305989/suppl/<filename>
gsm_url() { printf '%s/%snnn/%s/suppl/%s' "$GEO_SAMPLES" "${1:0:$((${#1}-3))}" "$1" "$2"; }

mkdir -p "$BIGWIG_DIR" "$RNASEQ_DIR"

# --- CUT&Tag bigWigs --------------------------------------------------------
mapfile -t rows < <(awk -F'\t' -v skip="$SKIP_ERCC" \
  'NR>1 && !(skip==1 && $8=="ercc") {print $1"\t"$2"\t"$3}' "$MANIFEST")
total_bytes=0
for row in "${rows[@]}"; do total_bytes=$(( total_bytes + $(cut -f3 <<<"$row") )); done
printf '\nCUT&Tag bigWigs: %d files, %.2f GB  ->  data/bigwig/\n' \
  "${#rows[@]}" "$(awk -v b="$total_bytes" 'BEGIN{printf "%.2f", b/1073741824}')"
for row in "${rows[@]}"; do
  IFS=$'\t' read -r gsm fname bytes <<<"$row"
  fetch "$(gsm_url "$gsm" "$fname")" "$BIGWIG_DIR/$fname" "$bytes" || true
done

# --- RNA-seq count table ----------------------------------------------------
printf '\nRNA-seq count table (GSE269081)  ->  data/rnaseq/\n'
if [[ -f "$RNASEQ_DIR/GSE269081_count_table.tsv" ]]; then
  echo "  [skip] GSE269081_count_table.tsv (already present)"
else
  if fetch "$COUNTS_URL" "$RNASEQ_DIR/GSE269081_count_table.tsv.gz"; then
    [[ "$DRY_RUN" -eq 1 ]] || { gunzip -f "$RNASEQ_DIR/GSE269081_count_table.tsv.gz"
                                echo "  [ok  ] unpacked GSE269081_count_table.tsv"; }
  fi
fi

# --- Optional reference annotations ----------------------------------------
if [[ "$REFS" -eq 1 ]]; then
  printf '\nReference annotations  ->  data/\n'
  fetch "$GTF_URL"  "$REPO_ROOT/data/gencode.vM25.annotation.gtf.gz" || true
  fetch "$RMSK_URL" "$REPO_ROOT/data/rmsk_mm10.txt.gz"               || true
fi

# --- Upstream model repository ---------------------------------------------
printf '\nG4ShapePredictor (topology model)  ->  external/\n'
if [[ -d "$REPO_ROOT/external/G4ShapePredictor/.git" ]]; then
  echo "  [skip] external/G4ShapePredictor (already cloned)"
elif [[ "$DRY_RUN" -eq 1 ]]; then
  echo "  [plan] git clone --depth 1 https://github.com/donn-liew/G4ShapePredictor"
elif command -v git >/dev/null; then
  mkdir -p "$REPO_ROOT/external"
  git clone --depth 1 https://github.com/donn-liew/G4ShapePredictor \
    "$REPO_ROOT/external/G4ShapePredictor"
else
  echo "  git not found - clone https://github.com/donn-liew/G4ShapePredictor into external/ manually" >&2
fi

printf '\n'
if (( ${#failed[@]} > 0 )); then
  printf 'WARNING: %d item(s) did not complete: %s\n' "${#failed[@]}" "${failed[*]}" >&2
  echo "Re-run this script to resume; downloads continue where they stopped." >&2
  exit 1
fi
echo "Done. Next: conda activate G4 && snakemake -n"
