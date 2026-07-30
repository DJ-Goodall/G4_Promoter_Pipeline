# Changelog

All notable changes to this workflow are recorded here. Versions follow
[semantic versioning](https://semver.org/): the major version changes when the
scientific output changes, the minor version when rules or parameters are added,
and the patch version for documentation and fixes that leave results unchanged.

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
