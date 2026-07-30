# Licence for the data products in this repository

The MIT terms in [`LICENSE`](LICENSE) cover the **workflow code**: the `Snakefile`,
`config/`, `scripts/`, `envs/` and the documentation.

The **derived data products** distributed alongside the code are released under the
**Creative Commons Attribution 4.0 International** licence
([CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)):

- `results/tables/*.csv`, `results/tables/*.tsv` — summary tables from the v1.0.0 run
- `docs/figures/*.png` — the showcase figures reproduced in the README
- `data/known_topology_controls.tsv` — the hand-curated set of literature-validated
  G-quadruplex topologies (the underlying structural determinations belong to the
  primary publications cited in its `note` column)

You may share and adapt these with attribution to this repository and to the source
study (*Sato et al.*, GEO GSE269081 / GSE269082 / GSE269084).

**Not covered here.** Third-party software invoked by this workflow (DeepG4,
G4ShapePredictor, pqsfinder, the Bioconductor stack) and third-party reference data
(GENCODE, ENCODE SCREEN, UCSC RepeatMasker, JASPAR) are neither included nor
redistributed in this repository and carry their own licences — see
[`THIRD_PARTY.md`](THIRD_PARTY.md). The input CUT&Tag and RNA-seq data remain under
the terms of their GEO deposition.
