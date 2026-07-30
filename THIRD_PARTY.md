# Third-party software and reference data

This workflow orchestrates existing tools and public reference data; none of it is
redistributed in this repository. Everything listed below is fetched by the setup
scripts or by the workflow itself at run time, and each item keeps its own
licence and citation. Please cite the underlying tools alongside this workflow.

## Sequence models (cloned or installed separately)

| Tool | Source | Licence | Role here |
|---|---|---|---|
| **DeepG4** | [raphaelmourad/DeepG4](https://github.com/raphaelmourad/DeepG4) — installed as an R package with `remotes::install_github()` | see upstream repository | G-quadruplex formation-probability CNN; independent sanity check on called peaks and the ROC in script `22_sanity_check.R` |
| **G4ShapePredictor** | [donn-liew/G4ShapePredictor](https://github.com/donn-liew/G4ShapePredictor) — cloned into `external/` | MIT (© 2023 Donn Liew) | Topology classification (parallel / antiparallel / hybrid) of putative quadruplex sequences in scripts `20`, `26a`, `26b`, and at motif resolution in Stage K |

The shipped G4ShapePredictor model pickles were trained with scikit-learn 1.0.2.
The `envs/g4sp.yaml` pins (python 3.9, scikit-learn 1.0.2, numpy < 2) exist because
newer scikit-learn versions cannot unpickle those models — see the comments in that
file before changing them.

## R / Bioconductor packages

Installed by `envs/install_r_packages.R`. The scientifically load-bearing ones:

| Package | Role | Citation |
|---|---|---|
| **pqsfinder** | Putative quadruplex sequence detection and loop-length topology heuristic | Hon J, Martínek T, Zendulka J, Lexa M. *Bioinformatics* (2017) |
| **DESeq2** | Differential expression across the five genotype contrasts | Love MI, Huber W, Anders S. *Genome Biology* (2014) |
| **limma / edgeR** | voom-based differential peak signal (gain/loss, G4 × R-loop) | Ritchie ME et al. *NAR* (2015); Robinson MD et al. *Bioinformatics* (2010) |
| **clusterProfiler / fgsea / msigdbr / org.Mm.eg.db** | GO and GSEA functional enrichment | Wu T et al. *The Innovation* (2021); Korotkevich G et al. (2021) |
| **rtracklayer / GenomicRanges / GenomicFeatures / Biostrings** | bigWig I/O, interval arithmetic, annotation, sequence handling | Lawrence M et al. *PLoS Comp Biol* (2013) |
| **BSgenome.Mmusculus.UCSC.mm10 / TxDb.Mmusculus.UCSC.mm10.knownGene** | mm10 genome sequence and transcript models | Bioconductor annotation packages |
| **motifmatchr / TFBSTools / JASPAR2020** | Transcription-factor motif enrichment in peaks without a canonical PQS | Tan G, Lenhard B. *Bioinformatics* (2016); Fornes O et al. *NAR* (2020) |
| **Gviz** | Locus-level track figures | Hahne F, Ivanek R. *Methods Mol Biol* (2016) |
| **EnhancedVolcano, UpSetR, pheatmap, ggplot2, ggpubr, cowplot, ggrepel, pROC** | Figures and summary statistics | see individual packages |

Bioconductor packages are released under Artistic-2.0 or similar open licences;
`sessionInfo()` from any run records the exact versions used.

## Reference data (downloaded at run time)

| Resource | Source | Terms | Used by |
|---|---|---|---|
| **GENCODE mouse release M25 annotation** | `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/` | free for academic and commercial use with attribution | gene biotypes (`05`), lncRNA definitions (`38`) |
| **ENCODE SCREEN V3 mm10 cCREs** | `downloads.wenglab.org/Registry-V3/GRCm38-cCREs.bed` (ENCODE fallback `ENCFF846BKF`) | ENCODE data release policy — cite ENCODE | enhancer region set (`09`) |
| **UCSC mm10 RepeatMasker table** | `hgdownload.gi.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz` | UCSC Genome Browser conditions of use | repeat context of unclassified peaks (`36`) |
| **UCSC mm10 cytoband / ideogram** | UCSC via Gviz | UCSC Genome Browser conditions of use | ideogram track (`49`; skipped gracefully when offline) |

## Input datasets

The CUT&Tag and RNA-seq data analysed here are public in NCBI GEO under
**GSE269081** (RNA-seq), **GSE269082** (G4 and R-loop CUT&Tag) and **GSE269084**
(series), from *Sato et al.*, "RNA transcripts regulate G-quadruplex landscapes
through G-loop formation". Please cite the original study when using these data.

## Hand-curated content in this repository

`data/known_topology_controls.tsv` lists eleven experimentally determined G-quadruplex
topologies (NMR/crystallography, K⁺ conditions) compiled from the primary
literature — c-MYC Pu27, c-KIT1/2, VEGF, KRAS, HRAS, BCL2 Pu39, human telomeric
repeat, the thrombin-binding aptamer and the *Oxytricha* telomeric G4. It is used
to benchmark and calibrate the topology classifier (`26b`, `26c`); each row cites
its structural basis in the `note` column.
