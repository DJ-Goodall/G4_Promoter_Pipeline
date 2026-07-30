#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Install every R package the workflow needs.
#
#   Rscript envs/install_r_packages.R              # install what is missing
#   Rscript envs/install_r_packages.R --check      # report only, install nothing
#
# Developed against R 4.5.1 with Bioconductor 3.21. The package list below is
# the set actually referenced by scripts/ (library() / requireNamespace() /
# pkg::fun calls), not a superset.
#
# On native Windows this installs into your system R library — bioconda has no
# Windows builds of the Bioconductor genomics stack, so conda cannot provide it
# (see README, "Installation"). On Linux/macOS you may prefer the conda route:
#   conda env create -f envs/r_g4.yaml
# ---------------------------------------------------------------------------

args      <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

cran <- c(
  # data wrangling / IO
  "dplyr", "tidyr", "readr", "stringr", "yaml", "jsonlite", "R.utils", "remotes",
  # plotting
  "ggplot2", "ggpubr", "ggrepel", "cowplot", "pheatmap", "RColorBrewer", "scales",
  "UpSetR", "hexbin",
  # stats / modelling
  "pROC", "ashr", "msigdbr",
  # python bridge for the DeepG4 rule
  "reticulate"
)

bioc <- c(
  # core infrastructure
  "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges",
  "SummarizedExperiment", "AnnotationDbi",
  # sequence / annotation
  "Biostrings", "BSgenome", "rtracklayer", "GenomicFeatures", "pqsfinder",
  "BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db",
  # differential expression + enrichment
  "DESeq2", "edgeR", "limma", "EnhancedVolcano", "clusterProfiler", "fgsea",
  # motifs
  "TFBSTools", "motifmatchr", "JASPAR2020",
  # tracks
  "Gviz"
)

github <- c(DeepG4 = "raphaelmourad/DeepG4")

installed <- rownames(installed.packages())
missing_cran   <- setdiff(cran, installed)
missing_bioc   <- setdiff(bioc, installed)
missing_github <- names(github)[!names(github) %in% installed]

msg <- function(...) cat(..., "\n", sep = "")
msg("R ", getRversion(), " | library: ", .libPaths()[1])
msg("CRAN:        ", length(cran) - length(missing_cran), "/", length(cran), " present")
msg("Bioconductor: ", length(bioc) - length(missing_bioc), "/", length(bioc), " present")
msg("GitHub:      ", length(github) - length(missing_github), "/", length(github), " present")

if (length(missing_cran) + length(missing_bioc) + length(missing_github) == 0) {
  msg("\nAll packages present — nothing to do.")
  quit(status = 0)
}
if (length(missing_cran))   msg("\nMissing from CRAN:        ", paste(missing_cran, collapse = ", "))
if (length(missing_bioc))   msg("Missing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
if (length(missing_github)) msg("Missing from GitHub:      ", paste(missing_github, collapse = ", "))

if (check_only) {
  msg("\n--check given: no packages installed.")
  quit(status = 1)
}

if (length(missing_cran)) {
  msg("\nInstalling CRAN packages...")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  msg("\nInstalling Bioconductor packages (this pulls large annotation packages)...")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

if (length(missing_github)) {
  # DeepG4 is not on CRAN or Bioconductor. It needs TensorFlow/Keras through
  # reticulate at run time — see config/config.yaml: deepg4.reticulate_condaenv.
  # If this step or TensorFlow fails, only the ROC sanity check (script 22) is
  # lost; the topology branch still completes via pqsfinder + G4ShapePredictor.
  msg("\nInstalling GitHub packages...")
  for (pkg in missing_github) {
    msg("  ", pkg, " <- ", github[[pkg]])
    tryCatch(remotes::install_github(github[[pkg]], upgrade = "never"),
             error = function(e) warning("failed to install ", pkg, ": ",
                                         conditionMessage(e), call. = FALSE))
  }
}

still_missing <- setdiff(c(cran, bioc, names(github)), rownames(installed.packages()))
if (length(still_missing)) {
  msg("\nStill missing after installation: ", paste(still_missing, collapse = ", "))
  quit(status = 1)
}
msg("\nAll packages installed. Next: conda activate G4; snakemake -n")
