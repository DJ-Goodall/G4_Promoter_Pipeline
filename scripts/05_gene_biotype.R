#!/usr/bin/env Rscript
# ============================================================================
# 21_gene_biotype.R   (env: system R 4.5.1)   --- Phase 9 ---
#
# Attach a gene biotype to every gene so the RNA-seq volcano (rule 22) and the
# lncRNA G4 metaprofiles (rule 23) can separate lncRNAs from coding genes.
#
# The count table (and therefore deseq2_KO_vs_WT.csv) carries ENSEMBL gene_ids
# (ENSMUSG, unversioned). GENCODE vM25 is the final GRCm38/mm10 release; it is
# Ensembl-based (same ENSMUSG ids), chr-named like UCSC mm10 (no seqlevelsStyle
# headache), and uses a single consolidated gene_type "lncRNA". We download the
# GTF once, parse the gene-level records, strip the GENCODE version suffix so the
# ids join the count table, and cache:
#   - results/tables/gene_biotype.csv  (gene_id, gene_name, gene_type, coords)
#   - cache/gencode_genes.rds          (all genes GRanges, names = gene_id)
#   - cache/lncrna_genes.rds           (lncRNA genes only, UCSC std chroms)
#
# Inputs:   (network) GENCODE GTF -> data/gencode.vM25.annotation.gtf.gz
# Deps:     rtracklayer, GenomicRanges, GenomeInfoDb, S4Vectors, dplyr, readr
# ============================================================================

suppressPackageStartupMessages({
  library(rtracklayer); library(GenomicRanges); library(GenomeInfoDb)
  library(S4Vectors); library(dplyr); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("data", "cache", "results/tables")

std_chroms   <- cfg$std_chroms
# Mouse GENCODE vM25 uses legacy long-noncoding biotypes (no single "lncRNA"); the
# lncRNA set is their union. Fall back to scalar `biotype` for older configs.
lnc_biotypes <- as.character(unlist(cfg$lncrna$biotypes %||% cfg$lncrna$biotype %||% "lncRNA"))
gtf_gz       <- cfg$lncrna$gtf_local %||% "data/gencode.vM25.annotation.gtf.gz"

# --- Download GENCODE GTF once (multi-mirror fallback, as in rule 19) --------
if (!file.exists(gtf_gz) || file.info(gtf_gz)$size < 1e6) {
  urls <- unique(c(
    cfg$lncrna$gtf_url,
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
    "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz"))
  urls <- urls[!is.na(urls) & nzchar(urls)]
  old_to <- getOption("timeout"); options(timeout = 1200); on.exit(options(timeout = old_to), add = TRUE)
  ok <- FALSE
  for (url in urls) {
    message("Downloading GENCODE GTF: ", url)
    ok <- tryCatch({
      utils::download.file(url, gtf_gz, mode = "wb", quiet = TRUE)
      file.exists(gtf_gz) && file.info(gtf_gz)$size > 1e6
    }, error = function(e) { message("  failed: ", conditionMessage(e)); FALSE })
    if (ok) break
  }
  if (!ok) stop("Could not download GENCODE vM25 GTF from any mirror; download it manually to ", gtf_gz)
}

# --- Parse gene-level records ----------------------------------------------
message("Importing GTF (gene records only) ...")
gtf <- rtracklayer::import(gtf_gz, feature.type = "gene")
# GENCODE gene_ids are versioned (ENSMUSG00000051951.5) -> strip to join counts.
S4Vectors::mcols(gtf)$gene_id <- sub("\\..*$", "", as.character(S4Vectors::mcols(gtf)$gene_id))
names(gtf) <- S4Vectors::mcols(gtf)$gene_id

biotype_df <- data.frame(
  gene_id   = as.character(S4Vectors::mcols(gtf)$gene_id),
  gene_name = as.character(S4Vectors::mcols(gtf)$gene_name),
  gene_type = as.character(S4Vectors::mcols(gtf)$gene_type),
  chr       = as.character(GenomicRanges::seqnames(gtf)),
  start     = GenomicRanges::start(gtf),
  end       = GenomicRanges::end(gtf),
  strand    = as.character(GenomicRanges::strand(gtf)),
  stringsAsFactors = FALSE) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)
readr::write_csv(biotype_df, "results/tables/gene_biotype.csv")
message("Wrote gene_biotype.csv: ", nrow(biotype_df), " genes; biotype breakdown (top):")
print(utils::head(sort(table(biotype_df$gene_type), decreasing = TRUE), 8))

# --- Cache full gene GRanges + lncRNA subset (UCSC std chroms) ---------------
gtf_std <- GenomeInfoDb::keepSeqlevels(
  gtf, intersect(GenomeInfoDb::seqlevels(gtf), std_chroms), pruning.mode = "coarse")
saveRDS(gtf_std, "cache/gencode_genes.rds")

lncrna <- gtf_std[as.character(S4Vectors::mcols(gtf_std)$gene_type) %in% lnc_biotypes]
names(lncrna) <- as.character(S4Vectors::mcols(lncrna)$gene_id)
saveRDS(lncrna, "cache/lncrna_genes.rds")

message(sprintf("Done. %d genes total, %d lncRNA genes (biotypes: %s) on std chroms.",
                length(gtf_std), length(lncrna), paste(lnc_biotypes, collapse = ", ")))
