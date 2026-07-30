#!/usr/bin/env Rscript
# ============================================================================
# 08_build_regions.R   (system R)   --- Stage B ---
#
# Build TxDb-derived regulatory-region GRanges (promoter, 5'UTR, first intron,
# gene body) for mm10. Ported from 20260522_g4_gloop_extended_V3.Rmd chunk
# s2-txdb-regions (+ s2-save-regions individual files). The enhancer region and
# the combined regions_all.rds are built by 09_build_enhancers.R.
#
# Outputs (into cache/):
#   regions_txdb.rds     list(promoter, 5UTR, intron1, genebody)
#   regions_promoter.rds regions_5UTR.rds regions_intron1.rds regions_genebody.rds
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache")

std_chroms <- cfg$std_chroms
up   <- cfg$regions$promoter_up   %||% 2000
down <- cfg$regions$promoter_down %||% 2000

txdb_regions <- cache_or_build(file.path("cache", "regions_txdb.rds"), {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"

  genes_gr <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
  genes_gr <- GenomeInfoDb::keepSeqlevels(genes_gr, std_chroms, pruning.mode = "coarse")

  prom <- GenomicFeatures::promoters(genes_gr, upstream = up, downstream = down)
  names(prom) <- names(genes_gr)

  # 5'UTR by transcript -> per-gene union
  utr5_by_tx <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = FALSE)
  tx2gene    <- AnnotationDbi::select(txdb, keys = names(utr5_by_tx),
                                      columns = c("TXID", "GENEID"), keytype = "TXID")
  utr5_flat  <- unlist(utr5_by_tx, use.names = FALSE)
  utr5_flat$gene_id <- rep(tx2gene$GENEID[match(names(utr5_by_tx), tx2gene$TXID)],
                           lengths(utr5_by_tx))
  utr5_flat <- utr5_flat[!is.na(utr5_flat$gene_id)]
  utr5_flat <- GenomeInfoDb::keepSeqlevels(utr5_flat, std_chroms, pruning.mode = "coarse")

  # First intron per transcript (strand-aware): + strand smallest start, - largest end
  introns_by_tx <- GenomicFeatures::intronsByTranscript(txdb, use.names = FALSE)
  nonempty <- lengths(introns_by_tx) > 0
  first_intron_list <- lapply(which(nonempty), function(i) {
    x <- introns_by_tx[[i]]
    sel <- if (as.character(GenomicRanges::strand(x))[1] == "-")
      x[which.max(GenomicRanges::end(x))] else x[which.min(GenomicRanges::start(x))]
    sel$gene_id <- tx2gene$GENEID[match(names(introns_by_tx)[i], tx2gene$TXID)]
    sel
  })
  first_intron <- unlist(GenomicRanges::GRangesList(first_intron_list), use.names = FALSE)
  first_intron <- first_intron[!is.na(first_intron$gene_id)]
  first_intron <- GenomeInfoDb::keepSeqlevels(first_intron,
    intersect(GenomeInfoDb::seqlevels(first_intron), std_chroms), pruning.mode = "coarse")

  collapse_per_gene <- function(gr) {
    spl <- GenomicRanges::split(gr, gr$gene_id)
    red <- GenomicRanges::reduce(spl, ignore.strand = FALSE)
    out <- unlist(red, use.names = FALSE)
    out$gene_id <- rep(names(red), lengths(red))
    out
  }
  list(promoter = prom,
       `5UTR`    = collapse_per_gene(utr5_flat),
       intron1  = collapse_per_gene(first_intron),
       genebody = genes_gr)
})

saveRDS(txdb_regions$promoter, file.path("cache", "regions_promoter.rds"))
saveRDS(txdb_regions$`5UTR`,   file.path("cache", "regions_5UTR.rds"))
saveRDS(txdb_regions$intron1,  file.path("cache", "regions_intron1.rds"))
saveRDS(txdb_regions$genebody, file.path("cache", "regions_genebody.rds"))

message("Region counts:")
print(vapply(txdb_regions, length, integer(1)))
