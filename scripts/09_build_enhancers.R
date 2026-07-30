#!/usr/bin/env Rscript
# ============================================================================
# 09_build_enhancers.R   (system R)   --- Stage B ---
#
# Download ENCODE SCREEN V3 mm10 cCRE, keep distal enhancer-like (dELS), and
# assemble the combined regions_all.rds. Ported from
# 20260522_g4_gloop_extended_V3.Rmd chunks s2-enhancer-download / s2-save-regions.
# Network-dependent: tries URL BED download first, falls back to the UCSC REST
# JSON API per chromosome. Cached to cache/regions_enhancer.rds.
#
# Inputs:   cache/regions_{promoter,5UTR,intron1,genebody}.rds (from script 08)
# Outputs:  cache/regions_enhancer.rds, cache/regions_all.rds,
#           cache/encode_ccre_mm10_dELS.bed, cache/encode_ccre_download_info.txt
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("cache")

std_chroms <- cfg$std_chroms
enh_urls   <- as.character(unlist(cfg$regions$enhancer_urls))

enh_cache_rds <- file.path("cache", "regions_enhancer.rds")
enh_bed       <- file.path("cache", "encode_ccre_mm10_dELS.bed")
enh_info      <- file.path("cache", "encode_ccre_download_info.txt")

regions_enhancer <- cache_or_build(enh_cache_rds, {
  # --- Path 1: URL download -> BED file ---
  if (!file.exists(enh_bed)) {
    dl_status <- "fail"; used_url <- NA_character_
    for (u in enh_urls) {
      message("trying URL '", u, "'")
      target <- if (grepl("\\.gz$", u)) paste0(enh_bed, ".gz") else enh_bed
      ok <- tryCatch({ utils::download.file(u, target, mode = "wb", quiet = FALSE); TRUE },
                     error = function(e) { message("  -> ", conditionMessage(e)); FALSE })
      if (ok) {
        if (grepl("\\.gz$", target)) R.utils::gunzip(target, destname = enh_bed, overwrite = TRUE)
        dl_status <- "ok"; used_url <- u; break
      }
    }
    writeLines(c(paste("URL:", used_url), paste("Date:", Sys.time()),
                 paste("Status:", dl_status)), enh_info)
  }

  if (file.exists(enh_bed)) {
    bed <- read.table(enh_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                      comment.char = "#", quote = "")
    colnames(bed)[1:3] <- c("chr", "start", "end")
    bed$class <- bed[[ncol(bed)]]
    dels <- bed[grepl("dELS", bed$class), , drop = FALSE]
    message(sprintf("ENCODE cCRE rows: %d | dELS: %d", nrow(bed), nrow(dels)))
    if (nrow(dels) == 0) stop("No dELS rows found — inspect ", enh_bed)
    gr <- GenomicRanges::GRanges(seqnames = dels$chr,
      ranges = IRanges::IRanges(start = dels$start + 1, end = dels$end), strand = "*",
      accession = if (ncol(dels) >= 4) dels[[4]] else NA_character_, class = dels$class)
  } else {
    message("BED download failed — UCSC REST API fallback ...")
    query_chroms <- std_chroms
    ccre_rows <- lapply(query_chroms, function(chr) {
      url <- sprintf("https://api.genome.ucsc.edu/getData/track?genome=mm10&track=encodeCcreCombined&chrom=%s", chr)
      tryCatch({ resp <- jsonlite::fromJSON(url); df <- resp[["encodeCcreCombined"]]
        if (is.data.frame(df) && nrow(df) > 0) df else NULL },
        error = function(e) { message(chr, " FAILED: ", conditionMessage(e)); NULL })
    })
    ccre_df <- do.call(rbind, Filter(Negate(is.null), ccre_rows))
    if (is.null(ccre_df) || nrow(ccre_df) == 0) stop("UCSC REST returned no cCRE data.")
    dels_df <- ccre_df[grepl("dELS", ccre_df$ccre), ]
    if (nrow(dels_df) == 0) stop("No dELS cCREs found.")
    gr <- GenomicRanges::GRanges(seqnames = dels_df$chrom,
      ranges = IRanges::IRanges(start = dels_df$chromStart + 1L, end = dels_df$chromEnd),
      strand = "*", accession = dels_df$name, class = dels_df$ccre)
  }
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  GenomeInfoDb::keepSeqlevels(gr, intersect(GenomeInfoDb::seqlevels(gr), std_chroms),
                              pruning.mode = "coarse")
})

message("Enhancer (dELS) count: ", length(regions_enhancer))

# --- Assemble combined regions_all.rds (matches extended_V3 order) ----------
regions_all <- list(
  promoter = readRDS(file.path("cache", "regions_promoter.rds")),
  enhancer = regions_enhancer,
  `5UTR`    = readRDS(file.path("cache", "regions_5UTR.rds")),
  intron1  = readRDS(file.path("cache", "regions_intron1.rds")),
  genebody = readRDS(file.path("cache", "regions_genebody.rds")))
saveRDS(regions_all, file.path("cache", "regions_all.rds"))

message("regions_all counts:")
print(vapply(regions_all, length, integer(1)))
