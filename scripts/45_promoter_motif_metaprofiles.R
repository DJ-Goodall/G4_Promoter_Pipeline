#!/usr/bin/env Rscript
# ============================================================================
# 30_promoter_motif_metaprofiles.R   (env: r_g4)   --- Phase 13, Analysis 2 ---
#
# Motif-resolution promoter metaprofiles. Unlike rule 08 (one topology per peak,
# uniform bins over a window), here EVERY G4 motif (rule 27) within TSS +/- 2 kb is
# a data point carrying (a) its own topology, (b) a strand-aware signed distance to
# its nearest TSS, and (c) a per-genotype CUT&Tag signal (motif centre +/- 100 bp).
# Two views of how topology depends on position relative to the TSS:
#   (1) SIGNAL  : mean per-motif G4 signal vs distance-to-TSS, coloured by topology,
#                 faceted by genotype.
#   (2) DENSITY : number of motifs vs distance-to-TSS, coloured by topology.
#
# Also caches the promoter-motif signal table (cache/promoter_motif_signal.rds) that
# rule 31 (DEG split) reuses, so the bigWig signal extraction runs only once.
#
# Inputs:   cache_v2 regions_promoter.rds (Entrez names), cache/motifs_all.rds,
#           G4 bigWigs
# Outputs:  results/tables/promoter_motif_profiles.csv,
#           cache/promoter_motif_signal.rds
#           results/figures/10_metaprofiles/{motif_signal_by_topology_tss,
#             motif_density_by_topology_tss}.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(IRanges); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)
ensure_dirs("results/tables", "results/figures/10_metaprofiles", "cache")
fig_root <- "results/figures"

genotypes <- cfg$genotypes
mcfg      <- cfg$motif_analysis
hw        <- mcfg$promoter_half_width %||% 2000
sig_hw    <- mcfg$signal_half_width %||% 100
n_bins    <- mcfg$n_bins %||% cfg$metaprofile$n_bins %||% 80
definite  <- c("parallel", "antiparallel", "hybrid")

# --- Promoters -> strand-aware TSS points + windows -------------------------
prom <- readRDS(file.path(cfg$paths$cache_v2, cfg$feature$region_files$promoter))
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), cfg$std_chroms), pruning.mode = "coarse")
prom_entrez <- names(prom); if (is.null(prom_entrez)) prom_entrez <- as.character(S4Vectors::mcols(prom)$gene_id)
tss      <- GenomicRanges::resize(prom, width = 1, fix = "center")     # TSS, strand kept
tss_pos  <- GenomicRanges::start(tss)
tss_strand <- as.character(GenomicRanges::strand(tss))
prom_win <- GenomicRanges::resize(tss, width = 2 * hw, fix = "center")

# --- Motifs (centres) within a promoter window, assigned to NEAREST TSS -----
motifs_gr  <- readRDS("cache/motifs_all.rds")
motif_topo <- as.character(S4Vectors::mcols(motifs_gr)$topology)
motif_cen  <- as.integer(S4Vectors::mcols(motifs_gr)$motif_center)
motif_pts  <- GenomicRanges::GRanges(GenomicRanges::seqnames(motifs_gr),
                                     IRanges::IRanges(motif_cen, width = 1))

ov <- GenomicRanges::findOverlaps(motif_pts, prom_win, ignore.strand = TRUE)
qh <- S4Vectors::queryHits(ov); sh <- S4Vectors::subjectHits(ov)
strand_sh <- tss_strand[sh]
signed <- ifelse(strand_sh == "-", tss_pos[sh] - motif_cen[qh], motif_cen[qh] - tss_pos[sh])
# nearest TSS per motif: smallest |dist|, tie-break lowest Entrez
ord <- order(qh, abs(signed), suppressWarnings(as.integer(prom_entrez[sh])))
first <- !duplicated(qh[ord])
sel <- ord[first]

pm <- data.frame(
  motif_id   = names(motifs_gr)[qh[sel]],
  peak_id    = as.character(S4Vectors::mcols(motifs_gr)$peak_id)[qh[sel]],
  chr        = as.character(GenomicRanges::seqnames(motifs_gr))[qh[sel]],
  motif_center = motif_cen[qh[sel]],
  topology   = motif_topo[qh[sel]],
  prom_entrez = prom_entrez[sh[sel]],
  dist       = as.integer(signed[sel]),
  stringsAsFactors = FALSE)
message(sprintf("Motifs within TSS +/- %d bp: %d (assigned to nearest of %d promoters)",
                hw, nrow(pm), length(prom)))

# --- Per-genotype per-motif signal (motif centre +/- sig_hw), cached ---------
sig_win <- GenomicRanges::GRanges(pm$chr,
             IRanges::IRanges(pmax(1L, pm$motif_center - sig_hw), pm$motif_center + sig_hw))
bw_meta <- build_bw_meta(cfg$paths$bigwig_dir, genotypes = genotypes, assay = cfg$assay)
for (g in genotypes) {
  bws <- bw_meta$filepath[bw_meta$genotype == g]
  if (length(bws) == 0) { pm[[paste0("signal_", g)]] <- NA_real_; next }
  cache_f <- sprintf("cache/promoter_motif_sig_%s_n%d_hw%d.rds", g, nrow(pm), sig_hw)
  pm[[paste0("signal_", g)]] <- cache_or_build(cache_f, {
    message("Promoter-motif signal: ", g, " (", length(bws), " bigWigs)")
    v <- mean_replicate_signal(sig_win, bws); evict_bigwigs(bws); v
  })
}
saveRDS(pm, "cache/promoter_motif_signal.rds")

# --- Bin signed distance, build signal + density profiles -------------------
breaks    <- seq(-hw, hw, length.out = n_bins + 1)
centres   <- (breaks[-1] + breaks[-(n_bins + 1)]) / 2
pm$bin    <- cut(pm$dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
pm <- pm[!is.na(pm$bin) & pm$topology %in% definite, ]
pm$position <- centres[pm$bin]

sig_long <- bind_rows(lapply(genotypes, function(g) {
  col <- paste0("signal_", g)
  pm %>% transmute(genotype = g, topology, bin, position, signal = .data[[col]])
}))
profiles <- sig_long %>%
  group_by(genotype, topology, bin, position) %>%
  summarise(n_motifs = dplyr::n(),
            mean = mean(signal, na.rm = TRUE),
            sem  = sd(signal, na.rm = TRUE) / sqrt(sum(!is.na(signal))),
            .groups = "drop")
profiles$genotype <- factor(profiles$genotype, levels = genotypes)
profiles$topology <- factor(profiles$topology, levels = definite)
readr::write_csv(profiles, "results/tables/promoter_motif_profiles.csv")

topo_labs <- setNames(topology_labels[match(definite, topology_levels)], definite)

# --- (1) SIGNAL: mean per-motif signal vs distance to TSS, by topology -------
p_sig <- ggplot(profiles, aes(position, mean, colour = topology, fill = topology)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ genotype, nrow = 1) +
  scale_colour_manual(values = topology_palette, labels = topo_labs, name = "Topology") +
  scale_fill_manual(values = topology_palette, labels = topo_labs, name = "Topology") +
  labs(x = "Distance to TSS (bp)", y = "Mean per-motif G4 signal",
       title = "Promoter G4 signal at individual motifs, by topology and position",
       subtitle = sprintf("Every motif within TSS +/- %d bp; signal = motif centre +/- %d bp; strand-aware",
                          hw, sig_hw)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_sig, "motif_signal_by_topology_tss", "10_metaprofiles", fig_root, width = 13, height = 4.5)

# --- (2) DENSITY: number of motifs vs distance to TSS, by topology -----------
dens <- profiles %>% distinct(topology, bin, position, n_motifs)
p_den <- ggplot(dens, aes(position, n_motifs, colour = topology)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = topology_palette, labels = topo_labs, name = "Topology") +
  labs(x = "Distance to TSS (bp)", y = "G4 motifs (count per bin)",
       title = "Where G4 topologies sit relative to the TSS",
       subtitle = sprintf("Motif positional density, TSS +/- %d bp, %d bins (genotype-independent)",
                          hw, n_bins)) +
  theme_pub() + theme(legend.position = "bottom")
save_plot(p_den, "motif_density_by_topology_tss", "10_metaprofiles", fig_root, width = 8, height = 5)

message("Done. Analysis 2 (promoter motif metaprofiles: signal + density).")
