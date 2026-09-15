#!/usr/bin/env Rscript
# ============================================================================
# 53_gene_body_expression_topology_metaprofiles.R  (system R) -- Stage N, pt 2 --
#
# Rule 52's figure re-cut by G4 topology: does the position of the signal across
# the gene depend on WHICH G4 structure the promoter carries?
#
# READS RULE 52'S CACHE AND OPENS NO BIGWIG. Every stratum here is a row subset
# of the (gene x bin) matrix rule 52 already built. If this script ever grows a
# bigWig input, that split has been broken.
#
# SAME COMPOSITE AXIS AS RULE 52 (v2, 2026-09-15): two fixed-bp anchor blocks
# (TSS +/-W and TES +/-W, 50 bp bins, REAL bp) either side of a scaled interior.
# It comes from the same composite_positions() / composite_axis() helpers rule
# 52 uses, so the two rules cannot disagree about where the TSS and TES sit. The
# v1 pure scale-regions axis was replaced because its 50 bp flank bins next to
# median 323 bp body bins put a 6.5x resolution step at the TSS, which drew a
# cliff that was binning rather than biology -- see the header of rule 52.
#
# TOPOLOGY IS ASSIGNED ON THE PROMOTER, not the gene body, so the labels are
# identical to rule 51's and the TSS-window and whole-gene views can be compared
# directly.
#
# Inputs:   cache/gene_expr_repmean_<assay>_f<F>_fb<NF>_bb<NB>.rds  (rule 52),
#           results/tables/gene_body_expression_genes.csv           (rule 52),
#           cfg$paths$promoters_rds, cache/peaks.rds,
#           results/tables/peak_topology.csv
# Outputs:  results/tables/gene_body_expression_topology_profiles.csv
#           results/tables/gene_body_expression_topology_counts.csv
#           results/tables/gene_body_expression_topology_shape_stats.csv
#           results/figures/25_tss_tes_metaprofiles/
#             gene_signal_by_expression_topology.{pdf,png}
#             gene_signal_by_expression_topology_shape.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

fig_subdir <- "25_tss_tes_metaprofiles"
fig_root   <- "results/figures"
ensure_dirs("results/tables", file.path(fig_root, fig_subdir))

mp        <- cfg$metaprofile
genotypes <- cfg$genotypes
W         <- as.integer(mp$expr_gene_anchor_window %||% 2000L)
abin      <- as.integer(mp$expr_gene_anchor_bin    %||% 50L)
n_bbins   <- as.integer(mp$expr_gene_body_bins     %||% 60L)
n_q       <- as.integer(mp$expr_n_quantiles %||% 4)
incl_sil  <- isTRUE(mp$expr_include_silent %||% TRUE)
min_n     <- as.integer(mp$expr_min_regions_per_panel %||% 1L)
q_levels  <- paste0("Q", seq_len(n_q))
grp_levels <- c(if (incl_sil) "silent", q_levels)

cache_file <- sprintf("cache/gene_composite_repmean_%s_w%d_ab%d_bb%d.rds",
                      cfg$assay, W, abin, n_bbins)

# --- 1. Rule 52's cache, validated ------------------------------------------
if (!file.exists(cache_file))
  stop("Missing ", cache_file, ". Run rule gene_body_expression_metaprofiles (52) first.")
cached <- readRDS(cache_file)
if (!identical(cached$anchor_window, W) || !identical(cached$anchor_bin, abin) ||
    !identical(cached$n_body_bins, n_bbins) || !identical(cached$assay, cfg$assay) ||
    !setequal(names(cached$rep_mean), genotypes))
  stop("Cache ", cache_file, " was built for a different anchor window / binning / ",
       "assay / genotype set than the current config. Delete it and re-run rule 52.")
genes     <- cached$regions          # the exact GRanges profiled, in its own order
# The composite layout: two fixed-bp anchor blocks either side of a scaled
# interior. Taken from the same helper rule 52 used, so the two rules cannot
# disagree about where the TSS and the TES are drawn.
layout    <- composite_positions(W, abin, n_bbins)
positions <- layout$x
axis_spec <- composite_axis(W)
stopifnot(identical(names(genes), rownames(cached$rep_mean[[genotypes[1]]])))

# --- 2. Topology of the strongest union peak at each gene's PROMOTER --------
prom <- readRDS(cfg$paths$promoters_rds)
GenomeInfoDb::seqlevelsStyle(prom) <- "UCSC"
prom <- GenomeInfoDb::keepSeqlevels(
  prom, intersect(GenomeInfoDb::seqlevels(prom), cfg$std_chroms), pruning.mode = "coarse")

union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
prom <- assign_region_topology(prom, union)

pix <- match(names(genes), names(prom))
if (anyNA(pix))
  stop(sum(is.na(pix)), " profiled gene(s) have no promoter in ",
       cfg$paths$promoters_rds, ". The cache and the promoter set are out of step; ",
       "re-run rule 50 then rule 52.")
gene_topo <- as.character(prom$group)[pix]

# Deliberately NOT topology_keep_groups(): it returns c(definite, "mixed"?,
# "none") and never includes "no_canonical_PQS" (helpers.R:76-80), which would
# silently drop those genes from the denominator -- the exact failure class the
# v1.1 audit was about. Build the keep list from the levels actually present, so
# the counts table sums to the profiled gene set.
keep_groups <- intersect(c(topology_levels, "none"), unique(gene_topo))
topo_labs   <- topology_group_labels(keep_groups)
message("Topology strata present at these genes' promoters: ",
        paste(keep_groups, collapse = ", "))
print(table(factor(gene_topo, levels = keep_groups)))

# --- 3. Expression groups from rule 52 --------------------------------------
gt <- readr::read_csv("results/tables/gene_body_expression_genes.csv",
        col_types = readr::cols(entrez = readr::col_character(),
                                .default = readr::col_guess()))
ix <- match(names(genes), gt$entrez)
if (anyNA(ix))
  stop(sum(is.na(ix)), " cached genes are absent from ",
       "gene_body_expression_genes.csv. The cache and the gene table are out of ",
       "step; re-run rule 52.")
gt <- gt[ix, ]
group_col <- function(g) {
  cn <- paste0("group_", g)
  as.character(if (isTRUE(mp$expr_per_genotype %||% TRUE) && cn %in% names(gt)) gt[[cn]]
               else gt$group_fixed)
}

# --- 4. Profiles per (genotype x topology x expression group) ---------------
count_rows <- list(); prof_rows <- list()
for (g in genotypes) {
  grp <- group_col(g)
  for (tp in keep_groups) {
    in_topo <- !is.na(gene_topo) & gene_topo == tp
    for (lv in grp_levels) {
      idx <- which(in_topo & !is.na(grp) & grp == lv)
      # n recorded BEFORE the guard, so a stratum too small to profile still
      # reports its true size instead of vanishing from the counts table too.
      count_rows[[length(count_rows) + 1]] <- data.frame(
        genotype = g, topology = tp, expr_group = lv, n_genes = length(idx),
        profiled = length(idx) >= min_n, stringsAsFactors = FALSE)
      if (length(idx) < min_n) next
      p <- profile_from_rep_mean(cached$rep_mean[[g]], idx, positions = positions)
      if (is.null(p)) next
      p$genotype <- g; p$topology <- tp; p$expr_group <- lv
      prof_rows[[length(prof_rows) + 1]] <- p
    }
  }
}
counts <- dplyr::bind_rows(count_rows)
readr::write_csv(counts, "results/tables/gene_body_expression_topology_counts.csv")
if (!length(prof_rows))
  stop("No (genotype x topology x expression) stratum reached ",
       "metaprofile.expr_min_regions_per_panel = ", min_n, ".")

profiles <- dplyr::bind_rows(prof_rows) %>% dplyr::rename(n_genes = n_regions)
profiles$genotype   <- factor(profiles$genotype, levels = genotypes)
profiles$topology   <- factor(profiles$topology, levels = keep_groups)
profiles$expr_group <- factor(profiles$expr_group, levels = grp_levels)
profiles <- profiles %>%
  dplyr::group_by(genotype, topology, expr_group) %>%
  dplyr::mutate(window_mean = mean(mean, na.rm = TRUE),
                shape     = mean / window_mean,
                shape_sem = sem  / window_mean) %>%
  dplyr::ungroup()
# profile_from_rep_mean() calls its axis column `position`, but on a composite
# axis that value is PLOT SPACE, not a genomic position -- the real distance is
# `bp`, and it is NA across the interior.
profiles <- dplyr::left_join(dplyr::rename(profiles, x = position),
                             layout[, c("bin", "block", "bp", "frac")], by = "bin")
readr::write_csv(profiles, "results/tables/gene_body_expression_topology_profiles.csv")

# Block-aware shares. On the composite axis "position < 0" no longer means
# "upstream of the TSS" -- x is bp inside the anchor blocks and plot space in
# the interior -- so the split is taken from the block each bin belongs to.
# Bin-weighted, not bp-weighted (an interior bin covers more genome than an
# anchor bin), so read these as "share of the drawn curve".
block_frac <- function(d) {
  m <- pmax(d$mean, 0); m[!is.finite(m)] <- 0
  tot <- sum(m)
  if (tot <= 0) return(data.frame(frac_tss_block = NA_real_, frac_interior = NA_real_,
                                  frac_tes_block = NA_real_))
  data.frame(frac_tss_block = sum(m[d$block == "tss"]) / tot,
             frac_interior  = sum(m[d$block == "interior"]) / tot,
             frac_tes_block = sum(m[d$block == "tes"]) / tot)
}
shape_stats <- profiles %>%
  dplyr::group_by(genotype, topology, expr_group) %>%
  dplyr::group_modify(~ {
    d <- .x[order(.x$bin), ]
    tss_b <- d[d$block == "tss", ]; int_b <- d[d$block == "interior", ]
    tes_b <- d[d$block == "tes", ]
    floor_v <- mean(int_b$mean, na.rm = TRUE)
    pk <- which.max(tss_b$mean)
    pk_bp <- if (length(pk)) tss_b$bp[pk] else NA_real_
    pk_v  <- if (length(pk)) tss_b$mean[pk] else NA_real_
    half  <- floor_v + (pk_v - floor_v) / 2
    dn    <- tss_b[!is.na(tss_b$bp) & tss_b$bp > pk_bp, ]
    hd    <- if (nrow(dn) && any(dn$mean <= half, na.rm = TRUE))
               dn$bp[which(dn$mean <= half)[1]] else NA_real_
    ti <- which.max(tes_b$mean)
    cbind(data.frame(n_genes = d$n_genes[1], tss_peak_bp = pk_bp,
                     tss_peak_value = pk_v, tss_half_decay_bp = hd,
                     interior_mean = floor_v,
                     tes_peak_bp = if (length(ti)) tes_b$bp[ti] else NA_real_,
                     tes_over_interior = if (length(ti)) tes_b$mean[ti] / floor_v
                                         else NA_real_,
                     window_mean = mean(d$mean, na.rm = TRUE)),
          block_frac(d))
  }) %>%
  dplyr::ungroup()
readr::write_csv(shape_stats,
                 "results/tables/gene_body_expression_topology_shape_stats.csv")

# --- 5. Figures --------------------------------------------------------------
lty_vals <- c(setNames(rep("solid", n_q), q_levels),
              silent = as.character(mp$expr_silent_linetype %||% "22"))
expr_cols <- c(
  setNames(grDevices::colorRampPalette(
    as.character(unlist(mp$expr_palette %||%
      c("#2166AC", "#4393C3", "#D6604D", "#B2182B"))))(n_q), q_levels),
  silent = as.character(mp$expr_silent_colour %||% "grey55"))
scale_note <- paste0("x is linear in bp WITHIN each anchor block; only the shaded ",
                     "interior is scaled, and no distance may be read off it")

# Per-panel n, built from the COUNTS table (not the profiles) so a stratum that
# was too small to draw still contributes its label. With scales = "free_y" an
# empty panel has no y-scale, so a geom_text at y = Inf cannot render there --
# an unlabelled blank panel therefore means "no curve reached the threshold",
# which is what the counts table is for.
panel_lab <- counts %>%
  dplyr::mutate(genotype = factor(genotype, levels = genotypes),
                topology = factor(topology, levels = keep_groups)) %>%
  dplyr::filter(!is.na(topology)) %>%
  dplyr::group_by(genotype, topology) %>%
  dplyr::summarise(label = paste0("n=", format(sum(n_genes), big.mark = ",", trim = TRUE)),
                   .groups = "drop")

build_topo_fig <- function(ycol, semcol, fname, title, subtitle) {
  p <- ggplot(profiles, aes(x, .data[[ycol]], colour = expr_group,
                            fill = expr_group)) +
    annotate("rect", xmin = axis_spec$interior[1], xmax = axis_spec$interior[2],
             ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.45) +
    geom_vline(xintercept = axis_spec$dividers, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = .data[[ycol]] - .data[[semcol]],
                    ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = expr_group), linewidth = 0.6) +
    geom_vline(xintercept = axis_spec$anchors, linetype = "dashed", colour = "grey30") +
    geom_text(data = panel_lab, aes(x = -Inf, y = Inf, label = label),
              inherit.aes = FALSE, hjust = -0.15, vjust = 1.5, size = 2.2,
              colour = "grey25") +
    # facet_wrap, not facet_grid: facet_grid(scales = "free_y") frees y per ROW,
    # which would put a 4,000-gene "none" panel on the same scale as a 20-gene
    # antiparallel one and flatten the latter to a line.
    facet_wrap(vars(genotype, topology), scales = "free_y",
               ncol = length(keep_groups), drop = FALSE,
               labeller = labeller(topology = function(x)
                 ifelse(x %in% names(topo_labs), topo_labs[x], x))) +
    scale_x_continuous(breaks = axis_spec$breaks, labels = axis_spec$labels,
                       limits = axis_spec$limits, expand = ggplot2::expansion(mult = 0.004)) +
    scale_colour_manual(values = expr_cols, name = "Expression", drop = FALSE) +
    scale_fill_manual(values = expr_cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = lty_vals, guide = "none", drop = FALSE) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_pub() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text = element_text(size = 7))
  save_plot(p, fname, fig_subdir, fig_root,
            width = 3 + 2.2 * length(keep_groups), height = 3.2 * length(genotypes))
}

sub_base <- sprintf(paste0("anchors +/-%s kb in %d bp bins (REAL bp), interior scaled ",
                           "to %d bins; topology = strongest union G4 peak at the ",
                           "PROMOTER (identical to rule 51's label, so the two are ",
                           "comparable)"),
                    format(W / 1000, trim = TRUE), abin, n_bbins)

build_topo_fig("mean", "sem", "gene_signal_by_expression_topology",
  "G4 signal TSS to TES, by expression quartile and promoter G4 topology",
  sprintf("%s\n%s\nfree y PER PANEL; curves below n = %d are not drawn but are counted in gene_body_expression_topology_counts.csv",
          sub_base, scale_note, min_n))

if (isTRUE(mp$expr_shape_normalise %||% TRUE))
  build_topo_fig("shape", "shape_sem", "gene_signal_by_expression_topology_shape",
    "G4 signal SHAPE TSS to TES, by expression quartile and promoter G4 topology",
    sprintf("%s; each curve divided by its own window mean, so only POSITION differs\n%s",
            sub_base, scale_note))

# --- 6. Report ---------------------------------------------------------------
message("\nGenes per (genotype, topology, expression group):")
print(as.data.frame(counts))
message("\nStrata not drawn (below n = ", min_n, "): ",
        sum(!counts$profiled), " of ", nrow(counts))
message("Done (no bigWig was opened).")
