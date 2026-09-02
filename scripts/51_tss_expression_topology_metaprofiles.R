#!/usr/bin/env Rscript
# ============================================================================
# 51_tss_expression_topology_metaprofiles.R   (env: r_g4)  --- Stage M, point 2 -
#
# Rule 50's figure, additionally split by the G4SP topology of the strongest
# union G4 peak overlapping each promoter (mixed/unresolved included). Answers:
# does the expression-vs-position relationship hold within a G4 topology class?
#
# OPENS NO BIGWIG. It reads the (promoter x bin) matrices rule 50 cached, so a
# topology re-cut costs ~20 s instead of the ~30 min bigWig pass. If this script
# ever grows a bigWig input, the 50/51 split has been broken.
#
# Topology strata are SPARSE at promoters (measured 2026-09-01: parallel 4202,
# mixed 1342, hybrid 143, antiparallel 68, no_canonical_PQS 56, none 18653).
# Split five ways by expression that is ~14 promoters per antiparallel panel, so
# every panel carries its own n, empty panels stay visible (drop = FALSE), and
# curves below metaprofile.expr_min_regions_per_panel are not drawn at all --
# but their n is recorded BEFORE the guard, so a blank panel is explained rather
# than absent (the rule 31 precedent).
#
# Inputs:   cache/tss_expr_repmean_<assay>_hw<HW>_b<N>.rds  (rule 50),
#           results/tables/tss_expression_gene_groups.csv   (rule 50),
#           cache/peaks.rds, results/tables/peak_topology.csv
# Outputs:  results/tables/tss_expression_topology_profiles.csv
#           results/tables/tss_expression_topology_counts.csv
#           results/tables/tss_expression_topology_shape_stats.csv
#           results/figures/24_expression_metaprofiles/
#             tss_signal_by_expression_topology.{pdf,png}
#             tss_signal_by_expression_topology_shape.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(S4Vectors)
  library(dplyr); library(readr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yaml"
source(file.path("scripts", "_shared", "helpers.R"))
cfg <- load_config(config_path)

fig_subdir <- "24_expression_metaprofiles"
fig_root   <- "results/figures"
ensure_dirs("results/tables", file.path(fig_root, fig_subdir))

mp         <- cfg$metaprofile
genotypes  <- cfg$genotypes
half_width <- mp$half_width %||% 2000
n_bins     <- mp$n_bins %||% 80
n_q        <- as.integer(mp$expr_n_quantiles %||% 4)
min_n      <- as.integer(mp$expr_min_regions_per_panel %||% 20)
incl_silent <- isTRUE(mp$expr_include_silent %||% TRUE)
grp_levels <- c(if (incl_silent) "silent", paste0("Q", seq_len(n_q)))

cache_file <- sprintf("cache/tss_expr_repmean_%s_hw%d_b%d.rds",
                      cfg$assay, half_width, n_bins)

# --- 1. Rule 50's cache, validated ------------------------------------------
if (!file.exists(cache_file))
  stop("Missing ", cache_file, ". Run rule tss_expression_metaprofiles first.")
cached <- readRDS(cache_file)
if (!identical(cached$half_width, half_width) || !identical(cached$n_bins, n_bins) ||
    !identical(cached$assay, cfg$assay) || !setequal(names(cached$rep_mean), genotypes))
  stop("Cache ", cache_file, " was built for a different window / bin count / assay / ",
       "genotype set than the current config. Delete it and re-run rule 50.")
win       <- cached$regions          # the exact GRanges profiled, in its own order
positions <- cached$positions

# --- 2. Topology of the strongest overlapping union peak --------------------
union <- readRDS("cache/peaks.rds")$union
topo  <- readr::read_csv("results/tables/peak_topology.csv", show_col_types = FALSE)
union$topology <- as.character(topo$topology_final[match(names(union), topo$peak_id)])
win <- assign_region_topology(win, union)        # adds mcols$group, cache order kept
stopifnot(identical(names(win), rownames(cached$rep_mean[[genotypes[1]]])))

# Deliberately NOT topology_keep_groups(): it returns c(definite, "mixed"?, "none")
# and never includes "no_canonical_PQS" (helpers.R:76-80), which would silently drop
# those 56 promoters from the denominator -- the exact failure class the v1.1 audit
# was about. Build the keep list from the levels actually present instead, so the
# counts table sums to the promoter universe.
keep_groups <- intersect(c(topology_levels, "none"), unique(win$group))
topo_labs   <- topology_group_labels(keep_groups)
message("Topology strata present at promoters: ", paste(keep_groups, collapse = ", "))
print(table(factor(win$group, levels = keep_groups)))

# --- 3. Expression groups from rule 50 --------------------------------------
# entrez forced to character: readr would otherwise coerce the IDs to numeric and
# the match() against the matrix rownames would half-succeed.
gg <- readr::read_csv("results/tables/tss_expression_gene_groups.csv",
                      col_types = readr::cols(entrez = readr::col_character(),
                                              .default = readr::col_guess()))
ix <- match(names(win), gg$entrez)
if (anyNA(ix))
  stop(sum(is.na(ix)), " cached promoters are absent from tss_expression_gene_groups.csv. ",
       "The cache and the gene-group table are out of step; re-run rule 50.")
gg <- gg[ix, ]

# Rule 50's primary arm: per-genotype cuts if configured, else the fixed column.
group_col <- function(g) {
  cn <- paste0("group_", g)
  as.character(if (isTRUE(mp$expr_per_genotype %||% TRUE) && cn %in% names(gg)) gg[[cn]]
               else gg$group_fixed)
}

# --- 4. Profiles per (genotype x topology x expression group) ---------------
count_rows <- list(); prof_rows <- list()
for (g in genotypes) {
  grp <- group_col(g)
  for (tp in keep_groups) {
    in_topo <- !is.na(win$group) & win$group == tp
    for (lv in grp_levels) {
      idx <- which(in_topo & !is.na(grp) & grp == lv)
      # n recorded BEFORE the guard, so a stratum too small to profile still
      # reports its true size instead of vanishing from the counts table too.
      count_rows[[length(count_rows) + 1]] <- data.frame(
        genotype = g, topology = tp, expr_group = lv, n_promoters = length(idx),
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
counts$topology   <- factor(counts$topology, levels = keep_groups)
counts$genotype   <- factor(counts$genotype, levels = genotypes)
counts$expr_group <- factor(counts$expr_group, levels = grp_levels)
readr::write_csv(counts, "results/tables/tss_expression_topology_counts.csv")

if (length(prof_rows) == 0)
  stop("No topology x expression stratum reached metaprofile.expr_min_regions_per_panel = ", min_n)
profiles <- dplyr::bind_rows(prof_rows) %>% dplyr::rename(n_promoters = n_regions)
profiles$topology   <- factor(profiles$topology, levels = keep_groups)
profiles$genotype   <- factor(profiles$genotype, levels = genotypes)
profiles$expr_group <- factor(profiles$expr_group, levels = grp_levels)
profiles <- profiles %>%
  dplyr::group_by(genotype, topology, expr_group) %>%
  dplyr::mutate(window_mean = mean(mean, na.rm = TRUE),
                shape     = mean / window_mean,
                shape_sem = sem  / window_mean) %>%
  dplyr::ungroup()
readr::write_csv(profiles, "results/tables/tss_expression_topology_profiles.csv")

shape_stats <- profiles %>%
  dplyr::group_by(genotype, topology, expr_group) %>%
  dplyr::group_modify(~ cbind(data.frame(n_promoters = .x$n_promoters[1]),
                              profile_shape_stats(.x))) %>%
  dplyr::ungroup()
readr::write_csv(shape_stats, "results/tables/tss_expression_topology_shape_stats.csv")

# --- 5. Figures --------------------------------------------------------------
# facet_wrap, NOT facet_grid: facet_grid(scales = "free_y") frees y per ROW, but
# these panels differ in magnitude cell by cell (parallel n=4202 vs antiparallel
# n=68), so each needs its own scale.
lty_vals <- c(setNames(rep("solid", n_q), paste0("Q", seq_len(n_q))),
              silent = as.character(mp$expr_silent_linetype %||% "22"))
# Same ramp as rule 50, from the same config keys, so the two figures are read
# with one legend: blue = least active -> red = most active, "silent" off-ramp.
expr_cols <- c(
  setNames(grDevices::colorRampPalette(
    as.character(unlist(mp$expr_palette %||%
      c("#2166AC", "#4393C3", "#D6604D", "#B2182B"))))(n_q), paste0("Q", seq_len(n_q))),
  silent = as.character(mp$expr_silent_colour %||% "grey55"))
n_lab <- counts %>%
  dplyr::group_by(genotype, topology) %>%
  dplyr::summarise(label = paste0("n = ", paste(n_promoters, collapse = "/")), .groups = "drop")

build_topo_fig <- function(ycol, semcol, fname, title, subtitle) {
  p <- ggplot(profiles, aes(position, .data[[ycol]], colour = expr_group, fill = expr_group)) +
    geom_ribbon(aes(ymin = .data[[ycol]] - .data[[semcol]],
                    ymax = .data[[ycol]] + .data[[semcol]]), alpha = 0.15, colour = NA) +
    geom_line(aes(linetype = expr_group), linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.4, size = 2.2, colour = "grey25") +
    facet_wrap(~ genotype + topology, scales = "free_y", drop = FALSE,
               ncol = length(keep_groups),
               labeller = labeller(topology = topo_labs)) +
    scale_colour_manual(values = expr_cols, name = "Expression", drop = FALSE) +
    scale_fill_manual(values = expr_cols, guide = "none", drop = FALSE) +
    scale_linetype_manual(values = lty_vals, guide = "none", drop = FALSE) +
    labs(x = "Distance to TSS (bp)", y = NULL, title = title, subtitle = subtitle) +
    theme_pub() + theme(legend.position = "bottom")
  save_plot(p, fname, fig_subdir, fig_root,
            width = 3.0 * length(keep_groups) + 1, height = 3.1 * length(genotypes) + 1.5)
}

sub_base <- sprintf(paste("+/- %d bp, %d bins; topology = strongest overlapping union G4 peak;",
                          "n per panel is %s in group order; curves under n = %d are not drawn"),
                    half_width, n_bins, paste(grp_levels, collapse = "/"), min_n)
build_topo_fig("mean", "sem", "tss_signal_by_expression_topology",
  "Promoter G4 signal at the TSS, by expression quartile and G4 topology",
  paste0(sub_base, "\nfree y per panel; ribbon = SEM across promoters, not a test of the difference"))
if (isTRUE(mp$expr_shape_normalise %||% TRUE))
  build_topo_fig("shape", "shape_sem", "tss_signal_by_expression_topology_shape",
    "Promoter G4 signal SHAPE at the TSS, by expression quartile and G4 topology",
    paste0(sub_base, "\neach curve divided by its own window mean, so only POSITION differs"))

message("\nPromoters per (genotype, topology, expression group):")
print(as.data.frame(counts))
skipped <- counts %>% dplyr::filter(!profiled, n_promoters > 0)
if (nrow(skipped) > 0) {
  message("\nStrata below expr_min_regions_per_panel = ", min_n,
          " (panel shown, no curve drawn):")
  print(as.data.frame(skipped))
}
message("\nDone. No bigWigs were opened.")
