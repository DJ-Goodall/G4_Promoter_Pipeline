# ============================================================================
# G4_Promoter_Pipeline — one Snakemake workflow covering:
#   * RNA-seq QC / DESeq2 / enrichment / RNA x CUT&Tag integration
#   * G4 (BG4) and R-loop (S9.6) CUT&Tag peak calling, region building,
#     genome-wide and regional enrichment
#   * sequence-level G4 topology (pqsfinder + DeepG4 + G4ShapePredictor) and all
#     downstream metaprofile / motif / lncRNA / G4 x R-loop analyses
#   * strand asymmetry + Gviz locus tracks
#
# Peaks + regions are recomputed here (Stage B) into ./cache; everything else
# reads from ./cache. One rule = one numbered script in scripts/. R rules run in
# system R 4.5.1; the G4ShapePredictor python rules use `conda run -n g4sp`; the
# DeepG4 rule reaches TensorFlow via reticulate -> the tf-new conda env.
#
# Run from the repository root with the Snakemake runner env active:
#   conda activate G4
#   snakemake -n                 # dry run (DAG)
#   snakemake --cores 4          # full run  (NO --use-conda; see README)
# ============================================================================
import glob, os

configfile: "config/config.yaml"

CONFIG   = "config/config.yaml"
BW_DIR   = config["paths"]["bigwig_dir"]
G4_BW    = sorted(glob.glob(os.path.join(BW_DIR, "*_G4_CnT_*_batch2_R*.bw")))   # 12 G4 bigWigs
RLOOP_BW = sorted(glob.glob(os.path.join(BW_DIR, "*_Rloop_CnT_mES_*_R*.bw")))   # 8 R-loop bigWigs
G4SP_ENV = config.get("g4sp", {}).get("conda_env", "g4sp")
COUNTS   = config["rnaseq"]["count_table"]
GENO     = config["genotypes"]
PK_GENO  = config["peak_calling"]["genotypes"]
PK_ASSAY = config["peak_calling"]["assays"]
CONTRASTS = ["DHX36KO_vs_WT","FANCJKO_vs_WT","dKO_vs_WT","dKO_vs_DHX36KO","dKO_vs_FANCJKO"]
KOS = ["DHX36KO","FANCJKO","dKO"]

# Stage-B peak/region caches (recomputed, not hand-built). Precompute all
# expand() lists at module scope (multi-line expand inside a rule body confuses
# Snakemake's rule parser, so rules reference these single tokens instead).
PEAK_RDS = (expand("cache/peaks_{a}_{g}.rds", a=PK_ASSAY, g=PK_GENO)
            + expand("cache/peaks_{a}_union.rds", a=PK_ASSAY)
            + expand("cache/peaks_{a}_ercc_union.rds", a=PK_ASSAY))
REGION_RDS = expand("cache/regions_{x}.rds",
                    x=["promoter", "5UTR", "intron1", "genebody", "enhancer", "all"])
DESEQ_TABS = expand("results/tables/deseq2_results_{c}.tsv", c=CONTRASTS)
VOLCANOES  = expand("results/figures/02_deseq2/volcano_{c}.pdf", c=CONTRASTS)
PROM_BOX   = expand("results/figures/07_integration/promoter_signal_boxplot_{ko}_vs_WT.pdf", ko=KOS)
TSS_COMB   = expand("results/figures/07_integration/tss_profile_combined_{ko}_vs_WT.pdf", ko=KOS)


rule all:
    input:
        # Stage A — RNA-seq
        "results/figures/01_qc/pca_plot.pdf",
        VOLCANOES,
        "results/figures/02_deseq2/union_DEG_heatmap.pdf",
        "results/figures/02_deseq2/upset_DEGs.pdf",
        "results/tables/deseq2_KO_vs_WT.csv",
        "results/tables/go_enrichment_summary.tsv",
        "results/figures/16_volcano/volcano_DHX36KO_top50.pdf",
        # Stage B — peaks + regions
        "results/tables/peak_summary.csv",
        "results/figures/04_peaks/01_peak_width_distribution.pdf",
        "results/figures/04_peaks/01_regions_by_chromosome.pdf",
        # Stage C — genome-wide CUT&Tag QC
        "results/figures/05_cuttag/violin_signal_G4_BG4.pdf",
        "results/figures/05_cuttag/G4_vs_Rloop_correlation_panel.pdf",
        # Stage D — regional enrichment
        "results/figures/06_regional/02_peak_counts_stacked.pdf",
        "results/figures/06_regional/02_region_fold_enrichment.pdf",
        "results/figures/06_regional/02_region_signal_boxplot.pdf",
        "results/figures/06_regional/02_region_metaprofile.pdf",
        "results/figures/06_regional/02_region_metaprofile_full_peaks.pdf",
        # Stage E — integration
        PROM_BOX,
        TSS_COMB,
        # Stage F — topology core
        "results/tables/peak_topology.csv",
        "results/tables/sanity_metrics.csv",
        "results/tables/topology_composition.csv",
        "results/figures/08_topology/topology_counts.pdf",
        "results/figures/08_topology/deepg4_roc.pdf",
        "results/figures/09_gain_loss/differential_class_fraction.pdf",
        "results/figures/10_metaprofiles/tss_signal_by_topology.pdf",
        "results/figures/11_propensity/propensity_crossvalidation.pdf",
        "results/figures/08_topology/calibration_precision_sweep.pdf",
        "results/figures/12_feature/topology_by_feature.pdf",
        "results/figures/13_rloop/rloop_fraction_by_topology.pdf",
        "results/figures/14_expression/promoter_g4_expression.pdf",
        # Stage G — topology metaprofile expansion
        "results/figures/10_metaprofiles/tss_signal_by_topology_norm.pdf",
        "results/figures/10_metaprofiles/feature_signal_by_topology.pdf",
        "results/figures/10_metaprofiles/deg_promoter_topology_grid.pdf",
        "results/figures/10_metaprofiles/deg_5utr_fragment_topology_grid.pdf",
        "results/figures/10_metaprofiles/feature_signal_by_topology_1utr.pdf",
        "results/figures/10_metaprofiles/feature_signal_total.pdf",
        # Stage H — unclassified peaks
        "results/figures/15_unclassified/unclassified_peak_metaprofile.pdf",
        "results/figures/15_unclassified/unclassified_peak_feature.pdf",
        "results/figures/15_unclassified/unclassified_peak_motifs.pdf",
        # Stage I — lncRNA
        "results/figures/17_lncrna/lncrna_g4_global_tss.pdf",
        # Stage J — G4 x R-loop
        "results/figures/18_g4_rloop_correlation/g4_vs_rloop_promoter_mean.pdf",
        "results/figures/19_g4_rloop_cooccurrence/g4_rloop_promoter_heatmap.pdf",
        # Stage K — all-PQS / multi-G4
        "results/tables/motif_topology.csv",
        "results/figures/20_motif_distributions/topology_composition_motif_vs_peak.pdf",
        "results/figures/10_metaprofiles/motif_signal_by_topology_tss.pdf",
        "results/figures/21_multiG4_expression/promoter_g4_count.pdf",
        # Stage L — strand + Gviz
        "results/tables/strand_bias_summary.tsv",
        "results/figures/22_strand/strand_bias_comprehensive.pdf",
        "results/figures/22_strand/metaprofile_by_strand.pdf",
        "results/figures/23_gviz/.done",


# ===========================================================================
# Stage A — RNA-seq
# ===========================================================================
rule rnaseq_qc:
    input: counts = COUNTS
    output:
        dds = "cache/dds.rds", rld = "cache/rld.rds",
        sinfo = "cache/sample_info.rds", gmeta = "cache/gene_meta.rds",
        fig = "results/figures/01_qc/pca_plot.pdf",
    log: "logs/01_rnaseq_qc.log"
    shell: "Rscript scripts/01_rnaseq_qc.R {CONFIG} > {log} 2>&1"

rule deseq2:
    input: dds = "cache/dds.rds", gmeta = "cache/gene_meta.rds"
    output:
        res = "cache/deseq2_res_list.rds",
        ko = "results/tables/deseq2_KO_vs_WT.csv",
        tabs = DESEQ_TABS,
    log: "logs/02_deseq2.log"
    shell: "Rscript scripts/02_deseq2.R {CONFIG} > {log} 2>&1"

rule deseq2_plots:
    input: res = "cache/deseq2_res_list.rds", rld = "cache/rld.rds", sinfo = "cache/sample_info.rds"
    output:
        union = "results/tables/significant_DEG_union.tsv",
        vol = VOLCANOES,
        heat = "results/figures/02_deseq2/union_DEG_heatmap.pdf",
        upset = "results/figures/02_deseq2/upset_DEGs.pdf",
    log: "logs/03_deseq2_plots.log"
    shell: "Rscript scripts/03_deseq2_plots.R {CONFIG} > {log} 2>&1"

rule enrichment:
    input: res = "cache/deseq2_res_list.rds"
    output: go = "results/tables/go_enrichment_summary.tsv"
    log: "logs/04_enrichment.log"
    shell: "Rscript scripts/04_enrichment.R {CONFIG} > {log} 2>&1"

rule gene_biotype:
    output:
        biotype = "results/tables/gene_biotype.csv",
        genes = "cache/gencode_genes.rds", lncrna = "cache/lncrna_genes.rds",
    log: "logs/05_gene_biotype.log"
    shell: "Rscript scripts/05_gene_biotype.R {CONFIG} > {log} 2>&1"

rule volcano_biotype:
    input: de = "results/tables/deseq2_KO_vs_WT.csv", biotype = "results/tables/gene_biotype.csv"
    output:
        labels = "results/tables/volcano_label_genes.csv",
        fig = "results/figures/16_volcano/volcano_DHX36KO_top50.pdf",
    log: "logs/06_volcano_biotype.log"
    shell: "Rscript scripts/06_volcano_biotype.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage B — peak calling + region building (REPLACES hand-built cache_V2)
# ===========================================================================
rule call_peaks:
    input: bw = G4_BW + RLOOP_BW
    output:
        peaks = PEAK_RDS,
        summary = "results/tables/peak_summary.csv",
    log: "logs/07_call_peaks.log"
    shell: "Rscript scripts/07_call_peaks.R {CONFIG} > {log} 2>&1"

rule build_regions:
    output:
        txdb = "cache/regions_txdb.rds",
        regs = expand("cache/regions_{x}.rds", x=["promoter","5UTR","intron1","genebody"]),
    log: "logs/08_build_regions.log"
    shell: "Rscript scripts/08_build_regions.R {CONFIG} > {log} 2>&1"

rule build_enhancers:
    input: regs = expand("cache/regions_{x}.rds", x=["promoter","5UTR","intron1","genebody"])
    output:
        enh = "cache/regions_enhancer.rds", all = "cache/regions_all.rds",
    log: "logs/09_build_enhancers.log"
    shell: "Rscript scripts/09_build_enhancers.R {CONFIG} > {log} 2>&1"

rule peak_qc:
    input: peaks = PEAK_RDS, regions = "cache/regions_all.rds"
    output:
        fig1 = "results/figures/04_peaks/01_peak_width_distribution.pdf",
        fig2 = "results/figures/04_peaks/01_regions_by_chromosome.pdf",
    log: "logs/10_peak_qc.log"
    shell: "Rscript scripts/10_peak_qc.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage C — genome-wide CUT&Tag QC
# ===========================================================================
rule genomewide_signal:
    input: bw = G4_BW + RLOOP_BW
    output:
        cache = "cache/genomewide_signal_mat.rds",
        fig1 = "results/figures/05_cuttag/violin_signal_G4_BG4.pdf",
        fig2 = "results/figures/05_cuttag/G4_vs_Rloop_correlation_panel.pdf",
    log: "logs/11_genomewide_signal.log"
    shell: "Rscript scripts/11_genomewide_signal.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage D — regional enrichment
# ===========================================================================
rule peak_region_enrichment:
    input: peaks = PEAK_RDS, regions = "cache/regions_all.rds"
    output:
        t1 = "results/tables/peak_counts_per_region.tsv",
        t2 = "results/tables/region_fold_enrichment.tsv",
        f1 = "results/figures/06_regional/02_peak_counts_stacked.pdf",
        f2 = "results/figures/06_regional/02_region_fold_enrichment.pdf",
    log: "logs/12_peak_region_enrichment.log"
    shell: "Rscript scripts/12_peak_region_enrichment.R {CONFIG} > {log} 2>&1"

rule region_signal:
    input: regions = "cache/regions_all.rds", bw = G4_BW + RLOOP_BW
    output:
        cache = "cache/region_signal_df.rds",
        stats = "results/tables/region_signal_stats.tsv",
        fig = "results/figures/06_regional/02_region_signal_boxplot.pdf",
    log: "logs/13_region_signal.log"
    shell: "Rscript scripts/13_region_signal.R {CONFIG} > {log} 2>&1"

rule region_metaprofiles:
    input: regions = "cache/regions_all.rds", bw = G4_BW + RLOOP_BW
    output:
        c1 = "cache/region_profile_df.rds", c2 = "cache/region_profile_full_df.rds",
        f1 = "results/figures/06_regional/02_region_metaprofile.pdf",
        f2 = "results/figures/06_regional/02_region_metaprofile_full_peaks.pdf",
    log: "logs/14_region_metaprofiles.log"
    shell: "Rscript scripts/14_region_metaprofiles.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage E — RNA x CUT&Tag integration
# ===========================================================================
rule integration_promoter_signal:
    input: res = "cache/deseq2_res_list.rds", gmeta = "cache/gene_meta.rds", bw = G4_BW + RLOOP_BW
    output: PROM_BOX
    log: "logs/15_integration_promoter_signal.log"
    shell: "Rscript scripts/15_integration_promoter_signal.R {CONFIG} > {log} 2>&1"

rule integration_tss_profiles:
    input: res = "cache/deseq2_res_list.rds", gmeta = "cache/gene_meta.rds", bw = G4_BW + RLOOP_BW
    output: TSS_COMB
    log: "logs/16_integration_tss_profiles.log"
    shell: "Rscript scripts/16_integration_tss_profiles.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage F — sequence-topology core (DeepG4 01-13, rewired to Stage-B cache)
# ===========================================================================
rule prepare_topology_inputs:
    input:
        peaks = expand("cache/peaks_G4_BG4_{g}.rds", g=GENO + ["union"]),
        promoters = "cache/regions_promoter.rds",
    output:
        peaks = "cache/peaks.rds", promoters = "cache/promoters.rds",
        catalog = "results/tables/peak_catalog.csv",
    log: "logs/17_prepare_topology_inputs.log"
    shell: "Rscript scripts/17_prepare_topology_inputs.R {CONFIG} > {log} 2>&1"

rule extract_sequences:
    input: peaks = "cache/peaks.rds"
    output:
        peak_fa = "cache/peak_seqs_201bp.fa", bg_fa = "cache/background_seqs.fa",
        pqs = "results/tables/peak_pqs.csv",
    log: "logs/18_extract_sequences.log"
    shell: "Rscript scripts/18_extract_sequences.R {CONFIG} > {log} 2>&1"

rule deepg4_predict:
    input: peak_fa = "cache/peak_seqs_201bp.fa", bg_fa = "cache/background_seqs.fa"
    output: scores = "results/tables/deepg4_scores.csv"
    log: "logs/19_deepg4_predict.log"
    shell: "Rscript scripts/19_deepg4_predict.R {CONFIG} > {log} 2>&1"

rule g4sp_topology:
    input: pqs = "results/tables/peak_pqs.csv"
    output: topo = "results/tables/g4sp_topology.csv"
    log: "logs/20_g4sp_topology.log"
    shell:
        "conda run -n {G4SP_ENV} python scripts/20_g4sp_topology.py "
        "--config {CONFIG} --pqs {input.pqs} --out {output.topo} > {log} 2>&1"

rule assign_topology:
    input: pqs = "results/tables/peak_pqs.csv", g4sp = "results/tables/g4sp_topology.csv"
    output:
        topo = "results/tables/peak_topology.csv",
        conc = "results/tables/topology_concordance.csv",
        fig = "results/figures/08_topology/topology_counts.pdf",
    log: "logs/21_assign_topology.log"
    shell: "Rscript scripts/21_assign_topology.R {CONFIG} > {log} 2>&1"

rule sanity_check:
    input:
        scores = "results/tables/deepg4_scores.csv",
        catalog = "results/tables/peak_catalog.csv",
        topo = "results/tables/peak_topology.csv",
    output:
        metrics = "results/tables/sanity_metrics.csv",
        fig = "results/figures/08_topology/deepg4_roc.pdf",
    log: "logs/22_sanity_check.log"
    shell: "Rscript scripts/22_sanity_check.R {CONFIG} > {log} 2>&1"

rule gain_loss:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv", bw = G4_BW
    output:
        comp = "results/tables/topology_composition.csv",
        diff = "results/tables/topology_differential.csv",
        diff_sum = "results/tables/topology_differential_summary.csv",
        enrich = "results/tables/topology_lost_enrichment.csv",
        fig = "results/figures/09_gain_loss/differential_class_fraction.pdf",
        fig2 = "results/figures/09_gain_loss/lost_enrichment_by_topology.pdf",
    log: "logs/23_gain_loss.log"
    shell: "Rscript scripts/23_gain_loss.R {CONFIG} > {log} 2>&1"

rule tss_metaprofiles:
    input: promoters = "cache/promoters.rds", peaks = "cache/peaks.rds",
           topo = "results/tables/peak_topology.csv", bw = G4_BW
    output:
        profiles = "results/tables/tss_topology_profiles.csv",
        profiles_norm = "results/tables/tss_topology_profiles_norm.csv",
        fig = "results/figures/10_metaprofiles/tss_signal_by_topology.pdf",
        fig_norm = "results/figures/10_metaprofiles/tss_signal_by_topology_norm.pdf",
    log: "logs/24_tss_metaprofiles.log"
    shell: "Rscript scripts/24_tss_metaprofiles.R {CONFIG} > {log} 2>&1"

rule propensity_metrics:
    input: scores = "results/tables/deepg4_scores.csv", pqs = "results/tables/peak_pqs.csv",
           topo = "results/tables/peak_topology.csv", diff = "results/tables/topology_differential.csv",
           peak_fa = "cache/peak_seqs_201bp.fa"
    output:
        prop = "results/tables/propensity_metrics.csv",
        cv = "results/tables/propensity_crossvalidation.csv",
        dep = "results/tables/propensity_by_dependence.csv",
        fig = "results/figures/11_propensity/propensity_crossvalidation.pdf",
    log: "logs/25_propensity_metrics.log"
    shell: "Rscript scripts/25_propensity_metrics.R {CONFIG} > {log} 2>&1"

rule threshold_sweep:
    input: g4sp = "results/tables/g4sp_topology.csv"
    output: sweep = "results/tables/topology_precision_sweep.csv"
    log: "logs/26a_threshold_sweep.log"
    shell:
        "conda run -n {G4SP_ENV} python scripts/26a_threshold_sweep.py "
        "--config {CONFIG} --g4sp {input.g4sp} --out {output.sweep} > {log} 2>&1"

rule control_benchmark:
    input: controls = config["calibration"]["controls_tsv"]
    output: bench = "results/tables/g4sp_control_benchmark.csv"
    log: "logs/26b_control_benchmark.log"
    shell:
        "conda run -n {G4SP_ENV} python scripts/26b_control_benchmark.py "
        "--config {CONFIG} --out {output.bench} > {log} 2>&1"

rule calibration_plots:
    input: sweep = "results/tables/topology_precision_sweep.csv",
           g4sp = "results/tables/g4sp_topology.csv",
           bench = "results/tables/g4sp_control_benchmark.csv"
    output: fig = "results/figures/08_topology/calibration_precision_sweep.pdf"
    log: "logs/26c_calibration_plots.log"
    shell: "Rscript scripts/26c_calibration_plots.R {CONFIG} > {log} 2>&1"

rule feature_partition:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           prop = "results/tables/propensity_metrics.csv"
    output:
        by_feature = "results/tables/topology_by_feature.csv",
        tests = "results/tables/topology_feature_tests.csv",
        prop_feat = "results/tables/propensity_by_feature.csv",
        fig = "results/figures/12_feature/topology_by_feature.pdf",
    log: "logs/27_feature_partition.log"
    shell: "Rscript scripts/27_feature_partition.R {CONFIG} > {log} 2>&1"

rule rloop_coupling:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           prop = "results/tables/propensity_metrics.csv", bw = RLOOP_BW
    output:
        overlap = "results/tables/g4_rloop_overlap.csv",
        signal = "results/tables/rloop_signal_by_topology.csv",
        diff = "results/tables/rloop_differential_at_g4.csv",
        fig = "results/figures/13_rloop/rloop_fraction_by_topology.pdf",
    log: "logs/28_rloop_coupling.log"
    shell: "Rscript scripts/28_rloop_coupling.R {CONFIG} > {log} 2>&1"

rule promoter_g4_expression:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           de = "results/tables/deseq2_KO_vs_WT.csv"
    output:
        expr = "results/tables/promoter_g4_expression.csv",
        tests = "results/tables/promoter_g4_expression_tests.csv",
        fig = "results/figures/14_expression/promoter_g4_expression.pdf",
    log: "logs/29_promoter_g4_expression.log"
    shell: "Rscript scripts/29_promoter_g4_expression.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage G — topology metaprofile expansion (DeepG4 14-17, 24)
# ===========================================================================
rule feature_metaprofiles:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv", bw = G4_BW
    output:
        profiles = "results/tables/feature_topology_profiles.csv",
        fig = "results/figures/10_metaprofiles/feature_signal_by_topology.pdf",
        fig_fixed = "results/figures/10_metaprofiles/feature_signal_by_topology_fixedY.pdf",
    log: "logs/30_feature_metaprofiles.log"
    shell: "Rscript scripts/30_feature_metaprofiles.R {CONFIG} > {log} 2>&1"

rule deg_promoter_metaprofiles:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           de = "results/tables/deseq2_KO_vs_WT.csv", bw = G4_BW
    output:
        profiles = "results/tables/deg_promoter_profiles.csv",
        fig_grid = "results/figures/10_metaprofiles/deg_promoter_topology_grid.pdf",
        fig_pool = "results/figures/10_metaprofiles/deg_promoter_pooled.pdf",
        fig_pool_unf = "results/figures/10_metaprofiles/deg_promoter_pooled_unfiltered.pdf",
        csv_counts = "results/tables/deg_promoter_counts.csv",
        fig_counts = "results/figures/10_metaprofiles/deg_promoter_counts.pdf",
    log: "logs/31_deg_promoter_metaprofiles.log"
    shell: "Rscript scripts/31_deg_promoter_metaprofiles.R {CONFIG} > {log} 2>&1"

rule deg_5utr_metaprofiles:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           de = "results/tables/deseq2_KO_vs_WT.csv", bw = G4_BW
    output:
        profiles = "results/tables/deg_5utr_profiles.csv",
        frag_grid = "results/figures/10_metaprofiles/deg_5utr_fragment_topology_grid.pdf",
        frag_pool = "results/figures/10_metaprofiles/deg_5utr_fragment_pooled.pdf",
        frag_pool_unf = "results/figures/10_metaprofiles/deg_5utr_fragment_pooled_unfiltered.pdf",
        gene_grid = "results/figures/10_metaprofiles/deg_5utr_gene_topology_grid.pdf",
        gene_pool = "results/figures/10_metaprofiles/deg_5utr_gene_pooled.pdf",
        gene_pool_unf = "results/figures/10_metaprofiles/deg_5utr_gene_pooled_unfiltered.pdf",
        csv_counts = "results/tables/deg_5utr_counts.csv",
        fig_counts = "results/figures/10_metaprofiles/deg_5utr_counts.pdf",
    log: "logs/32_deg_5utr_metaprofiles.log"
    shell: "Rscript scripts/32_deg_5utr_metaprofiles.R {CONFIG} > {log} 2>&1"

rule feature_metaprofiles_1utr:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           feat14 = "results/tables/feature_topology_profiles.csv", bw = G4_BW
    output:
        profiles = "results/tables/feature_topology_profiles_1utr.csv",
        fig = "results/figures/10_metaprofiles/feature_signal_by_topology_1utr.pdf",
        fig_fixed = "results/figures/10_metaprofiles/feature_signal_by_topology_1utr_fixedY.pdf",
    log: "logs/33_feature_metaprofiles_1utr.log"
    shell: "Rscript scripts/33_feature_metaprofiles_1utr.R {CONFIG} > {log} 2>&1"

rule feature_total_metaprofiles:
    input: feat = expand("cache/regions_{x}.rds", x=["promoter","5UTR","intron1","enhancer"]),
           lncrna = "cache/lncrna_genes.rds", bw = G4_BW
    output:
        profiles = "results/tables/feature_total_profiles.csv",
        fig = "results/figures/10_metaprofiles/feature_signal_total.pdf",
        fig_fixed = "results/figures/10_metaprofiles/feature_signal_total_fixedY.pdf",
    log: "logs/34_feature_total_metaprofiles.log"
    shell: "Rscript scripts/34_feature_total_metaprofiles.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage H — unclassified (no-canonical-PQS) peaks (DeepG4 18-20)
# ===========================================================================
rule unclassified_peak_profiles:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           prop = "results/tables/propensity_metrics.csv", deepg4 = "results/tables/deepg4_scores.csv",
           pk_fa = "cache/peak_seqs_201bp.fa", bg_fa = "cache/background_seqs.fa", bw = G4_BW
    output:
        profiles = "results/tables/unclassified_peak_profiles.csv",
        seqstats = "results/tables/unclassified_peak_seqstats.csv",
        fig_meta = "results/figures/15_unclassified/unclassified_peak_metaprofile.pdf",
        fig_prop = "results/figures/15_unclassified/unclassified_peak_propensity.pdf",
        fig_comp = "results/figures/15_unclassified/unclassified_peak_composition.pdf",
    log: "logs/35_unclassified_peak_profiles.log"
    shell: "Rscript scripts/35_unclassified_peak_profiles.R {CONFIG} > {log} 2>&1"

rule unclassified_peak_context:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           pk_fa = "cache/peak_seqs_201bp.fa", bg_fa = "cache/background_seqs.fa"
    output:
        feat = "results/tables/unclassified_feature_enrichment.csv",
        rep = "results/tables/unclassified_repeat_enrichment.csv",
        lowc = "results/tables/unclassified_lowcomplexity.csv",
        fig_ft = "results/figures/15_unclassified/unclassified_peak_feature.pdf",
        fig_rp = "results/figures/15_unclassified/unclassified_peak_repeat.pdf",
        fig_lc = "results/figures/15_unclassified/unclassified_peak_lowcomplexity.pdf",
    log: "logs/36_unclassified_peak_context.log"
    shell: "Rscript scripts/36_unclassified_peak_context.R {CONFIG} > {log} 2>&1"

rule unclassified_peak_motifs:
    input: peaks = "cache/peaks.rds", topo = "results/tables/peak_topology.csv",
           pk_fa = "cache/peak_seqs_201bp.fa", bg_fa = "cache/background_seqs.fa"
    output:
        motifs = "results/tables/unclassified_motif_enrichment.csv",
        fig = "results/figures/15_unclassified/unclassified_peak_motifs.pdf",
    log: "logs/37_unclassified_peak_motifs.log"
    shell: "Rscript scripts/37_unclassified_peak_motifs.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage I — lncRNA G4 metaprofiles (DeepG4 23)
# ===========================================================================
rule lncrna_g4_metaprofiles:
    input: lncrna = "cache/lncrna_genes.rds", peaks = "cache/peaks.rds",
           topo = "results/tables/peak_topology.csv", diff = "results/tables/topology_differential.csv",
           de = "results/tables/deseq2_KO_vs_WT.csv", bw = G4_BW
    output:
        gl_prof = "results/tables/lncrna_g4_gainloss_profiles.csv",
        glob_prof = "results/tables/lncrna_g4_global_profiles.csv",
        counts = "results/tables/lncrna_g4_counts.csv",
        fig_glob = "results/figures/17_lncrna/lncrna_g4_global_tss.pdf",
        fig_gl = "results/figures/17_lncrna/lncrna_g4_gainloss.pdf",
        fig_deg = "results/figures/17_lncrna/lncrna_g4_deg_pooled.pdf",
        fig_cnt = "results/figures/17_lncrna/lncrna_g4_counts.pdf",
    log: "logs/38_lncrna_g4_metaprofiles.log"
    shell: "Rscript scripts/38_lncrna_g4_metaprofiles.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage J — G4 x R-loop (DeepG4 25-26)
# ===========================================================================
rule g4_rloop_promoter_correlation:
    input: promoters = "cache/regions_promoter.rds", g4_bw = G4_BW, rloop_bw = RLOOP_BW
    output:
        per_prom = "results/tables/g4_rloop_promoter_signal.csv",
        stats = "results/tables/g4_rloop_promoter_correlation_stats.csv",
        fig_mean = "results/figures/18_g4_rloop_correlation/g4_vs_rloop_promoter_mean.pdf",
        fig_bin = "results/figures/18_g4_rloop_correlation/g4_vs_rloop_promoter_binned.pdf",
        fig_hist = "results/figures/18_g4_rloop_correlation/g4_vs_rloop_binned_per_promoter_hist.pdf",
        fig_summ = "results/figures/18_g4_rloop_correlation/g4_vs_rloop_correlation_summary.pdf",
    log: "logs/39_g4_rloop_promoter_correlation.log"
    shell: "Rscript scripts/39_g4_rloop_promoter_correlation.R {CONFIG} > {log} 2>&1"

rule g4_rloop_cooccurrence:
    input: promoters = "cache/regions_promoter.rds", peaks = "cache/peaks.rds",
           g4_bw = G4_BW, rloop_bw = RLOOP_BW
    output:
        strat_tbl = "results/tables/g4_rloop_stratified_profiles.csv",
        violin_tbl = "results/tables/g4_rloop_violin_stats.csv",
        ko_tbl = "results/tables/g4_rloop_ko_differential.csv",
        fig_heat = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_promoter_heatmap.pdf",
        fig_strat_r = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_stratified_rloopByG4.pdf",
        fig_strat_g = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_stratified_g4ByRloop.pdf",
        fig_violin = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_promoter_violin.pdf",
        fig_delta = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_delta_scatter.pdf",
        fig_volcano = "results/figures/19_g4_rloop_cooccurrence/g4_rloop_rloop_volcano.pdf",
    log: "logs/40_g4_rloop_cooccurrence.log"
    shell: "Rscript scripts/40_g4_rloop_cooccurrence.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage K — all-PQS / multi-G4 motif analysis (DeepG4 27-32)
# ===========================================================================
rule extract_all_pqs:
    input: peaks = "cache/peaks.rds"
    output:
        motifs = "cache/motifs_all_base.rds",
        catalog = "results/tables/motif_catalog.csv",
        g4sp_in = "results/tables/motif_pqs_for_g4sp.csv",
    log: "logs/41_extract_all_pqs.log"
    shell: "Rscript scripts/41_extract_all_pqs.R {CONFIG} > {log} 2>&1"

rule g4sp_topology_all:
    input: pqs = "results/tables/motif_pqs_for_g4sp.csv"
    output: topo = "results/tables/motif_g4sp_topology.csv"
    log: "logs/42_g4sp_topology_all.log"
    shell:
        "conda run -n {G4SP_ENV} python scripts/20_g4sp_topology.py "
        "--config {CONFIG} --pqs {input.pqs} --out {output.topo} > {log} 2>&1"

rule assign_motif_topology:
    input: catalog = "results/tables/motif_catalog.csv",
           g4sp = "results/tables/motif_g4sp_topology.csv", motifs = "cache/motifs_all_base.rds"
    output:
        topo = "results/tables/motif_topology.csv", motifs = "cache/motifs_all.rds",
    log: "logs/43_assign_motif_topology.log"
    shell: "Rscript scripts/43_assign_motif_topology.R {CONFIG} > {log} 2>&1"

rule peak_motif_distributions:
    input: peaks = "cache/peaks.rds", motifs = "cache/motifs_all.rds",
           mtopo = "results/tables/motif_topology.csv", ptopo = "results/tables/peak_topology.csv"
    output:
        summ = "results/tables/peak_motif_summary.csv",
        comp = "results/tables/topology_composition_compare.csv",
        fig1 = "results/figures/20_motif_distributions/peak_width_distribution.pdf",
        fig2 = "results/figures/20_motif_distributions/motifs_per_peak.pdf",
        fig3 = "results/figures/20_motif_distributions/topology_composition_motif_vs_peak.pdf",
    log: "logs/44_peak_motif_distributions.log"
    shell: "Rscript scripts/44_peak_motif_distributions.R {CONFIG} > {log} 2>&1"

rule promoter_motif_metaprofiles:
    input: promoters = "cache/regions_promoter.rds", motifs = "cache/motifs_all.rds", bw = G4_BW
    output:
        profiles = "results/tables/promoter_motif_profiles.csv",
        signal = "cache/promoter_motif_signal.rds",
        fig_sig = "results/figures/10_metaprofiles/motif_signal_by_topology_tss.pdf",
        fig_den = "results/figures/10_metaprofiles/motif_density_by_topology_tss.pdf",
    log: "logs/45_promoter_motif_metaprofiles.log"
    shell: "Rscript scripts/45_promoter_motif_metaprofiles.R {CONFIG} > {log} 2>&1"

rule deg_promoter_motif_metaprofiles:
    input: signal = "cache/promoter_motif_signal.rds", de = "results/tables/deseq2_KO_vs_WT.csv"
    output:
        profiles = "results/tables/deg_promoter_motif_profiles.csv",
        fig_grid = "results/figures/10_metaprofiles/deg_motif_topology_grid.pdf",
        fig_pool = "results/figures/10_metaprofiles/deg_motif_pooled.pdf",
        fig_den = "results/figures/10_metaprofiles/deg_motif_density.pdf",
    log: "logs/46_deg_promoter_motif_metaprofiles.log"
    shell: "Rscript scripts/46_deg_promoter_motif_metaprofiles.R {CONFIG} > {log} 2>&1"

rule multiG4_promoter_expression:
    input: promoters = "cache/regions_promoter.rds", motifs = "cache/motifs_all.rds",
           de = "results/tables/deseq2_KO_vs_WT.csv"
    output:
        counts = "results/tables/multiG4_promoter_counts.csv",
        stats = "results/tables/multiG4_expression_stats.csv",
        fig1 = "results/figures/21_multiG4_expression/promoter_g4_count.pdf",
        fig2 = "results/figures/21_multiG4_expression/g4count_vs_baseline_expression.pdf",
        fig3 = "results/figures/21_multiG4_expression/g4count_vs_ko_response.pdf",
        fig4 = "results/figures/21_multiG4_expression/g4count_topology_composition.pdf",
    log: "logs/47_multiG4_promoter_expression.log"
    shell: "Rscript scripts/47_multiG4_promoter_expression.R {CONFIG} > {log} 2>&1"


# ===========================================================================
# Stage L — strand asymmetry + Gviz locus tracks
# ===========================================================================
rule strand_asymmetry:
    input: peaks = expand("cache/peaks_G4_BG4_{g}.rds", g=GENO),
           promoter = "cache/regions_promoter.rds", utr5 = "cache/regions_5UTR.rds", bw = G4_BW
    output:
        scored = "cache/sequences_scored.rds",
        t1 = "results/tables/strand_bias_summary.tsv",
        t2 = "results/tables/strand_bias_relative.tsv",
        t3 = "results/tables/g4_detection_by_strand.tsv",
        t4 = "results/tables/g4_score_stats.tsv",
        f1 = "results/figures/22_strand/strand_bias_comprehensive.pdf",
        f2 = "results/figures/22_strand/metaprofile_by_strand.pdf",
        f3 = "results/figures/22_strand/mean_signal_by_strand.pdf",
        f4 = "results/figures/22_strand/metaprofile_template_vs_coding.pdf",
    log: "logs/48_strand_asymmetry.log"
    shell: "Rscript scripts/48_strand_asymmetry.R {CONFIG} > {log} 2>&1"

rule gviz_locus_tracks:
    input: scored = "cache/sequences_scored.rds", bw = G4_BW
    output: touch("results/figures/23_gviz/.done")
    log: "logs/49_gviz_locus_tracks.log"
    shell: "Rscript scripts/49_gviz_locus_tracks.R {CONFIG} > {log} 2>&1"
