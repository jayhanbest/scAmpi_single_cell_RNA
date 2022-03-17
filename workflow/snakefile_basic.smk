from snakemake.utils import min_version
from snakemake.utils import validate
import pandas as pd
import os

wildcard_constraints:
    sample = "[^/]+",
    i = "[^/]+"

# Define minimum Snakemake version
min_version("6.12.1")

# Include Config file
configfile: "config/config.yaml"

# Include report functionality
#report: "../report/workflow.rst"

# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"

# Include rules
include: "rules/scAmpi_basic_rules.smk"
include: "rules/scAmpi_clinical_rules.smk"
include: "rules/velocity_rules_snake.smk"

#outdir = checkpoints.diff_exp_analysis.get(sample = sample_ids).output[0]
#print(outdir)

# This rule defines which files should be created
localrules: scAmpi_basic
rule scAmpi_basic:
    input:
        expand('results/cellranger_run/{sample}.features.tsv', sample = sample_ids),
        expand('results/counts_raw/{sample}.h5', sample = sample_ids),
        expand('results/counts_filtered/{sample}.doublet_barcodes.txt', sample = sample_ids),
        expand('results/counts_raw/{sample}.h5.histogram_library_sizes.png', sample = sample_ids),
        expand('results/counts_filtered/{sample}.genes_cells_filtered.h5.histogram_library_sizes.png', sample = sample_ids),
        expand('results/counts_corrected/{sample}.corrected.RDS', sample = sample_ids),
        expand('results/clustering/{sample}.clusters_phenograph.csv', sample = sample_ids),
        expand('results/atypical_removed/{sample}.atypical_removed.RDS', sample = sample_ids),
        expand('results/clustering/{sample}.clusters_cell_count_percent.txt', sample = sample_ids),
        expand('results/gene_exp/{sample}.gene_expression_per_cluster.tsv', sample = sample_ids),
        expand('results/plotting/{sample}.celltype_barplot.png', sample = sample_ids),
        expand('results/gsva/{sample}.gsetscore_hm.png', sample = sample_ids),
        expand('results/diff_exp_analysis/{sample}/', sample = sample_ids),
#        # trigger clinical part of the pipeline
#        expand('results/aggregated/{sample}.aggregated.txt', sample = sample_ids),
#        # plot_drug_prediction is also aggregation rule (as is aggregate)
#        expand('results/plot_drug_prediction/{sample}.drug_prediction_umap.png', sample = sample_ids),
#        # plot gene set enrichment heatmap (is also aggregation rule)
#        expand('results/gene_set_enrichment/{sample}.heatmap_enrichment.png', sample = sample_ids),
#        expand('results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.heatmap_enrichment.png', sample = sample_ids),
#        # parse_for_minSetCover (is also aggregation rule)
#        expand('results/drug_combination/{sample}.drugToCluster.allDrugs.txt', sample = sample_ids),
#        expand('results/drug_combination/{sample}.drugToCluster.filteredDrugs.txt', sample = sample_ids),
#        # preprocess for upsetR plot
#        expand('results/upsetr_plot/{sample}.drugToCluster.allDrugs.processedForUpSetR.txt', sample = sample_ids),
#        expand('results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.processedForUpSetR.txt', sample = sample_ids),
#        # plot upset plot
#        expand('results/upsetr_plot/{sample}.drugToCluster.allDrugs.vennplot.png', sample = sample_ids),
#        expand('results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.vennplot.png', sample = sample_ids),
#        # find minSetCover
#        expand('results/drug_combination/{sample}.drugCombination.allDrugs.txt', sample = sample_ids),
#        expand('results/drug_combination/{sample}.drugCombination.filteredDrugs.txt', sample = sample_ids),
#        # druglist_to_subclones
#        expand('results/drug_combination/{sample}.full_druglist_to_subclones.txt', sample = sample_ids),
#        # velocity analysis (2 steps)
        expand('results/velocity/velocyto/{sample}/{sample}.loom', sample = sample_ids),
#        expand('results/velocity/{sample}/{sample}.velocity_success.txt', sample = sample_ids),
    output:
        'results/complete_scAmpi_basic.txt'
    params:
        lsfoutfile = 'results/complete_scAmpi_basic.lsfout.log',
        lsferrfile = 'results/complete_scAmpi_basic.lsferr.log',
        mem = '1000',
        scratch = '1000',
        time = '1'
    benchmark:
        'results/complete_scAmpi_basic.txt.benchmark'
    shell:
        'date > {output}'

