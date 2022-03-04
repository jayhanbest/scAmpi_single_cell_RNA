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

#config["inputOutput"]["sample_map"]

# old framework
#SAMPLEMAPPING = config['inputOutput']['sample_map']


# Include the rules
include: "rules/scAmpi_basic_rules.smk"

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
#        expand('results/diff_exp_analysis/{sample}/{sample}.diff_exp_analysis_success.txt', sample = sample_ids),
        expand('results/diff_exp_analysis/{sample}/', sample = sample_ids),
#        expand('results/databaseQuery/{sample}.{clusterid}.dgidb.txt', sample = sample_ids, clusterid = cluster_id),
#        expand('results/intermediate/{sample}.txt', sample = sample_ids),
        expand('results/query_dgidb/{sample}.aggregated.txt', sample = sample_ids),
#        expand('results/parse_diff_exp/{sample}.{i}.txt', sample = sample_ids, i = get_cluster_ids),
#        expand('results/intermediate/{sample}/{i}.txt', sample = sample_ids),
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


#def check_checkpoint(wildcards):
#    checkpoints.diff_exp_analysis.get(**wildcards).output[0]
#print(TEST)
