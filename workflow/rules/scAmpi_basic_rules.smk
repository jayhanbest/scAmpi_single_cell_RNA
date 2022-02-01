# cellranger call to process the raw samples
rule cellranger_count: 
    input:
        fastqs_dir = config['inputOutput']['input_fastqs'],
        reference = config['resources']['reference_transcriptome']
    output:
        features_file = 'results/cellranger_run/{sample}.features.tsv',
        matrix_file = 'results/cellranger_run/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_run/{sample}.barcodes.tsv'
    params:
        cr_out = 'results/cellranger_run/',
        local_cores = config['tools']['cellranger_count']['local_cores'],
        variousParams = config['tools']['cellranger_count']['variousParams'],
        metrics_summary = 'results/cellranger_run/{sample}.metrics_summary.csv',
        web_summary = 'results/cellranger_run/{sample}.web_summary.html',
	mySample = '{sample}' # needs to be the prefix of all fastq files that belong to this sample. NOTE: no dots are allowed in sample names!
    resources:
        mem_mb = config['computingResources']['highRequirements']['mem'],
        time_min = config['computingResources']['highRequirements']['time']
    threads:
        config['computingResources']['highRequirements']['threads']
    benchmark:
        'results/cellranger_run/{sample}.cellranger_count.benchmark'
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    shell:
	    '(cd {params.cr_out}; {config[tools][cellranger_count][call]} count --id={params.mySample} --sample={params.mySample} --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary {params.variousParams}); gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{params.mySample}/outs/web_summary.html" "{params.web_summary}" ; ln -s "{params.cr_out}{params.mySample}/outs/metrics_summary.csv" "{params.metrics_summary}"'

# create hdf5 from count matrix, genes, and cell barcodes file
rule create_hdf5:
    input:
        genes_file = 'results/cellranger_run/{sample}.features.tsv',
        matrix_file = 'results/cellranger_run/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_run/{sample}.barcodes.tsv'
    output:
        outfile = 'results/rawCounts/{sample}.h5'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/rawCounts/{sample}.create_hd5.benchmark'
    shell:
        'python ../scripts/create_hdf5.py -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output.outfile}'

# identify doublets with scDblFinder
rule identify_doublets:
    input:
        infile = 'results/rawCounts/{sample}.h5'
    output:
        outfile = 'results/filteredCounts/{sample}.doublet_barcodes.txt'
    params:
        sample = '{sample}',
        outdir = 'results/filteredCounts/',
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/filteredCounts/{sample}.identify_doublets.benchmark'
    shell:
        "Rscript ../scripts/identify_doublets.R " +
        "--hdf5File {input.infile} " +
        "--sample {params.sample} " +
        "--outdir {params.outdir}"


# filter out whole genes and cells.
# Filter out cells if the fraction of reads mapped to MT genes is too high or the NODG is too low.
# Filter out genes if they are not protein coding, if they are mitochondrial genes, or if they encode for ribosomal proteins.
rule filter_genes_and_cells:
    input:
        infile = 'results/rawCounts/{sample}.h5',
        doublets = 'results/filteredCounts/{sample}.doublet_barcodes.txt'
    output:
        outfile = 'results/filteredCounts/{sample}.genes_cells_filtered.h5'
    params:
        nmads_fractionMT = config['tools']['filter_genes_and_cells']['nmads_fractionMT'],
        nmads_NODG = config['tools']['filter_genes_and_cells']['nmads_NODG'],
        threshold_fractionMT = config['tools']['filter_genes_and_cells']['threshold_fractionMT'],
        threshold_NODG = config['tools']['filter_genes_and_cells']['threshold_NODG'],
        remove_doublets = config['tools']['filter_genes_and_cells']['remove_doublets'],
        outDir = 'results/filteredCounts/',
        genomeVersion = config['tools']['filter_genes_and_cells']['genomeVersion'],
        sample = '{sample}'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/filteredCounts/{sample}.filter_genes_and_cells.benchmark'
    shell:
            'Rscript ../scripts/filter_genes_and_cells.R ' +
            '--hdf5File {input.infile} ' +
            '--nmads_NODG {params.nmads_NODG} ' +
            '--nmads_fractionMT {params.nmads_fractionMT} ' +
            '--threshold_NODG {params.threshold_NODG} ' +
            '--threshold_fractionMT {params.threshold_fractionMT} ' +
            '--genomeVersion {params.genomeVersion} ' +
            '--doublet_barcodes {input.doublets} ' +
            '--remove_doublets {params.remove_doublets} ' +
            '--sample {params.sample} ' +
            '--outDir {params.outDir}'



# perform normalisation, cell cycle correction and other preprocessing using sctransform
rule sctransform_preprocessing:
    input:
        hdf5_file =  'results/filteredCounts/{sample}.h5',
    output:
        outfile = 'results/counts_corrected/{sample}.corrected.RDS',
        highly_variable = 'results/counts_corrected/{sample}.corrected.variable_genes.h5',
    params:
        sample = '{sample}',
        number_genes = config['tools']['sctransform_preprocessing']['number_genes'],
        min_var = config['tools']['sctransform_preprocessing']['min_var'],
        n_nn = config['tools']['sctransform_preprocessing']['n_nn'],
        outDir = 'results/counts_corrected/',
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/counts_corrected/{sample}.corrected.benchmark'
    shell:
        "Rscript ../scripts/sctransform_preprocessing.R --inHDF5 {input.hdf5_file} --sample {params.sample} --number_genes {params.number_genes} --min_var {params.min_var} --n_nn {params.n_nn} --outdir {params.outDir} "


# perform clustering with phenograph
rule phenograph:
    input:
        infile = 'results/counts_corrected/{sample}.variable_genes.h5'
    output:
        outfile = 'results/clustering/{sample}.clusters_phenograph.csv',
        distance_matrix = 'results/clustering/{sample}.distance_matrix.tsv',
        modularity_score = 'results/clustering/{sample}.modularity_score.txt'
    params:
        n_neighbours = config['tools']['clustering']['phenograph']['n_neighbours'],
        min_cluster_size = config['tools']['clustering']['phenograph']['min_cluster_size'],
        log_normalize = config['tools']['clustering']['phenograph']['log_normalize'],
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/clustering/{sample}.phenograph.benchmark'
    shell:
        'python ../scripts/apply_phenograph.py ' +
        '--input_file {input.infile} ' +
        '--output_file {output.outfile} ' +
        '--distance_matrix {output.distance_matrix} ' +
        '--modularity_score {output.modularity_score} ' +
        '--n_neighbours {params.n_neighbours} ' +
        '--min_size {params.min_cluster_size} ' +
        '--log_normalize {params.log_normalize} ' +
        '--n_threads {threads}'


# prepare sce object in RDS file for cell type classification
rule prepare_celltyping:
    input:
        RDS_file = 'results/counts_corrected/{sample}.RDS',
        cluster = 'results/clustering/{sample}.clusters_phenograph.csv',
        distanceMatrix = 'results/clustering/{sample}.distance_matrix.tsv',
        modularity_score = 'results/clustering/{sample}.modularity_score.txt',
    output:
        outfile = 'results/prep_celltyping/{sample}.RDS'
    params:
        outputDirec = 'results/prep_celltyping/',
        sampleName = '{sample}',
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/prep_celltyping/{sample}.prepare_celltyping.benchmark'
    shell:
        "Rscript ../scripts/prepare_celltyping.R --in_sce {input.RDS_file} --phenograph_cluster {input.cluster} --outputDirec {params.outputDirec} --sampleName {params.sampleName} --distanceMatrix {input.distanceMatrix} --modularity_score {input.modularity_score} "


# perform cell type classification
rule cell_type_classification:
    input:
        infile = 'results/prep_celltyping/{sample}.RDS',
    output:
        outfile = 'results/celltype_classification/{sample}.phenograph_celltype_association.txt',
        out_sce = 'results/celltype_classification/{sample}.RDS'
    params:
        min_genes = config['tools']['cell_type_classification']['min_genes'],
        celltype_lists = config['resources']['celltype_lists'],
        celltype_config = config['resources']['celltype_config'],
        outputDirec = 'results/celltype_classification/',
        sampleName = '{sample}',
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/celltype_classification/{sample}.cell_type_classification.benchmark'
    shell:
        'Rscript ../scripts/celltyping.r ' +
        '--SCE {input.infile} ' +
        '--celltype_lists {params.celltype_lists} ' +
        '--celltype_config {params.celltype_config} ' +
        '--sampleName {params.sampleName} ' +
        '--min_genes {params.min_genes} ' +
        '--outputDirec {params.outputDirec} '


# filter out atypical cells from sce object
rule remove_atypical:
    input:
        infile = 'results/celltype_classification/{sample}.RDS',
        cluster_table = 'results/celltype_classification/{sample}.phenograph_celltype_association.txt',
    output:
        out_sce = 'results/atypical_removed/{sample}.atypical_removed.RDS',
        out_table = 'results/atypical_removed/{sample}.atypical_removed.phenograph_celltype_association.txt'
    params:
        celltype_config = config['resources']['celltype_config'],
        outputDirec = 'results/atypical_removed/',
        sample_name = '{sample}',
        threshold_filter = config['tools']['remove_atypical']['threshold_filter'],
        min_threshold = config['tools']['remove_atypical']['min_threshold'],
        threshold_type = config['tools']['remove_atypical']['threshold_type'],
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/atypical_removed/{sample}.atypical_removed.benchmark'
    shell:
        "Rscripts ../scripts/filter_out_atypical_cells.R --sce_in {input.infile} --cluster_table {input.cluster_table} --celltype_config {params.celltype_config} --threshold_filter {params.threshold_filter} --min_threshold {params.min_threshold} --threshold_type {params.threshold_type} --outDir {params.outputDirec} --sample_name {params.sample_name} "


# perform gsva gene set analysis
rule gsva:
    input:
        infile = 'results/atypical_removed/{sample}.RDS',
    output:
        outfile = 'results/gsva/{sample}.gsetscore_hm.png',
    params:
        outputDirec = 'results/gsva/',
        sampleName = '{sample}',
        genesets = config['resources']['genesets'],
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gsva/{sample}.gsva.benchmark'
    shell:
        "Rscript ../scripts/gsva.r --SCE {input.infile} --geneset {params.genesets} --outputDirec {params.outputDirec} --sampleName {params.sampleName} "


# generate plots about sample composition and gene expression
rule plotting:
    input:
        infile = 'results/atypical_removed/{sample}.RDS',
    output:
        outfile = 'results/plotting/{sample}.celltype_barplot.png',
    params:
        outputDirec = 'results/plotting/',
        sampleName = '{sample}',
        genes_of_interest = config['resources']['priority_genes'],
        colour_config = config['resources']['colour_config'],
        use_alias = config['tools']['plotting']['use_alias']
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/plotting/{sample}.plotting.benchmark'
    shell:
        "Rscript ../scripts/scRNA_pipeline_plotting.R  --sce_in {input.infile} --genelist {params.genes_of_interest} --outDir {params.outputDirec} --sampleName {params.sampleName} --colour_config {params.colour_config} --toggle_label {params.use_alias}"


# adapt this directory in master snake file to prevent recomputing the cohort in each analysis
### THIS SHOULD PROBABLY BE CHANGED WITH NEW STRUCTURE. CLEARIFY! TODO 

# perform assembly of the nonmalignant reference cohort
rule assemble_nonmalignant_cohort:
    input:
        inputDir_hdf5 = config['tools']['assemble_non_malignant_reference']['hdf5_dir'],
        inputDir_classify = config['tools']['assemble_non_malignant_reference']['celltype_dir']
    output:
        outfile = 'results/non_malignant_reference/nonmalignant_reference_cohort.h5'
    params:
        outDir = 'results/non_malignant_reference/',
        non_malignant_types = config['tools']['assemble_non_malignant_reference']['non_malignant_types'],
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/non_malignant_reference/nonmalignant_reference_cohort.benchmark'
    shell:
        "{config[tools][assemble_non_malignant_reference][call]} --hdf5_dir {input.inputDir_hdf5} --celltype_dir {input.inputDir_classify} --non_malignant_celltypes {params.non_malignant_types} --out_dir {params.outDir} "

# plot reference cohort tSNEs
rule plot_tSNEs_nonmalignant_cohort:
    input:
        inputHDF5 = 'results/non_malignant_reference/nonmalignant_reference_cohort.h5'
    output:
        outfile_batch = 'results/non_malignant_reference/nonmalignant_reference_cohort.tSNE_batch.png',
	outfile_ct = 'results/non_malignant_reference/nonmalignant_reference_cohort.tSNE_celltype.png'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time']
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/non_malignant_reference/plot_tSNEs_nonmalignant.benchmark'
    shell:
        "{config[tools][plot_tSNE_nonmalignant][call]} --hdf5File {input.inputHDF5} --outFile_batch {output.outfile_batch} --outFile_ct {output.outfile_ct} "


# perform the differential expression analysis using a Wilcoxon test
rule diff_exp_genes:
    input:
        sce_in = 'results/atypical_removed/{sample}.RDS',
        cell_types = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
        #outpath = dynamic('results/diff_exp/' + "{sample}.{clusterid}.DEgenes.tsv"),
        success = 'results/diff_exp/{sample}.diffExp_success.txt'
    params:
        sampleName = '{sample}',
        malignant = config['inputOutput']['malignant_cell_type'],
        threshold_comparison = config['tools']['diff_exp']['threshold_comparison'],
        fdr_cut = config['tools']['diff_exp']['fdr_cut'],
        fc_cut = config['tools']['diff_exp']['fc_cut'],
        mindiff2second = config['tools']['diff_exp']['mindiff2second'],
        minNumberNonMalignant = config['tools']['diff_exp']['minNumberNonMalignant'],
        outpath = 'results/diff_exp/'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/diff_exp/{sample}.diff_exp.benchmark'
    shell:
        "Rscript ../scripts/apply_DE_analysis.R --sample_data {input.sce_in} " +
        "--sampleName {params.sampleName} " +
        "--cluster_table {input.cell_types} " +
        "--malignant_tag {params.malignant} " +
        "--fdr_cut {params.fdr_cut} " +
        "--fc_cut {params.fc_cut} " +
        "--mindiff2second {params.mindiff2second} " +
        "--threshold_comparison {params.threshold_comparison} " +
        "--minNumberNonMalignant {params.minNumberNonMalignant} " +
        "--outdir {params.outpath} ; " +
        "date > {output.success} "


# give out gene expression values per cluster
rule gene_exp:
    input:
        sce_in = 'results/atypical_removed/{sample}.RDS',
    output:
        out = 'results/gene_exp/{sample}.gene_expression_per_cluster.tsv'
    params:
        sampleName = '{sample}',
        outpath = 'results/gene_exp/',
        threshold_sample = config['tools']['gene_exp']['threshold_sample'],
        type_sample = config['tools']['gene_exp']['type_sample'],
        priority_genes = config['resources']['priority_genes'],
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gene_exp/{sample}.gene_exp.benchmark'
    shell:
        "Rscript ../scripts/get_cluster_gene_expression.R --sce_in {input.sce_in} --priority_genes {params.priority_genes} --filtering_threshold_sample {params.threshold_sample} --filter_type_sample {params.type_sample} --outDir {params.outpath} --sample_name {params.sampleName} "


# This rule generates general quality control plots to hdf5 expression files
rule generate_qc_plots :
    input:
        infile = '{sample}.h5'
    output:
        out = '{sample}.h5.histogram_library_sizes.png'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        '{sample}.generate_qc_plots.benchmark'
    shell:
        'Rscript ../scripts/generate_QC_plots.R --hdf5File {input.infile} '


# This rule creates a box plot comparing cell type fractions across samples
rule generate_cell_type_boxplot:
    input:
        previous_samples = config['resources']['previous_samples'],
        sample_cell_types = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
        out = 'results/plotting/{sample}.boxplot_cell_types_cohort.png'
    params:
        sampleName = '{sample}',
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
        outDir = 'results/plotting/',
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/plotting/{sample}.boxplot_cell_types_cohort.benchmark'
    shell:
        'Rscript ../scripts/generate_boxplot_fractions_celltypes.R --previous_samples {input.previous_samples} --current_sample {input.sample_cell_types} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --outDir {params.outDir}'



# This rule integrates samples of the cohort and visualizes the integration with UMAPs
rule sample_integration:
    input:
        previous_samples = config['resources']['previous_samples_counts'],
        current_sample = 'results/atypical_removed/{sample}.RDS'
    output:
        out = 'results/plotting/{sample}.sample_integration_highlight_current.png'
    params:
        sampleName = '{sample}',
        outDir = 'results/plotting/',
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
	colour_config = config['resources']['colour_config']
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/plotting/{sample}.sample_integration.benchmark'
    shell:
        'Rscript ../scripts/sample_integration.R  --cohort_list {input.previous_samples} --sample_data {input.current_sample} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --colour_config {params.colour_config} --outdir {params.outDir}'



# calculate for each cluster the number of cells it countains and the percentage of all cells
rule cellPercentInCluster:
    input:
        clusterCsv = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
#    print(drugID)
        out = 'results/clusterpercent/{sample}.clusters_cell_count_percent.txt'
    params:
        variousParams = config['tools']['cellPercentInCluster']['variousParams']
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/clusterpercent/{sample}.clusterPercent.benchmark'
    shell:
        'python ../scripts/count_cells_in_clusters.py  --inputTable {input.clusterCsv} --outFile {output.out} {params.variousParams}'
