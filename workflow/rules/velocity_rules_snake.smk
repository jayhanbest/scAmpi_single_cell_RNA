# perform preliminary unspliced/splice mRNA evaluation necessary for the velocity analysis
rule velocyto:
    input:
        cellranger_possorted_bam = 'results/cellranger_run/{widcards.sample}/outs/possorted_genome_bam.bam',
        barcode_file = "results/cellranger_run/{widcards.sample}.barcodes.tsv"
    output:
        outdir = 'results/velocity/velocyto/{widcards.sample}/',
        outfile = 'results/velocity/velocyto/{widcards.sample}/{widcards.sample}.loom'
    conda:
        '../envs/velocyto_v0.17.17.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
        mask_file = config['resources']['mask_file'],
        genome_annotation = config['resources']['genome_annotation']
    threads:
        config['computingResources']['highRequirements']['threads']
    benchmark:
        'results/velocity/velocyto/{widcards.sample}.velocity_velocyto.benchmark'
    shell:
        'velocyto run --bcfile {input.barcode_file} --outputfolder {output.outdir} --sampleid {widcards.sample} --mask {resources.mask_file} -@ {threads} {input.cellranger_possorted_bam} {resources.genome_annotation}'

# perform a complete velocity analysis on all available samples
rule velocity_anaysis:
    input:
        loom_file = 'results/velocity/velocyto/{widcards.sample}/{widcards.sample}.loom',
        sce_file = 'results/atypical_removed/{widcards.sample}.atypical_removed.RDS'
    output:
        outdir = 'results/velocity/{widcards.sample}/'
    conda:
        '../envs/scvelo_0.2.4.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
        hex_map_file = config['resources']['hex_map_file']
    params:
        n_neighbors = config['tools']['scvelo']['n_neighbors'],
        n_top_gene_plots = config['tools']['scvelo']['n_top_gene_plots'],
        n_top_genes = config['tools']['scvelo']['n_top_genes'],
        min_shared_counts = config['tools']['scvelo']['min_shared_counts'],
        embedding = config['tools']['scvelo']['embedding'],
        sce_variable = config['tools']['scvelo']['sce_variable']
    threads:
        config['computingResources']['highRequirements']['threads']
    benchmark:
        'results/velocity/{widcards.sample}/{widcards.sample}.velocity.benchmark'
    shell:
        'python workflow/scripts/run_scvelo.py --loom_file {input.loom_file} --rds_file {input.sce_file} --sample_name {widcards.sample} --hex_map_path {resources.hex_map_file} --embedding {params.embedding} --myvariable {params.sce_variable} --out_dir {output.outdir} --n-jobs {threads} --n_neighbors {params.n_neighbors} --n_top_gene_plots {params.n_top_gene_plots} --n_top_genes {params.n_top_genes} --min_shared_counts {params.min_shared_counts}'

