#! /usr/bin/env python
# -*- coding: utf-8 -*-
#a='a'
"""
Implementation of scvelo for RNA velocity estimation & analysis of scRNA-seq data w/ spliced & unspliced counts data.
https://scvelo.readthedocs.io/

Original structure: Marcus Lindberg, March 2020
Matteo Carrara November 2021

Requirements:
module load new r/3.6.0
conda activate rna_velocity

Important remarks:
    Dependencies:
        igraph is necessary for PAGA and for velocity clustering. The correct conda package is python-igraph
	louvain is necessary for velocity clustering. The correct conda package is louvain. WARNING: there is a package called python-louvain, that does not work. This is due to the fact that the package python-louvain is imported within the velocity clustering script as -import louvain- but the package should be imported as -community-. This cannot be changed as the import happens within the velocity clustering script from scVelo. The conda package -louvain- instead, allows to call directly -import louvain- 
    Special variables:
        The analysis includes one special variable called 'mysamplename'. Adding it at script call (--myvariable mysamplename) will trigger UMAP plots divided by sample name. This is useful when working with cohorts, in order to highlight the position of samples in the UMAPs and compare it with all other results
"""

import argparse
import sys
import scvelo as scv
from os import listdir, chdir, mkdir, getcwd
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
import pickle
from pathlib import Path
import matplotlib
import numpy as np
from difflib import SequenceMatcher
import copy
import anndata

parser = argparse.ArgumentParser(description = 'Estimates RNA velocity of spliced/unspliced counts data to generate velocity plots.')
parser.add_argument('--loom_file', dest = 'ldata_path', required = True, help = 'Path to .loom file containing spliced & unspliced counts matrices. If working with a cohort, it should be a list of paths of all loom files involved, one per line.')
parser.add_argument('--rds_file', dest = 'adata_path', required = True, help = 'Path to .RDS file containing processed SCE object w/ cell type annotations & UMAP coordinates. If working with a cohort, it should be a list of paths of all SCE object files involved, one per line, with the samples in the same order provided in the loom file list.')
parser.add_argument('--out_dir', dest = 'out_dir', required = True, help = 'Path to output directory for all plots.')
parser.add_argument('--sample_name', dest = 'sample_name', required = True, help = 'Name of sample to be used as prefix for all output files.')
parser.add_argument('--cohort_sample_names', dest = 'cohort_sample_names', required = False, help = 'Path to a list of sample names to assign to the provided cohort samplesi, one per line, with the samples in the same order provided in the loom file list.')
parser.add_argument('--colour_config_final_file', dest = 'colour_config_final_file', required = False, default = None, help = 'Path to file specifying colours to be used for specific final cell types in plots. (Default = none)')
parser.add_argument('--colour_config_major_file', dest = 'colour_config_major_file', required = False, default = None, help = 'Path to file specifying colours to be used for specific major cell types in plots. (Default = none)')
parser.add_argument('--hex_map_path', dest = 'hex_map_path', required = False, help = 'Path to file specifying mapping of colour names to hexadecimal values.')
parser.add_argument('--n_neighbours', dest = 'n_n', required = False, type = int, default = 20, help = 'Number of neighbors to be used for calculating 1st & 2nd order moments, as well as for plotting. (Default: 20)')
parser.add_argument('--n_top_gene_plots', dest = 'n_top_gene_plots', required = False, type = int, default = 10, help = 'Number of top velocity genes to generate individual velocity plots for. (Default: 10)')
parser.add_argument('--myvariable', action='append', dest = 'myvariable', required = True, help = ' Variable in the SingleCellExperiment object to be used for UMAP computation. This option can be provided multiple time and the script will loop through all of them. (Default: celltype_major) WARNING: The first element will be used for differential kinetic computation and correction')
parser.add_argument('--embedding', dest = 'myembedding', required = False, default = "umap", help = 'Embedding to use to plot the umaps. It must be an emebedding that exists in the input SingleCellEperiment object. (Default: umap)')
parser.add_argument('--min_shared_counts', dest = 'min_shared_counts', required = False, default = 30, help = 'Minimum number of counts (both unspliced and spliced) required for a gene during normalization')
parser.add_argument('--n_top_genes', dest = 'n_top_genes', required = False, default = 2000, help = 'Number of genes to consider during normalization. Not to be confused with option "n_top_gene_plot"')
parser.add_argument('--n_jobs', dest = 'n_jobs', required = False, default = 1, type = int, help = 'Number of concurrent jobs for parallelized steps')

args = parser.parse_args()
#scv.logging.print_version()

def print_break():
    print("#"*60)

###############################################################################
################################ I N P U T ####################################
###############################################################################

################################ LOAD DATA ####################################
'''
Load necessary data files:
  1. loom file w/ spliced/unspliced counts
  2. RDS file w/ cell type annotations & UMAP coordinates.
Files will be merged into one data object (AnnData h5ad format).
'''
def load_data(ldata_path, adata_path, sample_name):
    print_break()
    print("Reading RDS with processed counts & metadata: ", adata_path)
    anndata2ri.activate()
    r(f'temp_data <- readRDS("{adata_path}")')
    data_format = r('class(temp_data)')[0] #SingleCellExperiment only (not Seurat)
    print("Data class: ", data_format)
    if data_format != "SingleCellExperiment":
        sys.exit("Please provide an RDS containing a SingleCellExperiment object. Exiting...")
    adata = r(f'as(temp_data, "SingleCellExperiment")')
    anndata2ri.deactivate()
    adata.obs['mysamplename'] = sample_name
    print("Reading loom file: ", ldata_path)
    ldata = scv.read(ldata_path)
    if len((ldata.var["Accession"])) > len(set(ldata.var["Accession"])):
        "WARNING: Loom file contains duplicated genes."
    #The check below is dependent on loom file having barcodes in the following format:
    #  {sample_name}:ATCG...GTCAX (starting w/ sample name, ending in 'X')
    #Sanity checking that the barcodes have the expected format
    #This is true only if the loom file is not merged. When merging, the batch and the final X are already stripped
    #FIXME: find a way to check if the loom is merged or not, so that you know if you need to run this part or not
    #FIXME: a try-catch makes sense, but I have to find the correct element. Probably trying to access ldata.obs.sample_batch,
    #FIXME: which gives an AttributeError is the loom file is single sample
    #FIX IN PLACE: needs just some testing
#    try:
#        len(ldata.obs.sample_batch)
#        print("The loom file appears to contain multiple concatenated samples. Proceeding accordingly...")
#        ldata_barcodes = set(ldata.obs_names)
#    except AttributeError:
#    print("The loom file appears to contain a single sample. Proceeding accordingly...")
    barcodes_sc_last = [x[-1] for x in list(ldata.obs_names)]
    if any(not( (barcode.endswith('x') or barcode.endswith('X')) ) for barcode in barcodes_sc_last):
        sys.exit("Loom file contains one or more barcodes not ending with 'x' or 'X'. Exiting...")
    barcodes_sc_sn = [barcode[0:barcode.rfind(':')] for barcode in list(ldata.obs_names)]
    if any(len(barcode) == 0 for barcode in barcodes_sc_sn):
        sys.exit("Loom file contains one or more barcodes with no sample name. Exiting...")
    ldata_barcodes = set([barcode[barcode.rfind(':')+1:-1] for barcode in list(ldata.obs_names)])
#        barcodes_sc_last = [x[-1] for x in list(ldata.obs_names)]
#    if any(not( (barcode.endswith('x') or barcode.endswith('X')) ) for barcode in barcodes_sc_last):
#        sys.exit("Loom file contains one or more barcodes not ending with 'x' or 'X'. Exiting...")
#    barcodes_sc_sn = [barcode[0:barcode.rfind(':')] for barcode in list(ldata.obs_names)]
#    if any(len(barcode) == 0 for barcode in barcodes_sc_sn):
#        sys.exit("Loom file contains one or more barcodes with no sample name. Exiting...")
#
#    ldata_barcodes = set([barcode[barcode.rfind(':')+1:-1] for barcode in list(ldata.obs_names)])
    print("Number of cells in loom file (spliced & unspliced raw counts): " + str(len(ldata_barcodes)))
    #if the adata comes from an integrated file, the barcodes come with additional strings to be stripped
#    nucleotides = ['A', 'C', 'G', 'T']
#    adata_barcodes = ["".join([i for i in element if i in nucleotides]) for element in adata.obs_names]
#    adata_barcodes = set(list(adata_barcodes))
    adata_barcodes = set(list(adata.obs_names))
    print("Number of cells in processed SingleCellExperiment RDS: " + str(len(adata_barcodes)))
    matching_barcode_counts = len(ldata_barcodes.intersection(adata_barcodes))
    print(str(matching_barcode_counts) + " of " + str(len(ldata_barcodes)) + " barcodes in loom file found in processed RDS.")
    print("Corresponding to " + str("{:.2f}".format(matching_barcode_counts/len(ldata_barcodes)*100)) + "% of loom barcodes and " + str("{:.2f}".format(matching_barcode_counts/len(adata_barcodes)*100)) + "% of SingleCellExperiment RDS barcodes")
    # TODO: There are seurat wrappers that helps working with integrated data. Skip the chunk completely and use a simple check for now
    if matching_barcode_counts == 0:
        print("No matching barcodes found. Double check input files. If data is integrated, make sure to pre-process combined loom file so that it adheres to the barcode format in the SCE RDS.")
    if matching_barcode_counts < len(adata_barcodes)/10:
        sys.exit("WARNING: the provided data files seem to have few (< 10%) matching barcodes. Double check the data files provided.Exiting...")
#
#    sufficient_barcode_match = False
#    while(sufficient_barcode_match == False):
#        if matching_barcode_counts == 0:
#            print(f"WARNING: the provided data files seem to have no matching barcodes.\nDouble check the data files provided and make sure they are of the same sample.\n\nExample loom barcode: {list(ldata.obs_names)[0]} |  Example RDS barcode: {list(adata.obs_names)[0]}\n")
#            print("Data might be integrated. Attempting direct match of barcodes...")
#            #Re=attempt matching by doing direct match.
#            ldata_barcodes = set(list(adata.obs_names))
#            matching_barcode_counts = len(ldata_barcodes.intersection(adata_barcodes))
#            print(str(matching_barcode_counts) + " of " + str(len(ldata_barcodes)) + " barcodes in loom file found in processed RDS.")
#            if matching_barcode_counts == 0:
#                print("No matching barcodes found. Double check input files. If data is integrated, make sure to pre-process combined loom file so that it adheres to the barcode format in the SCE RDS.")
#                sys.exit("Exiting...")
#        elif matching_barcode_counts < len(adata_barcodes)/10:
#            sys.exit("WARNING: the provided data files seem to have few (< 10%) matching barcodes. Double check the data files provided. Exiting...")
#        elif matching_barcode_counts < len(adata_barcodes)/2:
#            print("WARNING: less than half of barcodes in the processed RDS were found in the spliced counts data. Double check your input files.\n")
#            sufficient_barcode_match = True
#            print("Proceeding...")
#        else:
#            sufficient_barcode_match = True
#
    for barcode in ldata_barcodes.symmetric_difference(adata_barcodes):
        if barcode not in ldata_barcodes:
            print("WARNING: there are barcodes in the processed data that are not in the unprocessed data. Double check the input files.")
            break
    print_break() 
    print("Merging loom file and RDS into single file.")
    data_merged = scv.utils.merge(adata, ldata)
    del adata
    del ldata
    print_break()
    return data_merged


'''
Load preprocessed, merged h5ad data file from previous run. Missing data not stored in object loaded from pickled files.

def load_merged(data_path):
    print_break()
    data = scv.read(data_path)
    try:
        data.uns = pickle.load(open(f"{data_path[:-12]}.uns.pickle", "rb"))
        data.obsm = pickle.load(open(f"{data_path[:-12]}.obsm.pickle", "rb"))
    except FileNotFoundError:
        scv.logging.warn("Processed data detected, but slots missing. Re-running analysis...")
        data = estimate_velocity(data, out_path, n_neighbors)
        pickle.dump(data.uns, open(f"{data_path[:-12]}.uns.pickle", "wb"))
        pickle.dump(data.obsm, open(f"{data_path[:-12]}.obsm.pickle", "wb"))
        print("Saved fully processed velocity data.")
    else:
        print("Successfully loaded fully processed velocity data.")
    print_break()
    return data
'''
########################## CALCULATE RNA VELOCITY ##############################
'''
Preprocess & process data for RNA velocity estimation (using dynamical model) w/ awareness for cell type-specific kinetics.
'''
def estimate_velocity(data, out_path, n_neighbors, min_shared_counts, n_top_genes, sample_name, myembedding, myjobs):
    print_break()
    print("Number of cells in whole sample: ", len(data.obs['barcodes']))
    
    #If the number of cells are limited (less than 20), bypass the user-provided number of neighbors and use n_cells-1. Send a warning if that happens
    n_n = get_n_neighbors(data, n_neighbors)
    
    #Print and plot the proportion of spliced and unspliced RNA
    scv.utils.show_proportions(data)
    plt = scv.pl.proportions(data)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.proportions.png")

    try:
        scv.pp.filter_and_normalize(data, min_shared_counts, n_top_genes, log = True)
    #TODO: Understand why we exclude filter_genes_dispersion() if we fail to run the wrapper above
    except KeyError:
        print("Error in filtering step. Re-attempting without normalizing gene dispersion...")
        scv.pp.filter_genes(data)
        scv.pp.normalize_per_cell(data)
        #scv.pp.filter_genes_dispersion(data)
        scv.pp.log1p(data)
    #Compute first and second order moments (means and uncentered variances) for deterministic and stochastic velocity estimation respectively
    scv.pp.moments(data, n_pcs=30, n_neighbors=n_n)
    
    #Depending on the method used and on the sample being cohort or not, the umap may be called differently
    #We search for a perfect match between all keys and myembedding. If no match is available, find a match that contains the full string of myembedding
    #Report any other partial or full matches in any case. Fail if there are no perfect matches and more than 1 partial match
    #We use iterable unpacking with the * operator
    res = [i for i in [*data.obsm] if myembedding in i]
    res2 = [i for i in [*data.obsm] if myembedding == i]
    if len(res) == 0:
        sys.exit("FATAL: no embedding matching (completely nor partially) the provided embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) == 0:
        sys.exit("More than 1 embedding found partially matching the selected embedding and no perfect match. Cannot decide which one to select. Please provide a more precise embedding.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) > 1:
        sys.exit("Multiple perfect matches of the selected embedding. This should never happen, please check original objects.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
    if (len(res) > 1) and len(res2) == 1:
        print("One or more partial matches of the selected embedding detected alongside one perfect match. Continuing with the perfect match.\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res2[0]
    if len(res2) == 1:
        print("One perfect match found of the selected embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res2[0]
    if len(res) == 1 and len(res2) ==0:
        print("One partial match found of the selected embedding\nProvided: " + myembedding + " - Available keys: " + str([*data.obsm]))
        myembedding = res[0]
    #create a dedicated X_umap for the analysis where the embedding is stored
    data.obsm["X_umap"] = data.obsm[myembedding]
    myembedding="umap"
    #if "umap_hvg" in data.obsm:
    #    data.obsm["X_umap"] = data.obsm["umap_hvg"]
    #elif "X_umap" in data.obsm:
    #    print("Loading UMAP coordinates.")
    #else:
    #    sys.exit("UMAP coordinates missing. Please make sure RDS is fully processed.")
    scv.tl.recover_dynamics(data, n_jobs=myjobs)  # this step required for dynamical mode; time-intensive
    #Start by computing the velocity without differential kinetics. This can be done afterwards
    scv.tl.velocity(data, mode='dynamical', diff_kinetics=False)
    scv.tl.velocity_graph(data, n_jobs=myjobs)
    scv.tl.recover_latent_time(data)
    print_break()
    return data


##################### SUBSET DATA BY VARIABLE #######################
def subset_data(data, myvar):
    print("\n")
    print(f"Subsetting main object by {myvar}")
    try:
        cell_type_subset = data.obs.index[(data.obs[myvar]==cell_type)]
    except KeyError:
        sys.exit(f"The selected variable {myvar} does not exist in the object. Please check the available columns in your input SingleCellExperiment object and retry. Exiting.")
    data_subset = data[cell_type_subset, :].copy()
    print("Number of cells in cell type: ", len(data_subset.obs['barcodes']), "\n")
    return data_subset

'''
Subset data based on major cell type label & re-compute RNA analysis.
'''

#def estimate_velocity_by_celltype(data, cell_type):
#    print("\n")
#    cell_type_subset = data.obs.index[(data.obs['celltype_major']==cell_type)]
#    data_subset = data[cell_type_subset, :].copy()
#    print("Number of cells in cell type: ", len(data_subset.obs['barcodes']))
#    n_n = get_n_neighbors(data_subset, args.n_n)
#    scv.utils.show_proportions(data_subset)
#    #scv.pp.neighbors(data_subset, n_pcs=30, n_neighbors=n_n)
#    scv.pp.moments(data_subset, n_pcs=30, n_neighbors=n_n)
#    try:
#        scv.tl.recover_dynamics(data_subset)
#    except KeyError:
#        scv.logging.warn("Failed to recover dynamics for subtype; skipping.")
#        pass
#    scv.tl.velocity(data_subset, mode='dynamical', diff_kinetics=False)
#    scv.tl.velocity_graph(data_subset)
#    scv.tl.recover_latent_time(data_subset)
#    print("\n")
#    return data_subset


###############################################################################
############################### O U T P U T ###################################
###############################################################################

################################ SAVE DATA ####################################
'''
Save processed data in AnnData (.h5ad - custom hdf5) formatted file.
Slots in object not compatible w/ .h5ad file (e.g. data.uns slot) will be saved (pickled) separately.
'''
def save_h5ad(out_path, sample_name, data):
    try:
        data.write_h5ad(f'{out_path}/{sample_name}_merged.h5ad')
    except NotImplementedError as e:
        scv.logging.warn(e)
        print("Failed to save processed data into one file. Saving data separately...")
        print(f"Processed .h5ad file saved as {out_path}/{sample_name}_merged.h5ad.")
        pickle.dump(data.uns, open(f"{out_path}/{sample_name}.uns.pickle", "wb"))
        print(f"Saved unsorted velocity annotations as {out_path}/{sample_name}.uns.pickle.")
        pickle.dump(data.obsm, open(f"{out_path}/{sample_name}.obsm.pickle", "wb"))
        print(f"Saved velocity UMAP coordinates as {out_path}/{sample_name}.obsm.pickle.")
        pass
    else:
        print(f"Full data saved successfully. Processed .h5ad file saved as {out_path}/{sample_name}_merged.h5ad.")


################################## PLOTS ######################################
'''
Generate & save velocity plots: UMAP embedding stream, heatmap of top genes,
latent-time scatter plots, & individual gene plots.
Uses color config & labels from config passed by plot_scvelo fxn.
Only to be called by plot_scvelo fxn.
'''
def save_plots(sample_name, data, config, n_top_gene_plots, n_neighbors, myembedding, out_path, myvariable):
    top_genes = data.var['fit_likelihood'].sort_values(ascending=False).index
    if(n_top_gene_plots > len(top_genes)):
        print(f'WARNING: Asked to plot more genes ({n_top_genes_plots}) than the maximum number of top variable genes detected ' + len(top_genes) + ". Dropped the number of genes to plot to the total amount of top variable genes available. Please consider reducing the value of the variable n_top_genes_plots")
        max_genes = len(top_genes)
    else:
        max_genes = n_top_gene_plots

    matplotlib.use('Agg')
    n_n = get_n_neighbors(data, n_neighbors)
    if n_n < 20:
        print("Sample size too small. Skipping plots...")
        pass
    #TODO: This part is for when you are NOT working single sample. Skipping until fully implemented
    #if(scv.pl.utils.settings.plot_prefix[:-1] == sample_name):
    #    sample_name_extended = sample_name + ".all"
    #    sample_name = "all"
    #else:
    #    sample_name_extended = scv.pl.utils.settings.plot_prefix[:-1] + "." + sample_name

    if len(top_genes) > len(set(top_genes)):
        print("WARNING: possible duplicated genes:")
        print(set([gene for gene in top_genes if top_genes.count(gene) > 1]))
    #try:

        #scv.pl.velocity_embedding_stream(data, basis='umap', color='phenograph_clusters',
        #    save=f"{sample_name}.phenograph_clusters.png", fontsize=20,
        #    figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
        #    title=config['title']+" velocity - phenograph clusters", palette=config['phenograph_palette'], n_neighbors=n_n)

        #FIXME: Right now this part is added in the try, but no exception is handled. The exception will be completely reworked after these basic features are in
        #Plot the velocity vector field at single-cell level. This works as a more fine-grained view of the previous two plots.
        #If there are too many vectors, you can play with the option "density"

        #scv.pl.velocity_embedding(data, arrow_length=3, arrow_size=2, dpi=200, color='phenograph_clusters',
        #        save=f"{sample_name}.sc_velocity_phenograph_clusters.png", fontsize=20,
        #        figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
        #        title=config['title']+" single-cell velocity vectors", palette=config['phenograph_palette'])

    plt = scv.pl.scatter(data, basis=myembedding, color='latent_time', #size=100,
        color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0, 1],
        show=False, legend_loc=config['legend_loc'],
        title=config['title'] + " latent time", figsize=(14,10))
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.latent_time.png")
    gene_path = out_path + "top_genes/"
    Path(gene_path).mkdir(parents=True, exist_ok=True)
    for myvar in myvariable:
        print(f'Plotting for variable {myvar}')
        plt = scv.pl.velocity_embedding_stream(data, basis=myembedding, color=config[myvar+'_colours'],fontsize=20,
            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
            title=config['title']+" velocity", palette=config[myvar+'_palette'], n_neighbors=n_n)
        matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity.{myvar}.png")
        scv.pl.velocity_embedding(data, arrow_length=3, arrow_size=2, dpi=200, color=config[myvar+'_colours'], fontsize=20,
            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
            title=config['title']+" single-cell velocity vectors", palette=config[myvar+'_palette'])
        matplotlib_save(plt, fname=f"{out_path}/QC_{sample_name}.sc_velocity_{myvar}.png")
        plot_speed_coherence(data, myvar, out_path, sample_name, myembedding)

        #FIXME: fails internally, cannot find 'transition_confidence'
        #plot_paga_graph(data, myvar, out_path, sample_name, myembedding)
        #TODO: The kinetic rate parameter function is set but not called
        #kinetic_rate_params()
            
        # The heatmap forces a prefix called "heatmap" that goes BEFORE the value set for scv.pl.utils.settings.figdir (but AFTER the content of scv.pl.utils.settings.plot_prefix). THIS IS A BUG OF SCVELO.
	    # The workaround is to save the plot using matplotlib. For consistency, all plots are saved with matplotlib
        plt = scv.pl.heatmap(data, var_names=top_genes[0:300], tkey='latent_time', n_convolve=100, col_color=myvar, show=False, font_scale=0.5, xticklabels=False)
        matplotlib_save(plt, fname=f"{out_path}/{sample_name}.top_driver_genes_heatmap.{myvar}.png")

        plot_velocity_clustering(data, myvar, out_path, sample_name, myembedding, config)

        print('Generating plots for top variable genes.')

        for gene in top_genes[0:max_genes]:
            print(f'Working on gene {gene}')
            plt = scv.pl.scatter(data, basis=myembedding, var_names=gene,color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, linewidth=2, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{sample_name}.{gene}.phase_portrait.{myvar}.png")

            plt = scv.pl.scatter(data, basis=myembedding, x='latent_time', y=gene, color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, linewidth=2, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{sample_name}.{gene}.phase_latent_time.{myvar}.png")

            plt = scv.pl.velocity(data, basis=myembedding, var_names=gene,color=config[myvar+'_colours'], color_map = "RdBu_r", colorbar=True, legend_loc='none', size = 15, linewidth=2, figsize=(20,15), show=False)
            matplotlib_save(plt, fname=f"{gene_path}/{sample_name}.{gene}.multiplot.{myvar}.png")

        #FIXME: check the whole try-catch statement. You can integrate it as an exception in the normal workflow, instead of repeating the workflow if the exception arises
    #except ValueError:
    #    print(f"{sample_name} data set too small ({data.n_obs} cells). Velocity plots not saved.")
    #    print("Attempting workaround...")
    #    # For edge cases, usually w/ small number of cells. Need to scale values down.
    #    data.obsm['X_umap'] /= np.max(np.abs(data.obsm['X_umap']))
    #    try:
    #        scv.pl.velocity_embedding_stream(data, basis='X_umap', color=config['colour'],
    #            save=f"{sample_name}_velocity.png", fontsize=20, legend_fontsize = 10,
    #            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
    #            title=config['title']+" velocity", palette=config['palette'], n_neighbors=n_n)
    #    except ValueError:
    #        print(f"Work around failed. Skipping {sample_name} plots.")
    #        pass
    #else:
    #    scv.pl.scatter(data, color='latent_time', #size=100,
    #        color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0, 1],
    #        save=f"{sample_name}.latent_time.png", show=False, legend_loc=config['legend_loc'],
    #        title=config['title'] + " latent time", palette=config['palette'], figsize=(14,10))
    #
    #    try:
    #        scv.pl.heatmap(data, var_names=top_genes, tkey='latent_time', n_convolve=100, col_color="phenograph_clusters",save=f"{sample_name}.topdrivergenes-heatmap.phenograph.png", show=False, xticklabels=False, yticklabels=False)
    #       #Save heatmap values: heatmap.data.to_csv(); no barcode information
    #        if sample_name == "all":
    #            scv.pl.heatmap(data, var_names=top_genes, tkey='latent_time', n_convolve=100, col_color=config['col_colour'],save=f"{sample_name}.topdrivergenes-heatmap.png", show=False, xticklabels=False, yticklabels=False)
    #            
    #    except ValueError:
    #        print("Too few cells. Skipping heatmap.")
    #        pass
    #
    #    try:
    #        scv.pl.velocity_embedding_stream(data, basis='umap', color="seurat_clusters",
    #            save=f"{sample_name}.seurat_clusters.png", fontsize=20,
    #            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
    #            title=config['title']+" velocity - Seurat clusters", n_neighbors=n_n)
    #    except KeyError:
    #        pass
    #    except ValueError:
    #        print("No Seurat cluster colour config detected; skipping...")
    #        pass

'''
Configure plot parameters (colors, margin, labels, size, title).
Passes config to save_plots.
'''
def plot_configuration(sample_name, data, out_path, myvariable, hex_map_path = None, colour_config_final_path = None, colour_config_major_path = None):
    print("Setting up the plot configuration")
    plot_config = {  # "legend_loc":"upper right",
        "legend_loc": "right",
        #"colour": "celltype_final",
        #"col_colour": "celltype_final",
        "min_mass": "0",  # 0 - all trajectories; 100 - large only
        "title": sample_name.replace('_', '-'),
        "layer": ""
        }
    phenograph_clusters_colours = ["red", "green", "blue", "cyan", "yellow", "purple", "brown","chocolate", "chartreuse", "darkgoldenrod", "steelblue", "slateblue", "olivedrab", "gold", "violetred", "darkcyan","orchid", "darksalmon", "darkslategrey", "khaki", "indianred","magenta", "slategray", "darkolivegreen", "mediumaquamarine", "hotpink", "darkorange", "bisque", "darkseagreen", "dodgerblue", "deeppink", "tan", "mediumorchid"]
    for myvar in myvariable:
        plot_config[myvar+'_colours'] = myvar
        plot_config[myvar+'_palette'] = phenograph_clusters_colours
    if(colour_config_final_path): 
        cell_types = list(data.obs["celltype_final"].cat.categories)
        colour_config = read_colour_config(colour_config_path = colour_config_final_path, hex_map_path = hex_map_path)
        data.uns['celltype_final_colours'] = []    
        #match each celltype with the more similar entry in the color config available and map the color to it
        #data.uns['celltype_final_colours'] = map_type_colours(cell_types, colour_config)
        #min_path = out_path + sample_name
        colour_config_map = {}
        #Path(min_path).mkdir(parents=True, exist_ok=True)
        #FIXME: chdir are NOT suggested as you easily lose track of where you are in complex pipelines. Remove this part and use absolute paths
        #chdir(min_path)
        #plot_config["celltype_final_palette"] = data.uns['celltype_final_colours']
        plot_config["celltype_final_palette"] = map_type_colours(cell_types, colour_config)
    else:
        #Create and use a standard palette
        data.uns['celltype_final'] = phenograph_clusters_colours
    if(colour_config_major_path): 
        cell_types = list(data.obs["celltype_major"].cat.categories)
        colour_config = read_colour_config(colour_config_path = colour_config_major_path, hex_map_path = hex_map_path)
        data.uns['celltype_major_colours'] = []    
        #match each celltype with the more similar entry in the color config available and map the color to it
        #data.uns['celltype_major_colours'] = map_type_colours(cell_types, colour_config)
        #min_path = out_path + sample_name
        colour_config_map = {}
        #Path(min_path).mkdir(parents=True, exist_ok=True)
        #FIXME: chdir are NOT suggested as you easily lose track of where you are in complex pipelines. Remove this part and use absolute paths
        #chdir(min_path)
        #plot_config["celltype_major_palette"] = data.uns['celltype_major_colours']
        plot_config["celltype_major_palette"] = map_type_colours(cell_types, colour_config)
    else:
        plot_config["celltype_major_palette"] = phenograph_clusters_colours
    plot_config["phenograph_clusters_palette"] = read_colour_config(hex_map_path = hex_map_path, manual_colour_list = phenograph_clusters_colours)
    #scv.settings.set_figure_params(frameon=True, fontsize=20,
    plot_config['generic_palette'] = read_colour_config(hex_map_path = hex_map_path, manual_colour_list = phenograph_clusters_colours)
    print(plot_config)
    #save_plots(sample_name, data, plot_config, args.n_top_gene_plots)
    return(plot_config)

########################## Plotting helper fxns ###############################
'''
Sets additional plot settings. (Removes default saved plot file prefix & sets a custom one.)

def set_plot_settings(out_path, prefix = None):
    scv.pl.utils.settings.figdir = out_path
    scv.pl.utils.settings.file_format_figs="png"
    if(prefix):
        scv.pl.utils.settings.plot_prefix = f'{prefix}.'
    else:
        scv.pl.utils.settings.plot_prefix = ''
    scv.pl.utils.settings.plot_suffix = ''
'''

'''
Helper fxn to set appropriate number of neighbours if sample size is smallcelltype_final.
Default is 20. (Can be passed as command-line arg.)
'''
def get_n_neighbors(data, n_neighbors):
    n_cells = data.n_obs
    if n_cells < 20:
        print(f"Warning: The sample size is below 20. The value of variable n_n passed at command line is bypassed and set to {n_cells -1}")
        n_n = n_cells - 1
    else:
        n_n = n_neighbors
    return n_n

'''
Helper function for calculating similarities between cell type labels.
Matches are determined once maximal similarity scores are found.
'''
def match_cell_type(type_to_match, reference_list):
    from difflib import SequenceMatcher
    def similar(a, b):
        return SequenceMatcher(None, a, b).ratio()
    cell_type = ''
    previous_sim_score = 0
    for each in reference_list:
        sim_score = similar(each, type_to_match)
        if (sim_score < previous_sim_score):
            continue
        else:
            cell_type = each
            previous_sim_score = sim_score
    return cell_type

'''
Maps color names from color config file to color hexadecimal values (required by plotting functions).
'''
def read_colour_config(hex_map_path, colour_config_path = None, manual_colour_list = None):
    colour_map = {}
    print("Colour config used: ", colour_config_path)
    print("Hex colour mapping used: ", hex_map_path)
    if colour_config_path:
        with open(colour_config_path) as f:
            for line in f:
                t = line.split()
                colour_map[t[0]] = t[1]
    elif manual_colour_list:
        for i in range(len(manual_colour_list)):
            colour_map[i+1] = manual_colour_list[i]
            colour_map[str(i+1)] = manual_colour_list[i]
    else:
        #If no colour config is available, but also no manual color palette, create a default one for the existing values
        #TODO: complete this part
        sys.exit("ERROR: No colour config file provided for function read_colour_config call")
    if('cell_type' in colour_map.keys()):
        colour_map.pop('cell_type')
    hex_map = open(hex_map_path, "r").read()
    for line in hex_map.split('\n')[:-1]:
        name = line.split(',')[0]
        colour = line.split(',')[1]
        for k, v in colour_map.items():
            if v.lower() == name.lower() or v.lower() == name.replace(" ", "").lower():
                colour_map[k] = colour
            elif v.lower() == name.lower()[:-1] or v.lower() == name.replace(" ", "").lower()[:-1]:
                colour_map[k] = colour
    colour_map['unknown'] = "#000000"
    return colour_map

'''
Assigns colors to cell type labels based on hexadecimal value mapping.
Similar cell types (but w/ arbitrarily different annotation) will be assigned the same color.
'''
def map_type_colours(types_to_map, hex_colour_config):
    colour_map = {}
    print(types_to_map)
    for cell_type in types_to_map:
        matched = match_cell_type(cell_type, list(hex_colour_config.keys()))
        colour_map[cell_type] = hex_colour_config[matched]
    return colour_map


############################### MATPLOTLIB SAVE ##############################
def matplotlib_save(plt, fname):
    matplotlib.pyplot.savefig(fname, bbox_inches='tight')
    print(f"Saving figure to {fname}")
    matplotlib.pyplot.close()


############################### GENE LISTS ###################################
'''
Save top 300 most-likely velocity driver genes to text file.
Each gene is separated by newline.
'''
def save_top_genes_list(out_path, sample_name, data, save=False):
    #TODO: sample_name_extended is used in case of cohort analysis. Commented out at the moment
    #if sample_name == scv.pl.utils.settings.plot_prefix[:-1]:
    #  sample_name_extended = sample_name
    #else:
    #  sample_name_extended = scv.pl.utils.settings.plot_prefix[:-1] + "." + sample_name
    top_likelihood_genes = data.var['fit_likelihood'].sort_values(ascending=False)[:300].to_frame()
    top_likelihood_genes = top_likelihood_genes.rename_axis("gene").reset_index()
    top_genes = top_likelihood_genes['gene'].values
    print(f"{sample_name} - top velocity genes: ", list(top_genes)[:11])
    if save:
        top_likelihood_genes.to_csv(f'{out_path}/{sample_name}.top_likelihood-ranked_genes.tsv', sep = '\t')
        print(f"{sample_name} velocity gene list saved.")


############################### GENE LISTS BY VARIABLE ######################
def save_top_genes_list_byvar(out_path, sample_name, data, myvar, save=False):
    scv.tl.rank_dynamical_genes(data, groupby=myvar)
    df = scv.get_df(data, 'rank_dynamical_genes/names')
    if save:
        df.to_csv(f"{out_path}/{sample_name}.top_likelyhood-ranked_genes_by_{myvar}.tsv", index = False, sep = '\t')
        print(f"{sample_name} velocity gene list by {myvar} saved.")


############################### SPEED AND COHERENCE #########################

def plot_speed_coherence(data, myvar, out_path, sample_name, myembedding):
    print('Plotting speed of rate of differentiation and coherence of the vector field')
    scv.tl.velocity_confidence(data)
    keys = 'velocity_length', 'velocity_confidence'
    plt = scv.pl.scatter(data, basis=myembedding, c=keys,cmap='coolwarm', perc=[2,98], size = 30, linewidth=1, legend_loc="right", legend_fontsize = 10, figsize=(13,10), show=False)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity-confidence.png")
    print (f'Saving the speed and coherence table divided by {myvar}')
    df = data.obs.groupby(myvar)[keys].mean().T
    df.style.background_gradient(cmap='coolwarm',axis=1)
    print(f'Saving the velocity length confidence table divided by {myvar}')
    df.to_csv(f'{out_path}/{sample_name}_velocity_length_confidence_{myvar}.tsv', sep = '\t')
    return(data)


############################## PAGA VELOCITY GRAPH #########################
## Additional trajectory inference method by showing the connectivity between clusters and velocity-inferred directionality
## TODO: iGraph installation not recognized. The function is not currently working under the rna_velocity environment
def plot_paga_graph(data, myvar, out_path, sample_name, myembedding):
    data.uns['neighbors']['distances'] = data.obsp['distances']
    data.uns['neighbors']['connectivities'] = data.obsp['connectivities']

    scv.tl.paga(data, groups=myvar)
    df = scv.get_df(data, 'paga/transition_confidence', precision=2).T
    df.style.background_gradient(cmap='Blues').format('{:.2g}')
    df.to_csv(f'{out_path}/{sample_name}.paga_transition_confidence.{myvar}.tsv', sep = '\t')
    plt = scv.pl.paga(data, basis=myembedding, size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, color=myvar, palette=config[myvar+'_palette'], linewidth=1, legend_loc="right", legend_fontsize = 10, figsize=(13,10), show=False)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.paga_graph.{myvar}.png")

    return data


############################# KINETIC RATE PARAMETERS ######################
## TODO: If it's useful, add it to the main plot call. The results in the example are underwhelming at best, as there is no distribution in the histogram, only isolated bars
def kinetic_rate_params(data, out_path, sample_name):
    df = data.var
    df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

    kwargs = dict(xscale='log', fontsize=16)
    with scv.GridSpec(ncols=3) as pl:
        plt = pl.hist(df['fit_alpha'], xlabel='transcription_rate', **kwargs)
        plt = pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
        plt = pl.hist(df['fit_gamma'], xlabel='degradation_rate', xticks=[.1, .4, 1], **kwargs)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.kinetic_rate_params.png")

    df = scv.get_df(data, 'fit*', dropna=True).head()
    df.to_csv(f"{out_path}/{sample_name}.kinetic_rate_params_df.tsv", index = False, sep = '\t')

    return data


############################# DIFFERENTIAL KINETICS #################
def diff_kinetics(data, sample_name, myvar, n_top_gene_plots, out_path, config, myjobs):
    print(f'Starting computation of differential kinetics. Cells in UMAPS will be colored by {myvar}')
    if n_top_gene_plots > 30:
        print(f'WARNING: n_top_gene_plots variable set too high, the resulting plot may be too big or difficult to read. Please consider lowering the value to 30 or lower')
    kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
    #Extract the top 300 most varible genes and compute the differential kinetic test on those. 300 is chosen to reduce computational time and avoid to introduce noise by including genes with lower likelihoods
    top_genes = data.var['fit_likelihood'].sort_values(ascending=False).index[:300]
    scv.tl.differential_kinetic_test(data, var_names=top_genes, groupby=myvar)
    plt = scv.pl.scatter(data, basis=top_genes[:n_top_gene_plots], ncols=5, add_outline='fit_diff_kinetics', color=config[myvar+'_colours'], colorbar=True, palette=config[myvar+'_palette'], size = 30, legend_loc="right", legend_fontsize = 12, figsize=(13,10), show=False, fontsize=20, **kwargs)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.top_genes_diff_kinetics.{myvar}.png")

    print("Recomputing velocity based on the differential kinetics")
    scv.tl.velocity(data, diff_kinetics=True)
    scv.tl.velocity_graph(data, n_jobs=myjobs)
    plt = scv.pl.velocity_embedding(data, dpi=200, arrow_size=2, arrow_length=2, color=config[myvar+'_colours'], fontsize=20,
            figsize=(14, 10), show=False, legend_loc=config['legend_loc'],
            title=config['title']+" single-cell velocity vectors - kinetic corrected", palette=config[myvar+'_palette'])
    matplotlib_save(plt, fname=f"{out_path}/QC_{sample_name}.sc_velocity_kinetic_corrected.{myvar}.png")
    return data


############################### VELOCITY CLUSTERING #########################
def plot_velocity_clustering(data, myvar, out_path, sample_name, myembedding, config):
    print(f'Clustering velocities and plotting the clusters against the values of {myvar}')
    scv.tl.velocity_clusters(data, match_with=myvar)
    mycolours = copy.deepcopy(config['generic_palette'])
    j=1
    for i in set(data.obs['velocity_clusters']):
        mycolours[i] = mycolours[j]
        del mycolours[j]
        j=j+1
    plt = scv.pl.scatter(data, basis=myembedding, color='velocity_clusters', size=30, legend_loc="right", legend_fontsize=10, figsize=(13,10), show=False, palette=mycolours)
    matplotlib_save(plt, fname=f"{out_path}/{sample_name}.velocity_clusters.{myvar}.png")
    return(data)

############################### RUN ALL ######################################

'''
Run main scvelo preprocessing, processing, & plotting workflow.
Velocities are estimated for whole sample & by cell type.
'''
def run_scvelo(ldata_path, adata_path, out_path, sample_name, colour_config_final_path, colour_config_major_path, hex_map_path, n_top_gene_plots, myvariable, n_neighbors, min_shared_counts, n_top_genes, myembedding, myjobs, cohort_sample_names):
    #try:
    #    output = load_merged(f"{out_path}/{sample_name}_merged.h5ad")
    #except (ValueError, OSError) as error:
    #    print("No existing .h5ad file from prior run detected. Starting pre-processing...")
    #    print("Sample name: ", sample_name)
    out_path = out_path + "/"
    print("Output dir: ", out_path)
    print("Testing if the provided data are for single samples or a cohort")
    try:
        with open(ldata_path) as f:
            ldata = f.readlines()
        ldata = [ name.rstrip() for name in ldata ]
        with open(adata_path) as f:
            adata = f.readlines()
        adata = [ name.rstrip() for name in adata ]
        with open(cohort_sample_names) as f:
            cohort_sample_names = f.readlines()
        cohort_sample_names = [ name.rstrip() for name in cohort_sample_names ]
        print("Detected cohort data. Proceeding accordingly...")
        ii = 0
        while ii < len(ldata):
            l = ldata[ii]
            a = adata[ii]
            s = cohort_sample_names[ii]
            print("Loading data for sample " + s)
            print("Sample " + str(ii+1) + " of " + str(len(ldata)))
            data_merged_temp = load_data(l, a, s)
            #data_merged_temp.var['sample_name'] = s
            if ii == 0:
                data = data_merged_temp
            else:
                data = data.concatenate(data_merged_temp, join = "inner")
            ii = ii + 1
    except (UnicodeDecodeError) as error:
        print("Detected single sample. Proceeding accordingly...")
        data = load_data(ldata_path, adata_path, sample_name)
    
    print(data)
    output = estimate_velocity(data, out_path, n_neighbors, min_shared_counts, n_top_genes, sample_name, myembedding, myjobs)
    myembedding = "umap"
    #save_h5ad(out_path, sample_name, output)
    #else:
    #    print("Loaded existing processed .h5ad file.")
    #
    print("Number of cells in whole sample: ", len(output.obs['barcodes']))
   
    ##TODO: First complete the framework for single samples, then try the integrated
    #if "sample_name" in output.obs:
    #    print("Integrated sample detected.")
    #    sub_samples = list(output.obs['sample_name'].cat.categories)
    #    print("Subsamples: ", ", ".join(sub_samples))
    #else:
    #    print("Single sample detected: ", sample_name)
    
    print_break()
    print("Current sample: ", sample_name)

    #TODO: decide if it's necessary to separate the results in different folders. See mostly when working with cohorts
    #out_path = out_path + sample_name + "/"
    #Path(out_path).mkdir(parents=True, exist_ok=True)
    
    plot_path = out_path + "plots/"
    Path(plot_path).mkdir(parents=True, exist_ok=True)
    print(f"Saving plots (whole sample) to {out_path}/plots/...")
    #set_plot_settings(plot_path)
    plot_config = plot_configuration(sample_name, output, plot_path, myvariable, hex_map_path, colour_config_final_path, colour_config_major_path)
    save_plots(sample_name, output, plot_config, n_top_gene_plots, n_neighbors, myembedding, plot_path, myvariable)
    print("Plots saved.")

    print_break()
    top_genes_path = out_path + "top_genes/"
    Path(top_genes_path).mkdir(parents=True, exist_ok=True)
    save_top_genes_list(top_genes_path, sample_name, output, save=True)

    diff_kin_path = plot_path + "differential_kinectics/"
    Path(diff_kin_path).mkdir(parents=True, exist_ok=True)
    for myvar in myvariable:
        print(f"Saving differential kinetics for variable {myvar}...")
        save_top_genes_list_byvar(top_genes_path, sample_name, output, myvar, save=True)
        diff_kinetics(output, sample_name, myvar[0:], n_top_gene_plots, diff_kin_path, plot_config, myjobs)

    writeLines("Run completed successfully", out_path + "/" + sample_name + ".velocity_success.txt")
    print_break()
    print(f"Completed velocity analysis. Output is available in {out_path}")
    
#    #TODO: understand if this framework by celltype makes sense and extracts more information
#    #It could make sense to stratify sub-populations of cells behaving in a different way
#    cell_types = list(output.obs['celltype_major'].cat.categories)
#    print("Cell types: ",cell_types)
#   
#   # Calculate cell type-specific velocities
# 
#    for cell_type in cell_types:
#        print_break()
#        print("Estimating RNA velocity for " + sample_name + " - " + cell_type + " ...")
#        scvelo_celltype = estimate_velocity_by_celltype(output, cell_type)
#        if "gene_ids" in scvelo_celltype.var:
#            gene_format = "gene_ids"
#        else:
#            gene_format = "Accession"
#        
#        print("Saving plots...")
#        plot_scvelo(cell_type, scvelo_celltype, out_path, hex_map_path, colour_config)
#        print("Plots saved.")
#        save_top_genes_list(out_path, cell_type, scvelo_celltype, save=True)

##############################################################################################

run_scvelo(args.ldata_path, args.adata_path, args.out_dir, args.sample_name, args.colour_config_final_file, args.colour_config_major_file, args.hex_map_path, args.n_top_gene_plots, args.myvariable, args.n_n, args.min_shared_counts, args.n_top_genes, args.myembedding, args.n_jobs, args.cohort_sample_names)
