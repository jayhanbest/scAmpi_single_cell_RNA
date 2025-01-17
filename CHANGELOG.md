# Changelog

## [2.0.4] - 2022-12-01

### Fixed

- adapt R scripts using ggplot2 for plotting (`filter_genes_and_cells.R`) to changes introduced with ggplot v3.4.0. With new defaults in `scale_*_manual` only factor levels found in the data are shown in the legend.
- in rule `cellranger_count` use consistently full path given in `config[inputOutput][input_fastqs]`. Update respective README sections.
- fix relative path to conda env yaml file in rule inheritance of `create_hdf5_starsolo`.

## [2.0.3] - 2022-11-21

### Fixed

- fix bug in script `filter_genes_and_cells.R` that resulted in colour discrepancy between legend and plot in `{sample}.visualize_filtered_cells.png` in rare cases.

## [2.0.2] - 2022-08-31

### Changed

- change run time keyword in rules and config from "time_min" or "time" to "runtime" as is recommended.
- have template memory values per job and not per thread
- use "mem_mb" consistently in rule files and config files

## [2.0.1] - 2022-07-19

### Changed

- Update example Ovarian cancer cell types to include Mesothelial cells and pDCs

### Fixed

- Fix empty output plot from rule `plot_upsetr`

## [2.0.0] - 2022-07-05

The scAmpi pipeline framework was fully revised to follow the current Snakemake best practices.
No changes were made regarding the analysis steps.

### Changed

- the Snakemake pipeline framework was updated to follow current best practices
- using Snakemake checkpoints the scAmpi basic and clinical parts were joined together to be one workflow
- several samples can be run in parallel (previously that was only possible for scAmpi_basic)

## [1.0.0] - 2022-04-05

First full, publicly available version of the scAmpi single-cell RNA analysis pipeline.
