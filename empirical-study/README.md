# Empirical Study

> [!NOTE]
> 
> Some of the files and directories mentioned below are only available in the Dryad data repository and are not included in the GitHub repository.

This directory contains the data, analyses, results, and figures from the analysis of green plants by InPhyNet. The empirical study validates InPhyNet's performance on biological data and demonstrates its practical utility for large-scale phylogenetic network inference.

## Key Files

- **`analysis.ipynb`**: Main Jupyter notebook containing the complete empirical analysis workflow
- **`mnet.net`**: Final inferred species network
- **`mnet_coded.net`**: Species network with coded internal nodes
- **`model_selection.csv`**: Model selection results comparing different inference approaches
- **`snaq_networks.net`**: Networks inferred using SNaQ for comparison

## Directories

### `okp_data/`
Contains the primary data used for the empirical analyses.

### `all_okp_data/`
Contains all data from the original One Thousand Plant Transcriptome Initiative analyses.

### `subset_data/`
Subset decompositions and intermediate files for divide-and-conquer analysis.

### `ToB_inference/`
Tree-of-Blobs inference results and related analyses.

### `scripts/`
Analysis scripts and utilities specific to the empirical study.

### `fig/`
Figures and visualizations generated from the empirical analysis
