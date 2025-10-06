# Simulation Scripts

This directory contains the main workflow scripts that orchestrate the InPhyNet simulation study. These scripts implement the complete pipeline from true network generation through final performance evaluation.

## Main Pipeline Scripts

### Core Simulation Workflow
- **`perform_simulation.sh`**: Master script coordinating the entire simulation pipeline

### Method-Specific Scripts  
- **`perform_snaq.sh`**: SNaQ simulations for large (200-taxa) networks
- **`perform_phylonet_simulations.sh`**: PhyloNet constraint inference workflow
- **`perform_phylonet-ml_simulations.sh`**: Maximum likelihood PhyloNet analysis
- **`perform_squirrel_simulations.sh`**: Squirrel method simulations

### Individual Pipeline Steps
- **`generate_true_network.R`**: Generates ground truth species networks with reticulation
- **`simulate_empirical_data.sh`**: Simulates gene trees, sequences, and estimates gene trees
- **`subset_decomposition.sh`**: Performs taxon subset decomposition for divide-and-conquer
- **`infer_constraints.sh`**: Runs constraint network inference on subsets
- **`generate_empirical_data.sh`**: Alternative empirical data generation workflow

## Subscripts Directory
The `subscripts/` subdirectory contains specialized utility functions and helper scripts called by the main pipeline scripts. The details of these scripts are not enumerated here.

## Simulation Pipeline Overview

### 1. True Network Generation
- Simulates species networks using birth-death-hybridization processes
- Controls reticulation frequency and network topology
- Outputs networks in extended Newick format

### 2. Data Simulation  
- **Gene tree simulation**: Multi-species coalescent on species networks
- **Sequence simulation**: Molecular evolution along gene trees using Seq-Gen
- **Gene tree estimation**: Maximum likelihood inference using IQ-TREE

### 3. Subset Decomposition
- **Small networks (≤50 taxa)**: Direct TINNIK-based decomposition
- **Large networks (>50 taxa)**: Hierarchical decomposition with ASTRAL species trees prior to TINNIK-based decomposition
- Creates non-overlapping subsets for constraint inference

### 4. Constraint Inference
- Applies SNaQ, PhyloNet, or Squirrel to each subset
- Generates partial network constraints
- Records runtime and memory usage

### 5. Network Merging
- Uses InPhyNet algorithm to merge constraint networks
- Produces final species network estimate

## Parameter Combinations

Scripts systematically enumerate combinations of:
- **Network size**: 25, 50, 100, 200 taxa
- **Gene trees**: 100, 1000 trees
- **ILS level**: Low, high incomplete lineage sorting
- **Sequence length**: 100, 1000 base pairs
- **Subset size**: Maximum 10 or 20 taxa per subset (only 10 for PhyloNet-ML)
- **Replicates**: 10 independent runs per parameter set
- **Methods**: SNaQ, PhyloNet-MPL, PhyloNet-ML, Squirrel constraint inference

## Execution

The scripts are designed for distributed computing via HTCondor. They include:
- Automatic dependency checking
- Intermediate result caching
- Error handling and logging
