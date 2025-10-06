# InPhyNetSimulations

This repository contains the simulation study and empirical analysis for InPhyNet, a novel statistical method for inferring large phylogenetic networks. InPhyNet addresses the computational challenges of inferring species networks with hundreds of taxa by using a divide-and-conquer approach that decomposes the problem into smaller, manageable subproblems.

## Overview

Phylogenetic networks are mathematical models that represent evolutionary relationships among species, including both divergence (speciation) and reticulation (hybridization, horizontal gene transfer, introgression) events. Traditional network inference methods become computationally intractable for large datasets (>30 taxa). InPhyNet overcomes this limitation through:

1. **Subset decomposition**: Breaking the full taxon set into smaller, non-overlapping subsets
2. **Constraint network inference**: Inferring partial networks on each subset using existing methods
3. **Network merging**: Combining constraint networks into a comprehensive species network

## Repository Structure

- **`condor/`**: HTCondor cluster configuration files for distributed computing
- **`data/`**: Simulation results organized by parameter combinations (ntaxa/ngt/ILS/nbp/m/replicate/method)
- **`empirical-study/`**: Real data analysis on primate datasets
- **`figs/`**: R scripts for generating publication figures and result visualization
- **`helpers/`**: Julia utility functions for data simulation, analysis, and plotting
- **`scripts/`**: Main simulation pipeline scripts and workflow orchestration
- **`software/`**: External phylogenetic software dependencies (PhyloNet, ASTRAL, IQ-TREE, etc.)

## How to replicate simulations

Below are the set of parameter combinations that are enumerated in the original simulation experiment. To replicate the results of a given parameter combination, run the script `scripts/perform_simulation.sh <ntaxa> <ngt> <ILS> <nbp> <m> <rep number> <inference method>`. If data already exists from intermediate steps (e.g. estimated gene trees are already present), previous steps will be skipped. To clean the repo of all such intermediary data, run `scripts/remove_intermediary_data.sh`.

Original simulations were performed with SNaQ and PhyloNet simulations being performed on an HT Condor cluster, whereas all Squirrel simulations were performed on a single machine. This is because the Squirrel simulations are exceptionally quick to run, whereas the SNaQ and PhyloNet simulations are exceptionally time consuming.

## Simulation parameters

- [`n`] Number of taxa in the true network (25, 50, 100, 200)
- [`ngt`] Number of gene trees (100, 1000)
- [`ils`] ILS level (low, high)
- [`nbp`] MSA number of base pairs (100, 1000)
- [`m`] Maximum subset size (10, 20 (only 10 for PhyloNet-ML))
- [`r`] Replicate number (1-10)
- [`imethod`] Constraint network inference software (SNaQ, PhyloNet-MPL, PhyloNet-ML, Squirrel)

Each unique combation of parameters is used to generate a unique seed that is used as input to each algorithm where applicable. When software requires multiple seeds (e.g. we need to generate `ngt` MSAs--using the same seed for each would not yield productive results) we use this value `+i-1`, where `i` is the `i`th iteration of the software requiring a seed (e.g. the first MSA is generated with `seed`, then with `seed+1`, then `seed+2`, and so on). The function used to generate these random seeds is `generate_seed`, located in `scripts/perform_simualtion.sh`.

## Simulation pipeline breakdown

Base directory (`basedir`) for the output data of the simulations corresponding to each unique combination of parameters:

```
data/n/ngt/ils/nbp/m/r/
e.g.: data/50/100/low/100/10/1/
```

1. Ground truth network generation
    - Script: `scripts/generate_true_network.R`
    - Input: `number of taxa`, `replicate number`
    - Output: ground truth species networks are placed in `basedir/true.net`
2. Empirical data simulation (true gene tree --> MSA --> estimated gene tree)
    - Script: `scripts/simulate_empirical_data.sh`
    - Input: `number of taxa`, `replicate number`, `number of gene trees`, `ILS level`, `sequence length`, `inference method`
    - Output: for every method except Squirrel, estimated gene trees are placed in `basedir/estgts.tre` - for all methods the MSAs are placed in `basedir/msa.fasta` and the true gene trees are placed in `basedir/truegts.tre`
    - Previously generated data used: ground truth network
3. Subset decomposition
    - Script: `scripts/subset_decomposition.sh`
    - Input: `number of taxa`, `replicate number`, `number of gene trees`, `ILS level`, `sequence length`
    - Output: final subset decomposition is placed in `basedir/subsets`, a plaintext file where the file consists of `k` lines of comma-separated taxon names where each line is a subset
    - Previously generated data used: estimated gene trees
    - Steps if `n=50`:
      - Run TINNIK on the estimated gene tree data
      - Perform edge centroid decomposition while restricting subsets so that any taxa that belong to the same blob must stick together
    - Steps if `n>50`:
      - Run ASTRAL to infer a species tree
      - Run the above TINNIK step with subsets of taxa of size 50. These subsets consist of:
        - Sets of 50 taxa that are close to one another, i.e. if the true network were a completely unbalanced tree the decomposition would look like: `subset1 = {t1,t2,t3,...,t49,t50}`, `subset2 = {t51,t52,...,t100}`, ...
        - Sets of 50 taxa that are equidistant from one another, i.e. again if the true network were an unbalanced tree with 200 taxa, the decomposition would look like: `subset1 = {t1,t5,t9,...,t193,t197}`, `subset2 = {t2,t6,t10,...,t194,t198}`, ..., `subset4 = {t4,t8,t12,...,t196,t200}`.
4. Constraint network inference
    - Script: `scripts/infer_constraints.sh`
    - Input: all input from step 3 AND `inference method`
    - Output: constraint networks are output to `basedir/constraints-<method>.net`, total runtime output to `basedir/constraints-<method>.runtime`
    - Previously generated data used: estimated gene trees, subset decomposition
5. Full network construction
    - Script: `scripts/run_inphynet.jl`
    - Input: all parameter values
    - Output: final network output at `basedir/inphynet-<method>.net`, total runtime output to `basedir/inphynet-<method>.runtime`
    - Previously generated data used: constraint networks, estimated gene trees

