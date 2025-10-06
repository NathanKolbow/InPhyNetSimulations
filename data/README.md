# Simulation Data

This directory contains all simulation results organized in a hierarchical structure based on simulation parameters. The InPhyNet simulation study systematically evaluates performance across different network sizes, data characteristics, and inference methods.

## Directory Structure

Data is organized as: `ntaxa/ngt/ils/nbp/m/replicate/`

### Parameters

- **`ntaxa`**: Number of taxa in the true network (25, 50, 100, 200)
- **`ngt`**: Number of gene trees (100, 1000) 
- **`ils`**: Incomplete lineage sorting level (low, high)
- **`nbp`**: Number of base pairs in molecular sequence alignments (100, 1000)
- **`m`**: Maximum subset size for decomposition (typically 10)
- **`replicate`**: Replicate number (1-10)

### Example Path

`data/50/100/low/100/10/1/snaq/` contains results for:
- 50 taxa network
- 100 gene trees
- Low ILS level
- 100 base pair sequences
- Maximum subset size 10
- Replicate 1

## File Types Within Each Parameter Combination

- **`true.net`**: Ground truth species network in extended Newick format
- **`truegts.tre`**: True gene trees simulated from the species network
- **`msa.fasta`**: Molecular sequence alignments simulated from gene trees
- **`estgts.tre`**: Estimated gene trees inferred from sequence data
- **`subsets`**: Subset decomposition (comma-separated taxon names per line)
- **`<method>.net`**: Constraint networks inferred by <method>
- **`<method>.runtime`**: Runtime information for constraint network inferred with <method>
- **`inphynet-<method>.net`**: Output from InPhyNet using constraints from <method> as input
- **`inphynet-<method>.runtime`**: Runtime information for network merging

## Summary Files

- **`all.csv`**: Aggregated results across all parameter combinations
- **`estnets.netfile`**: Collection of all estimated networks
- **`truenets.netfile`**: Collection of all true networks
