# HTCondor Cluster Configuration

This directory contains configuration files and utilities for running InPhyNet simulations on an HTCondor distributed computing cluster. HTCondor enables parallelization of computationally intensive phylogenetic network inference across multiple machines.

## Key Files

- **`driver.sh`**: Main execution script that runs simulation jobs on cluster nodes
- **`submit`**: HTCondor job submission file with resource requirements and job parameters
- **`submit.table`**: Table of parameter combinations for batch job submission
- **`submit.table.jl`**: Julia script to generate parameter combination table

## Directories

- **`err/`**: Standard error output files from cluster jobs
- **`log/`**: HTCondor log files tracking job status and resource usage  
- **`out/`**: Standard output files from completed jobs

## Usage

HTCondor was primarily used for SNaQ and PhyloNet simulations, which are computationally expensive. Squirrel simulations were run locally due to their efficiency.

For more details about HTCondor, visit [htcondor.org](https://htcondor.org/).