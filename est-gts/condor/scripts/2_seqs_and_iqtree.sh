#!/bin/bash

ntaxa=$1
rep=$2
ils=$3
ngt=$4
m=$5
gt_idx=$6

# Software
seqgen="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/seq-gen"
iqtree="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/iqtree"

# IQTree output prefix
iqtree_prefix="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/iqtree/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}_${gt_idx}"

# Input file path
true_gt_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.treefile_${gt_idx}"

# Output file paths
seqgen_output_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/sequences/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.nexus_${gt_idx}"
consolidated_path="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.treefile"
output_tree="${iqtree_prefix}.treefile"
cleaned_file="${iqtree_prefix}.ckp.gz"

# 1. if consolidated file already exists then we're done.
if [ -f "${consolidated_path}" ] && [ $(wc -l < "${consolidated_path}") -gt 0 ]; then
    echo "Consolidated IQTree files already exist and contain data - skipping."
    exit 0
fi

# 2. if seqgen output does not exist, run seqgen
if [ ! -f "${seqgen_output_file}" ]; then
    echo "Running seqgen"
    ${seqgen} -q -s0.036 -n1 -f0.3,0.2,0.2,0.3 -mHKY -on -l1000 ${true_gt_file} -z${gt_idx} >> "${seqgen_output_file}"
fi

# 3. unless `output_tree` is present while `cleaned_file` is not, run IQTree
if [ ! -f "${output_tree}" ] || [ -f "${cleaned_file}" ]; then
    echo "Running IQTree"
    ${iqtree} -nt 8 -s "${seqgen_output_file}" -o "OUTGROUP" -pre "${iqtree_prefix}" -seed 0 -quiet

    # 4. Clean up
    echo "Cleaning up"
    rm "${seqgen_output_file}"
    rm "${iqtree_prefix}.model.gz"
    rm "${iqtree_prefix}.bionj"
    rm "${iqtree_prefix}.ckp.gz"
    rm "${iqtree_prefix}.log"
    rm "${iqtree_prefix}.mldist"
    rm "${iqtree_prefix}.iqtree"
fi


