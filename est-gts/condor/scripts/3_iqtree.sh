#!/bin/bash

ntaxa=$1
rep=$2
ils=$3
ngt=$4
m=$5
gt_idx=$6


iqtree="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/iqtree"

seq_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/sequences/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.nexus_${gt_idx}"

output_prefix="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/iqtree/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}_${gt_idx}"


output_tree="${output_prefix}.treefile"
cleaned_file="${output_prefix}.ckp.gz"

if [ -f "${output_tree}" ] && [ ! -f "${cleaned_file}" ]; then
    echo "${output_tree} exists but ${cleaned_file} does not - exiting."
    exit 0
fi



echo "iqtree: n${ntaxa}-r${rep}-${ils}, ${ngt}gt, m=${m}, gt#${gt_idx}"
${iqtree} -nt 8 -s "${seq_file}" -o "OUTGROUP" -pre "${output_prefix}" -seed 0 -quiet

# Clean up
rm "${output_prefix}.model.gz"
rm "${output_prefix}.bionj"
rm "${output_prefix}.ckp.gz"
rm "${output_prefix}.log"
rm "${output_prefix}.mldist"
rm "${output_prefix}.iqtree"
