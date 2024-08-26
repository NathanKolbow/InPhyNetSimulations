#!/bin/bash

ntaxa=$1
rep=$2
ils=$3
ngt=$4
m=$5
gt_idx=$6


iqtree="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/iqtree"

seq_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/sequences/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.nexus_${gt_idx}"

output_prefix="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/iqtree/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}_${gt_idx}"

echo "iqtree: n${ntaxa}-r${rep}-${ils}, ${ngt}gt, m=${m}, gt#${gt_idx}"
${iqtree} -nt 8 -s "${seq_file}" -o "OUTGROUP" -pre "${output_prefix}" -seed 0 -quiet