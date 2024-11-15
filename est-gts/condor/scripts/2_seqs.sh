#!/bin/bash

ntaxa=$1
rep=$2
ils=$3
ngt=$4
m=$5

seqgen="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/seq-gen"
consolidated_path="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.treefile"

if [ -f "${consolidated_path}" ]; then
    echo "Consolidated IQTree files already exist for this - skipping."
    exit 0
fi


for i in $(seq 1 ${ngt})
do
    echo -ne "\r${i}/${ngt}"
    true_gt_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}/${i}.treefile"

    seq_file="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/sequences/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.nexus_${i}"

    ${seqgen} -q -s0.036 -n1 -f0.3,0.2,0.2,0.3 -mHKY -on -l1000 ${true_gt_file} -z${i} >> "${seq_file}"
done
echo ""