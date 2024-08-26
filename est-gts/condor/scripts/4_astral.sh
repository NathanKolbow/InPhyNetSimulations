#!/bin/bash

ntaxa=$1
rep=$2
ils=$3
ngt=$4
m=$5


astral_input="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/est-gts/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.treefile"
output_path="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/astral/n${ntaxa}-r${rep}-${ils}-${ngt}gt-m${m}.treefile"
astral="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/software/Astral/astral.5.7.1.jar"


java -jar ${astral} -i "${astral_input}" -o "${output_path}" -t 0 -s 0