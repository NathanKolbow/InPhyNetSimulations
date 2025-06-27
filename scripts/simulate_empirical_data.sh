#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/${imethod}/"


#------------------------------------------------------------------------------
# 1. Simulate true gene trees
#------------------------------------------------------------------------------
export JULIA_DEPOT_PATH="${scriptdir}/../"
julia "${scriptdir}/subscripts/simulate_true_gene_trees.jl" "${basedir}/true.net" $ngt $ils $seed "${basedir}/truegts.tre"

if [[ -f "${basedir}/msa.fasta" ]]; then
    rm "${basedir}/msa.fasta"
fi
if [[ -f "${basedir}/estgts-incomplete.tre" ]]; then
    rm "${basedir}/estgts-incomplete.tre"
fi

#------------------------------------------------------------------------------
# 2. Simulate MSAs and infer gene trees
#------------------------------------------------------------------------------
seqgen="${scriptdir}/../software/seq-gen"
iqtree="${scriptdir}/../software/iqtree3"
count=0
while IFS= read -r line; do
    ((count++))

    # 1. Simulate MSAs
    echo "$line" > "${basedir}/temp-data/temp.tre"
    ${seqgen} -n1 -of -f0.3,0.2,0.2,0.3 -mHKY -l${nbp} -q -z${seed} "${basedir}/temp-data/temp.tre" > "${basedir}/temp-data/msa${count}.fasta"
    cat "${basedir}/temp-data/msa${count}.fasta" >> "${basedir}/msa.fasta"

    # 2. Infer gene tree
    ${iqtree} -s "${basedir}/temp-data/msa${count}.fasta" -pre "${basedir}/temp-data/iqtree${count}" -seed $seed -quiet &> /dev/null
    cat "${basedir}/temp-data/iqtree${count}.treefile" >> "${basedir}/estgts-incomplete.tre"

    ((seed++))
done < "${basedir}/truegts.tre"
mv estgts-incomplete.tre estgts.tre