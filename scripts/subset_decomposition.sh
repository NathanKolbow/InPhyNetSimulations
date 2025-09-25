#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/"


# Infer species tree with ASTRAL
if [[ ! -s "${basedir}/astral.tre" ]]; then
    echo "> Inferring ASTRAL tree"
    java -jar "${scriptdir}/../software/astral.5.7.1.jar" -i "${basedir}/estgts.tre" -t 0 -s ${seed} -o "${basedir}/astral.tre" &> /dev/null
else
    echo "> ASTRAL tree already exists."
fi

# Part 1
export JULIA_DEPOT_PATH="${scriptdir}/../"
if [[ ! -s "${basedir}/temp-data/subsets/nsubsets" ]]; then
    julia ${scriptdir}/subscripts/subset_decomp_part1.jl "${basedir}/estgts.tre" "${basedir}/temp-data/"
else
    echo "> Subset part 1 output already exists."
fi

# Part 2
while IFS=, read -r str num; do
    [ -z "$str" ] && continue
    if [[ -s "${basedir}/temp-data/subsets/${str}${num}.tob" ]]; then
        echo "    > ${str} ${num} already finished."
    else
        echo "    > ${str} ${num} running."
        Rscript "${scriptdir}/subscripts/subset_decomp_part2.R" "${basedir}/temp-data/subsets/" "$str" $num
    fi
done < "${basedir}/temp-data/subsets/nsubsets"

# Part 3
julia ${scriptdir}/subscripts/subset_decomp_part3.jl "${basedir}/estgts.tre" "${basedir}/temp-data/subsets/" "${basedir}/astral.tre" "${basedir}/subsets" $m
