#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"
subscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${subscriptdir}/../../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/"
tempdir="${basedir}/temp-data/"


# 1. Trim gene trees down to the subset and figure out numRetics
export JULIA_DEPOT_PATH="${subscriptdir}/../../"
julia "${subscriptdir}/phylonet_preproc.jl" "${basedir}/estgts.tre" "${basedir}/subsets" "${basedir}/temp-data/" "${basedir}/true.net"

# 2. Run each nexus file w/ the phylonet java jar
rm ${basedir}/phylonet.runtime
rm ${basedir}/phylonet.net
rm ${basedir}/phylonet.runtime-temp
rm ${basedir}/phylonet.net-temp

find ${basedir}/temp-data/phylonet*.nexus -type f -print0 |
while IFS= read -r -d '' FILE; do
    echo "Running PhyloNet on ${FILE}"
    start=$(date +%s.%N)
    java -jar "${subscriptdir}/../../software/PhyloNetv3_8_2.jar" "${FILE}"
    end=$(date +%s.%N)
    dur=$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.6f", e - s }')
    echo ${dur} >> ${basedir}/phylonet.runtime-temp
done

# 3. Concatenate output into a single constraints file
rm ${basedir}/phylonet.net
find ${basedir}/temp-data/phylonet*.net -print0 |
while IFS= read -r -d '' FILE; do
    sed -n '3p' "$FILE" >> ${basedir}/phylonet.net-temp
done

mv ${basedir}/phylonet.runtime-temp ${basedir}/phylonet.runtime
mv ${basedir}/phylonet.net-temp ${basedir}/phylonet.net

