#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"
subscriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${subscriptdir}/../../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/"


# 1. Trim gene trees down to the subset and figure out numRetics
export JULIA_DEPOT_PATH="${subscriptdir}/../../"
julia "${subscriptdir}/phylonet_preproc.jl" "${basedir}/estgts.tre" "${basedir}/subsets" "${basedir}/temp-data/" "${basedir}/true.net"

# 2. Run each nexus file w/ the phylonet java jar
{ time find ${basedir}/temp-data/phylonet*.nexus -type f -exec java -jar "$subscriptdir/../../software/PhyloNetv3_8_2.jar" {} \; ; } 2>> "${basedir}/phylonet.runtime"

# 3. Concatenate output into a single constraints file
find ${basedir}/temp-data/phylonet*.net \
    -type f \
    -exec cat {} >> "${basedir}/phylonet.net"
