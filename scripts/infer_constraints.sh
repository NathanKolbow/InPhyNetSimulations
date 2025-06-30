#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/"
python="${scriptdir}/../python-venv/bin/python3.13"

#--------------------
# SNaQ
#--------------------
if [[ "${imethod}" == "snaq" ]]; then
    if [[ -f "${basedir}/constraints-${imethod}.net" ]]; then
        rm "${basedir}/constraints-${imethod}.net"
    fi
    export JULIA_DEPOT_PATH="${scriptdir}/.."
    julia -t8 "${scriptdir}/subscripts/infer_snaq.jl" "${basedir}/subsets" "${basedir}/estgts.tre" "${basedir}/snaq.net" "${basedir}/snaq.runtime" "${basedir}/true.net" "${basedir}/temp-data/" "${seed}"
fi


#--------------------
# PhyloNet
#--------------------
# Tutorial: https://phylogenomics.rice.edu/html/commands/InferNetwork_MPL.html
if [[ "${imethod}" == "phylonet" ]]; then
    if [[ -f "${basedir}/constraints-${imethod}.net" ]]; then
        rm "${basedir}/constraints-${imethod}.net"
    fi
    ${scriptdir}/subscripts/infer_phylonet.sh "${basedir}/estgts.tre" "${basedir}/temp-data" "${basedir}/phylonet.net" "${basedir}/phylonet.runtime"
fi

#--------------------
# Squirrel
#--------------------
# Tutorial: https://github.com/nholtgrefe/squirrel/tree/main/physquirrel
if [[ "${imethod}" == "squirrel" ]]; then
    if [[ -f "${basedir}/constraints-${imethod}.net" ]]; then
        rm "${basedir}/constraints-${imethod}.net"
    fi
    ${python} "${scriptdir}/subscripts/infer_squirrel_constraints.py" "${basedir}/subsets" "${basedir}/msa.fasta" "${basedir}/squirrel.net" "${basedir}/squirrel.runtime" "${basedir}/temp-data/"
fi
