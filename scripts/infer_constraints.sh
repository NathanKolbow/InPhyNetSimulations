#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/${imethod}/"
python="${scriptdir}/../python-venv/bin/python3.13"

#--------------------
# SNaQ
#--------------------



#--------------------
# PhyloNet
#--------------------
# Tutorial: https://phylogenomics.rice.edu/html/commands/InferNetwork_MPL.html


#--------------------
# Squirrel
#--------------------
# Tutorial: https://github.com/nholtgrefe/squirrel/tree/main/physquirrel
if [[ "${imethod}" == "squirrel" ]]; then
    if [[ -f "${basedir}/constraints.net" ]]; then
        rm "${basedir}/constraints.net"
    fi
    ${python} "${scriptdir}/subscripts/infer_squirrel_constraints.py" "${basedir}/subsets" "${basedir}/msa.fasta" "${basedir}/constraints.net" "${basedir}/constraints.runtime" "${basedir}/temp-data/"
fi
