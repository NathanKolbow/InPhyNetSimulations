#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#for r in $(seq 1 10); do
for r in $(seq 1 10); do
for m in 10 20; do
for ngt in 100 1000; do
for ils in low high; do
for nbp in 1000 100; do
for n in 25 50 100 200; do
for imethod in snaq phylonet phylonet-ml squirrel; do
    # printf "\033[0;31m-----------------------------------------------------------------------------------\n"
    # printf "\033[0;31m>> ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r snaq <<\n"
    # printf "\033[0;31m-----------------------------------------------------------------------------------\n\n\033[0m"
    ${scriptdir}/subscripts/perform_inphynet.sh $n $ngt $ils $nbp $m $r $imethod
done
done
done
done
done
done
done
