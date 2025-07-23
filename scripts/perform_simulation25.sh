#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for r in $(seq 1 10); do
for m in 10; do
for ngt in 100 1000; do
for ils in low high; do
for nbp in 100 1000; do
for n in 25; do
for imethod in snaq squirrel phylonet phylonet-ml; do
    printf "\033[0;31m-----------------------------------------------------------------------------------\n"
    printf "\033[0;31m>> ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r $imethod <<\n"
    printf "\033[0;31m-----------------------------------------------------------------------------------\n\n\033[0m"
    ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r $imethod
done
done
done
done
done
done
done
