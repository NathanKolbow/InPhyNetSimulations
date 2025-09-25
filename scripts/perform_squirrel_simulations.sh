#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for n in 25 50 100 200; do
for m in 10 20; do
for r in $(seq 1 10); do
for ils in low high; do
for ngt in 100 1000; do
for nbp in 100 1000; do
    #if [[ $m -eq 20 && $n -eq 25 ]]; then
    #    break
    #fi
    #printf "\033[0;31m---------------------------------------------------------------------------------------\n"
    printf "\033[0;31m>> ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r squirrel <<\n"
    #printf "\033[0;31m---------------------------------------------------------------------------------------\n\n\033[0m"
    ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r squirrel
done
done
done
done
done
done
