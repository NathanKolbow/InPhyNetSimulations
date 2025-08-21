#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#for r in $(seq 1 10); do
for r in 1; do
for m in 20 10; do
for ngt in 100 1000; do
for ils in low high; do
for nbp in 1000 100; do
for n in 50 100 200; do
    printf "\033[0;31m-----------------------------------------------------------------------------------\n"
    printf "\033[0;31m>> ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r snaq <<\n"
    printf "\033[0;31m-----------------------------------------------------------------------------------\n\n\033[0m"
    ${scriptdir}/perform_simulation.sh $n $ngt $ils $nbp $m $r snaq
done
done
done
done
done
done
