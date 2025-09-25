#!/bin/bash

IFS="." read -ra ADDR <<< `hostname`
depot_path="/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/depot_paths/${ADDR}"

JULIA_DEPOT_PATH=${depot_path} $1 -t4 $2 $3 $4 $5 $6