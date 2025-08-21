#!/bin/bash

. /etc/profile.d/modules.sh

# Move to project dir
echo "Moving directories"
cd /mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNetSimulations/

# Setup Julia
echo "Exporting Julia path"
export PATH="$PWD/software/julia-1.11.6/bin:$PATH"
echo "Exporting Julia depot path"
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR
echo "Using Julia: `which julia`"

# Setup R
echo "Loading R-4.4.0"
module load R-4.4.0
echo "Exporting R_LIBS_PATH"
export R_LIBS_PATH="/mnt/home/nkolbow/R/x86_64-pc-linux-gnu-library/4.4"

# Run simulation
echo "Running . ./scripts/perform_simulation.sh $1 $2 $3 $4 $5 $6 $7"
. ./scripts/perform_simulation.sh $1 $2 $3 $4 $5 $6 $7
