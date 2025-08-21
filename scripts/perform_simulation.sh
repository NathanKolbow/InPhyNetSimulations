#!/usr/bin/env bash
set -euo pipefail

#------------------------------------------------------------------------------
# 1. Parse & validate args
#------------------------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 n ngt ils nbp m r imethod
  n       : 25, 50, 100, or 200
  ngt     : 100 or 1000
  ils     : "low" or "high"
  nbp     : 100 or 1000
  m       : 10 or 20
  r       : 1-10
  imethod : "snaq" or "squirrel" or "phylonet"
EOF
  exit 1
}
[[ $# -eq 7 ]] || usage
n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"


#------------------------------------------------------------------------------
# 2. Checkpoint files for each step
#------------------------------------------------------------------------------
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/"
mkdir -p "${basedir}"
mkdir -p "${basedir}/temp-data/"

step1_outputs=("${basedir}/true.net")
step2_outputs=("${basedir}/estgts.tre")
step3_outputs=("${basedir}/subsets")
  step4_outputs=("${basedir}/${imethod}.net" "${basedir}/${imethod}.runtime")
  step5_outputs=("${basedir}/inphynet-${imethod}.net" "${basedir}/inphynet-${imethod}.runtime")


#------------------------------------------------------------------------------
# 3. Find the first step that's incomplete
#------------------------------------------------------------------------------
outputs_exist() {
  local arr_name="$1"
  local -n files=$arr_name
  for f in "${files[@]}"; do
    [[ -e "$f" ]] || return 1
  done
  return 0
}

start_step=6
for i in {1..5}; do
  if ! outputs_exist "step${i}_outputs"; then
    start_step=$i
    break
  fi
done

if (( start_step > 5 )); then
  echo "All steps 1–5 appear complete. Nothing to do."
  exit 0
else
  echo "> Starting from step $start_step"
fi


generate_seed() {
  local n="$1" ngt="$2" ils="$3" nbp="$4" m="$5" r="$6" imethod="$7"
  if [ $n -eq 1000 ]; then
    echo "$r"
  else
    # 1. build a string key out of all seven args
    local key="${n}_${ngt}_${ils}_${nbp}_${m}_${r}_${imethod}"
    # 2. hash it (md5)
    local hash
    if command -v md5sum >/dev/null 2>&1; then
      hash=$(printf "%s" "$key" | md5sum | cut -d' ' -f1)
    elif command -v md5 >/dev/null 2>&1; then
      hash=$(printf "%s" "$key" | md5 -q)
    else
      echo "ERROR: no md5sum or md5 available" >&2
      exit 1
    fi
    local hex=${hash:0:8}
    # 4. convert hex to decimal
    local seed=$((16#$hex))
    seed=$(( seed % 2147483647 ))
    echo "${seed}"
  fi
}
seed=`generate_seed $n $ngt $ils $nbp $m $r $imethod`
echo "> Unique seed: ${seed}"


if (( start_step <= 1 )); then
  # 1. Generate ground truth network
  echo "> Generating ground truth network."
  module load R-4.4.0
  export R_LIBS_PATH="/mnt/home/nkolbow/R/x86_64-pc-linux-gnu-library/4.4"
  Rscript "${scriptdir}/generate_true_network.R" $n $ils $seed "${basedir}/true.net"
fi


if (( start_step <= 2 )); then
  # 2. Generate empirical data
  echo "> Generating empirical data."
  "${scriptdir}/simulate_empirical_data.sh" $n $ngt $ils $nbp $m $r $imethod $seed
fi


if (( start_step <= 3 )); then
  # 3. Subset decomposition
  echo "> Performing subset decomposition."
  "${scriptdir}/subset_decomposition.sh" $n $ngt $ils $nbp $m $r $imethod $seed
fi


if (( start_step <= 4 )); then
  # 4. Infer constraint networks
  echo "> Inferring constraint networks."
  "${scriptdir}/infer_constraints.sh" $n $ngt $ils $nbp $m $r $imethod $seed
fi


if (( start_step <= 5 )); then
  # 5. Construct full network with InPhyNet
  echo "> Constructing full network with InPhyNet."
  export JULIA_DEPOT_PATH="${scriptdir}/../"
  julia "${scriptdir}/combine_inphynet.jl" "${basedir}/${imethod}.net" "${basedir}/estgts.tre" "${basedir}/inphynet-${imethod}.net" "${basedir}/inphynet-${imethod}.runtime"
fi




