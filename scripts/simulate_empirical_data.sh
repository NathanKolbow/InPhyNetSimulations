#!/bin/bash

n="$1"; ngt="$2"; ils="$3"; nbp="$4"; m="$5"; r="$6"; imethod="$7"; seed="$8"

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="${scriptdir}/../data/${n}/${ngt}/${ils}/${nbp}/${m}/${r}/${imethod}/"


#------------------------------------------------------------------------------
# 1. Simulate true gene trees
#------------------------------------------------------------------------------
