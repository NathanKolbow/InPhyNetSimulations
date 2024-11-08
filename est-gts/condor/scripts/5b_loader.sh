#!/bin/bash

# R-4.4.0 module
export PATH="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/bin:$PATH"
export MANPATH="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/share/man:$MANPATH"
export CMAKE_PREFIX_PATH="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/.:$CMAKE_PREFIX_PATH"
export LD_LIBRARY_PATH="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/lib64/R/lib:$LD_LIBRARY_PATH"
export PKG_CONFIG_PATH="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/lib64/pkgconfig:$PKG_CONFIG_PATH"
export R_HOME="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/lib64/R"
export R_ROOT="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0"

# gcc-13.2.0 module
root="/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/gcc-13.2.0"
export CMAKE_LIBRARY_PATH="$root/lib64:$CMAKE_LIBRARY_PATH"
export CMAKE_PREFIX_PATH="$root:$CMAKE_PREFIX_PATH"
export LD_LIBRARY_PATH="$root/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$root/lib64:$LD_LIBRARY_PATH"
export MANPATH="$root/share/man:$MANPATH"
export PATH="$root/bin:$PATH"
export XDG_DATA_DIRS="$root/share:$XDG_DATA_DIRS"


Rscript $1 $2 $3 $4 $5 $6 $7 $8 $9