#!/usr/bin/env Rscript

# Parse args
args <- commandArgs(trailingOnly = TRUE)
ntaxa <- as.integer(args[1])
rep <- as.integer(args[2])
ils <- args[3]
ngt <- as.integer(args[4])
m <- as.integer(args[5])
which_subset <- args[6]         # string: "clade" or "ancestral"
index <- as.integer(args[7])
ncores <- as.integer(args[8])   # 5 or 6 seems to be where efficiency stops improving, so call it an even 4

subset_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/"
subset_dir <- paste0(subset_dir, "n", ntaxa, "-r", rep, "-", ils, "-", ngt, "gt-m", m, "/")

# Quit if data exists
if(file.exists(paste0(subset_dir, which_subset, index, ".tob"))) {
    cat("Already exists -- quitting.\n")
    quit()
}

# Load packages
library(ape)
library(MSCquartets)

# Read respective gene trees
gts <- read.tree(file=paste0(subset_dir, which_subset, index, ".tre"))
qtab <- quartetTableParallel(gts, numCores=ncores)

# Infer tree of blobs
tt <- TINNIK(qtab, alpha = 0.01, beta = 0.99, plot = FALSE)

# Save tree of blobs
write.tree(tt$ToB, file = paste0(subset_dir, which_subset, index, ".tob"))
# save.image(file = paste0(subset_dir, which_subset, index, ".tob.RData"))
