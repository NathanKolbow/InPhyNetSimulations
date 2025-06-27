#!/usr/bin/env Rscript

# Parse args
args <- commandArgs(trailingOnly = TRUE)
subset_dir <- args[1]
which_subset <- args[2]         # string: "clade" or "ancestral"
index <- as.integer(args[3])


# subset_dir <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/"
# subset_dir <- paste0(subset_dir, "n", ntaxa, "-r", rep, "-", ils, "-", ngt, "gt-m", m, "/")

# Quit if data exists
if(file.exists(paste0(subset_dir, "/", which_subset, index, ".tob"))) {
    cat("Already exists -- quitting.\n")
    quit()
}

# Load packages
library(ape)
library(MSCquartets)

# Read respective gene trees
gts <- read.tree(file=paste0(subset_dir, "/", which_subset, index, ".tre"))
qtab <- quartetTableParallel(gts, numCores=2)

# Infer tree of blobs
tt <- TINNIK(qtab, alpha = 0.01, beta = 0.99, plot = FALSE)

# Save tree of blobs
write.tree(tt$ToB, file = paste0(subset_dir, "/", which_subset, index, ".tob"))
# save.image(file = paste0(subset_dir, which_subset, index, ".tob.RData"))