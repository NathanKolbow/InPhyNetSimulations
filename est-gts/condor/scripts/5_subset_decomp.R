args <- commandArgs(trailingOnly = TRUE)
ntaxa <- as.integer(args[1])
rep <- as.integer(args[2])
ils <- args[3]
ngt <- as.integer(args[4])
m <- as.integer(args[5])

# File paths
astral_file <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/astral/"
astral_file <- paste0(astral_file, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
est_gt_file <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/"
est_gt_file <- paste0(est_gt_file, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
subset_folder <- "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/"
subset_folder <- paste0(subset_folder, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")

# If we've already done this iteration, quit
quit()
if(file.exists(paste0(subset_folder, ".nsubset"))) {
    quit()
}

# Load packages
library(ape)
library(MSCquartets)

# Load 