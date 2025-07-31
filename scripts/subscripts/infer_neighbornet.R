library(phangorn)
library(ape)

args <- commandArgs(trailingOnly = TRUE)
msa_path <- args[1]
output_path <- args[2]

alignment <- read.dna(msa_path, format="fasta")
d <- dist.ml(alignment)
net <- neighborNet(d)

write.nexus.splits(output_path)
