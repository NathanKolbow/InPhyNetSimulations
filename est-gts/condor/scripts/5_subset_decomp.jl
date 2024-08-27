using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])


astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/astral/"
astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

est_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/"
est_gt_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

subset_folder = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/"
subset_folder *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)"

astral_tree = readTopology(astral_file)
est_gts = readMultiTopology(est_gt_file)

subsets = sateIdecomp(astral_tree, m)

for (i, subset) in enumerate(subsets)
    tre0 = pruneTruthFromDecomp(astral_tree, subset)
    iter_gts = [pruneTruthFromDecomp(gt, subset) for gt in est_gts]
    iter_folder = "$(subset_folder)-subset$(i)/"
    
    mkdir(iter_folder)
    writeTopology(tre0, joinpath(iter_folder, "tre0.treefile"))
    writeMultiTopology(iter_gts, joinpath(iter_folder, "pruned_gts.treefile"))
end