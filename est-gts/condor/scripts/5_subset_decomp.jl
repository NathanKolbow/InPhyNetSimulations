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

@info "Reading trees"
astral_tree = readTopology(astral_file)
est_gts = readMultiTopology(est_gt_file)

@info "Decomposing"
subsets = sateIdecomp(astral_tree, m)

for (i, subset) in enumerate(subsets)
    @info "$(i)/$(length(subsets))"
    iter_folder = "$(subset_folder)-subset$(i)/"

    if isdir(iter_folder) && isfile(joinpath(iter_folder, "pruned_gts.treefile")) && length(readlines(joinpath(iter_folder, "pruned_gts.treefile"))) == length(est_gts)
        printstyled("\t] ALREADY COMPLETED\n", color=:red)
        continue
    end

    tre0 = pruneTruthFromDecomp(astral_tree, subset)
    iter_gts = [pruneTruthFromDecomp(gt, subset) for gt in est_gts]
    
    if !isdir(iter_folder) mkdir(iter_folder) end
    writeTopology(tre0, joinpath(iter_folder, "tre0.treefile"))
    writeMultiTopology(iter_gts, joinpath(iter_folder, "pruned_gts.treefile"))
end
open("$(subset_folder).nsubsets", "w+") do f
    for j = 1:length(subsets)
        for run_idx = 1:10
            write(f, "$(j),$(run_idx)\n")
        end
    end
end

printstyled("[COMPLETE] ", color = :green)
println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m), $(length(subsets)) subsets")
