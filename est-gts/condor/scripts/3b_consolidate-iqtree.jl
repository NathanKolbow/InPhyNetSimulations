using PhyloNetworks


ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])


iqtree_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/iqtree/"
consolidated_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/simulation-data/est-gts/"
consolidated_path *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

trees = Array{HybridNetwork}(undef, ngt)
for j = 1:ngt
    trees[j] = readTopology(joinpath(
        iqtree_path,
        "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)_$(j).treefile"
    ))
end

missing_idxs = findall(i -> !isassigned(trees, i), 1:ngt)
if length(missing_idxs) != 0
    throw(ErrorException("Missing $(length(missing_idxs)) estimated gene trees!"))
end

writeMultiTopology(trees, consolidated_path)