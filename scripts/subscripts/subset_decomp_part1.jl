using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
# Pkg.instantiate()
# Pkg.update()

# Parse ARGS
est_gt_file = ARGS[1]
@info est_gt_file
tempdir = ARGS[2]
k = 50

# File paths
subset_dir = joinpath(tempdir, "subsets")
if !isdir(subset_dir) mkdir(subset_dir) end

# Load packages
using PhyloNetworks, InPhyNet, DataFrames, StatsBase
import InPhyNet: prune_network

# Get NJ tree
est_gts = readMultiTopology(est_gt_file)
@info length(est_gts)
D, namelist = calculateAGID(est_gts)
@info namelist
nj_tre = inphynet(D, Vector{HybridNetwork}([]), namelist)   # PhyloNetworks.nj is SUPER slow for large `n`
k = min(k, length(namelist))
cladewiseorder!(nj_tre)
ordered_leaves = [node.name for node in nj_tre.node[nj_tre.vec_int1] if node.leaf && node.name != "OUTGROUP"]

# Perform "broad" subset decomp w/ NJ tree
clade_subsets = Vector{Vector{String}}([])
ancestral_subsets = Vector{Vector{String}}([])
n_subsets = max(1, length(ordered_leaves) ÷ k)

# Quit if data already exists
if isfile(joinpath(subset_dir, "nsubsets")) && length(readlines(joinpath(subset_dir, "nsubsets"))) >= 2*n_subsets exit() end

for j = 1:n_subsets
    start_idx = (j-1)*k+1
    end_idx = j*k
    push!(clade_subsets, ordered_leaves[start_idx:end_idx])
    clade_subsets[j] = sample(clade_subsets[j], k, replace=false)
end

for j = 1:n_subsets
    push!(ancestral_subsets, [])
    for i = 1:k
        push!(ancestral_subsets[j], ordered_leaves[(i-1)*n_subsets+j])
    end 
end

# Write the pruned gene trees to files
function smart_prune(net::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    if length(subsets) == 1
        return [prune_network(net, subsets[1])]
    elseif length(subsets) <= 4
        return prune_network(net, subsets)
    else
        i = length(subsets) ÷ 2
        split1 = subsets[1:i]
        split2 = subsets[(i+1):length(subsets)]
        n1 = prune_network(net, reduce(vcat, split1))
        n2 = prune_network(net, reduce(vcat, split2))
        return [smart_prune(n1, split1); smart_prune(n2, split2)]
    end
end



pruned_gts = Array{HybridNetwork}(undef, length(est_gts), n_subsets)
Threads.@threads for j = 1:length(est_gts)
    pruned_gts[j, :] .= smart_prune(est_gts[j], clade_subsets)
end
for i = 1:n_subsets
    writeMultiTopology(pruned_gts[:,i], joinpath(subset_dir, "clade$(i).tre"))
end

Threads.@threads for j = 1:length(est_gts)
    pruned_gts[j, :] .= smart_prune(est_gts[j], ancestral_subsets)
end
for i = 1:n_subsets
    writeMultiTopology(pruned_gts[:,i], joinpath(subset_dir, "ancestral$(i).tre"))
end

# Write the `nsubsets` file
open(joinpath(subset_dir, "nsubsets"), "w+") do f
    for i = 1:n_subsets
        for which_subset in ["clade", "ancestral"]
            write(f, "$(which_subset),$(i)\n")
        end
    end
end