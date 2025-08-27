using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
# Pkg.instantiate()
# Pkg.update()

using PhyloNetworks, InPhyNet
gt_file = ARGS[1]
subset_file = ARGS[2]
temp_dir = ARGS[3]
truenet_file = ARGS[4]

# Turn `temp_dir` into a linux-like file string so that
# PhyloNet doesn't get confused...
phylonet_temp_dir = deepcopy(temp_dir)
if temp_dir[1] != '/'
    phylonet_temp_dir = "/" * replace(rstrip(temp_dir, '\\'), "\\" => "/", ":" => "")
end

gts = readmultinewick(gt_file)
subsets = readlines(subset_file)
truenet = readnewick(truenet_file)

for (isubset, subset) in enumerate(subsets)
    pruned_gts = prune_networks(gts, split(subset, ","))
    pruned_truenet = prune_network(truenet, split(subset, ","))
    nretic = pruned_truenet.numhybrids

    open("$(temp_dir)/phylonet$(isubset).nexus", "w+") do f
        write(f, "#NEXUS\n\nBEGIN TREES;\n\n")
        for (igt, gt) in enumerate(pruned_gts)
            write(f, "Tree gt$(igt)=$(writenewick(gt, round=true))\n")
        end
        write(f, "\nEND;\n\n\nBEGIN PHYLONET;\n\nInferNetwork_MPL (all) $(nretic) -pl 5 \"$(temp_dir)/phylonet$(isubset).net\";\n\nEND;")  # $(phylonet_temp_dir)/phylonet$(isubset).net;\n\nEND;")
    end
end
