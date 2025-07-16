using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()
Pkg.update()

using PhyloNetworks, SNaQ, InPhyNet
subset_file = ARGS[1]
gt_file = ARGS[2]
net_output = ARGS[3]
runtime_output = ARGS[4]
truenet_file = ARGS[5]
temp_dir = ARGS[6]
seed = parse(Int, ARGS[7])

temp_net_output = "$(net_output)-incomplete"
temp_rt_output = "$(runtime_output)-incomplete"

starting_subset = 1
if isfile(temp_net_output) && !isfile(temp_rt_output)
    rm(temp_net_output)
elseif !isfile(temp_net_output) && isfile(temp_rt_output)
    rm(temp_rt_output)
elseif isfile(temp_net_output) && isfile(temp_rt_output)
    if length(readlines(temp_net_output)) != length(readlines(temp_rt_output))
        rm(temp_net_output)
        rm(temp_rt_output)
    else
        starting_subset = length(readlines(temp_net_output)) + 1
    end
end

truenet = readnewick(truenet_file)
full_gts = readmultinewick(gt_file)
subsets = readlines(subset_file)
for (isubset, subset) in enumerate(subsets)
    if isubset < starting_subset continue end
    subset = split(subset, ",")

    rt = time()
    subset_truenet = prune_network(truenet, subset)
    gts = prune_networks(full_gts, subset)
    df = readtrees2CF(gts, writeSummary=false)
    snaq_net = snaq!(gts[1], df; hmax=subset_truenet.numhybrids, runs=10, filename="$(temp_dir)/snaq$(isubset)", seed=seed+isubset)
    rt = time() - rt

    open(temp_net_output, "a+") do f
        write(f, "$(writenewick(snaq_net))\n")
    end
    open(temp_rt_output, "a+") do f
        write(f, "$(rt)\n")
    end
end
mv(temp_net_output, net_output)
mv(temp_rt_output, runtime_output)