# Read arguments
group = ARGS[1]
subset_number = parse(Int64, ARGS[2])
h = parse(Int64, ARGS[3])
run_number = parse(Int64, ARGS[4])


# Put together paths
subset_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/empirical-study/subset_data/$(group)/subset$(subset_number)"
snaq_prefix = "$(subset_dir)/snaq_s$(subset_number)_h$(h)"


# Check if we already completed this w/ the setup back when we were doing all 10 runs at once
if isfile("$(snaq_prefix).out") && length(readlines("$(snaq_prefix).out")) > 0 && occursin("Elapsed time", readlines("$(snaq_prefix).out")[3])
    @info "Already completed w/ all 10 runs."
    exit()
end

# Check if we already completed this w/ the new setup where each run is done individually
snaq_prefix = "$(subset_dir)/snaq_s$(subset_number)_h$(h)_run$(run_number)"
if isfile("$(snaq_prefix).out") && length(readlines("$(snaq_prefix).out")) > 0 && occursin("Elapsed time", readlines("$(snaq_prefix).out")[3])
    @info "Already completed this individual run."
    exit()
end

# Packages
# using Pkg
# Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/precompile-setup.jl")
using PhyloNetworks


# Read data
tre0 = readTopology("$(subset_dir)/best_species_tree.treefile")
est_gts = readMultiTopology("$(subset_dir)/estimated_gts.treefile")


q, t = countquartetsintrees(est_gts)
df = readTableCF(writeTableCF(q, t))

# Make output dir if it doesn't exist
if !isdir("$(subset_dir)") mkdir("$(subset_dir)") end
if !isfile("$(snaq_prefix).out") touch("$(snaq_prefix).out") end
if !isfile("$(snaq_prefix).log") touch("$(snaq_prefix).log") end
if !isfile("$(snaq_prefix).err") touch("$(snaq_prefix).err") end


snaq!(tre0, df, hmax = h, filename = snaq_prefix, runs = 1, probST = 0.5, seed = run_number)
