ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])


# 0. Check if the DF already has an entry for these params
using CSV, DataFrames

df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", DataFrame)
if nrow(filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, df)) > 0
    @info "Entry already exists, quitting."
    exit()
end


using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


# 1. Figure out how many subsets this network has
truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/astral/"
astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

tre0 = readTopology(astral_file)
subsets = sateIdecomp(tre0, m)

n_subsets = length(subsets)


# 2. Load data, all the while making sure everything actually exists
snaq_newicks = Array{String}(undef, n_subsets, 10)
runtimes = Array{Float64}(undef, n_subsets, 10)
nlls = Array{Float64}(undef, n_subsets, 10)

runtimes .= 0.0
nlls .= Inf

for subset_idx = 1:n_subsets
    for run_number = 1:10
        snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq/"
        snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-run$(run_number)"

        # isfile("$(snaq_prefix).runtime") || error("File $(snaq_prefix).runtime does not exist!")
        isfile("$(snaq_prefix).runtime") || continue
        
        out_info = split(readlines("$(snaq_prefix).out")[1], " -Ploglik = ")
        snaq_newicks[subset_idx, run_number] = out_info[1]
        nlls[subset_idx, run_number] = parse(Float64, out_info[2])
        runtimes[subset_idx, run_number] = parse(Float64, readlines("$(snaq_prefix).runtime")[1])
    end
end


# 3. Select constraints with lowest negative log likelihood
constraints = Vector{HybridNetwork}([
    readTopology(newick) for newick in snaq_newicks[findmin(nlls, dims=2)[2]]
][:,1])


# 4. Calculate D
set_seed!(ntaxa, rep, ils, ngt, m, 0)
gts = simulatecoalescent(truenet, ngt, 1)
D, namelist = calculateAGIC(gts)


# 5. Root all the constraints
for c in constraints
    while length(c.node[c.root].edge) > 2
        rootatnode!(c, getchildren(c.node[c.root])[1])
    end
end


# 5. Construct network with InPhyNet
inphynet_runtime = @elapsed mnet = inphynet(D, constraints, namelist)


# 6. Check and save results
rootatnode!(mnet, "OUTGROUP")
rootatnode!(truenet, "OUTGROUP")

major_hwcd = hardwiredClusterDistance(majorTree(mnet), majorTree(truenet), true)
if major_hwcd != 0
    println("Calculating major unrooted HWCD:")
    major_hwcd, _ = hardwiredClusterDistance_unrooted_mine(majorTree(mnet), majorTree(truenet))
end

hwcd = hardwiredClusterDistance(mnet, truenet, true)
if hwcd != 0
    println("Calculating unrooted HWCD:")
    hwcd, unrooted_hwcd_roots = hardwiredClusterDistance_unrooted_mine(mnet, truenet)

    rootatnode!(mnet, unrooted_hwcd_roots[1])
    rootatnode!(truenet, unrooted_hwcd_roots[2])
end

nretic_true = truenet.numHybrids
nretic_est = mnet.numHybrids

println("Calculating min retic subset HWCD:")
min_hwcd, min_hwcd_truenet = find_minimum_retic_subset_hwcd(truenet, mnet, verbose=true)

snaq_runtime_sum = sum(runtimes)
snaq_runtime_serial = maximum(runtimes)

pruned_constraints = pruneTruthFromDecomp(truenet, [tipLabels(c) for c in constraints])
sum_constraint_hwcd = sum(hardwiredClusterDistance(pruned_constraints[i], constraints[i], false) for i = 1:length(constraints))

@show sum_constraint_hwcd
@show hwcd
@show min_hwcd
@show major_hwcd
@show nretic_true
@show nretic_est
@show snaq_runtime_sum
@show snaq_runtime_serial

if any(nlls .== Inf)
    error("Some entries not finished, quitting before data is saved.")
    exit()
end

CSV.write("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/out.csv", 
DataFrame(
    ntaxa=ntaxa, rep=rep, ils=ils, ngt=ngt, m=m,
    snaq_runtime_sum=snaq_runtime_sum,
    snaq_runtime_serial=snaq_runtime_serial,
    inphynet_runtime=inphynet_runtime,
    nretic_true=nretic_true, nretic_est=nretic_est,
    hwcd=hwcd, major_hwcd=major_hwcd
), append=true)



