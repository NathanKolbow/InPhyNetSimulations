start_time = time()

@everywhere ntaxa = $(parse(Int64, ARGS[1]))
@everywhere rep = $(parse(Int64, ARGS[2]))
@everywhere ils = $(ARGS[3])
@everywhere ngt = $(parse(Int64, ARGS[4]))
@everywhere m = $(parse(Int64, ARGS[5]))
@everywhere subset_idx = $(parse(Int64, ARGS[6]))


snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq/"
snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)"

if isfile("$(snaq_prefix).runtime")
    throw(ErrorException("$(snaq_prefix).runtime already exists, quitting!"))
    exit()
end


@info "Loading packages and helpers"
@everywhere begin
    using Pkg
    Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

    using PhyloNetworks
    include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
end

pruned_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)/pruned_gts.treefile"

tre0_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)/tre0.treefile"


# Find the true number of reticulations
@info "Loading true net"
truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
@info "Loading est gts"
est_gts = readMultiTopology(pruned_gt_file)

leaf_names = [leaf.name for leaf in est_gts[1].leaf]
n_retic = pruneTruthFromDecomp(truenet, leaf_names).numHybrids

@info "Counting quartets"
q, t = countquartetsintrees(est_gts)

@info "Converting quartets to table"
df = readTableCF(writeTableCF(q, t))

@info "Loading astral tree"
tre0 = readTopology(tre0_file)


@info "Starting SNaQ"
seed = set_seed!(ntaxa, rep, ils, ngt, m, 0)
snaq!(tre0, df, hmax = n_retic, seed = seed, filename = snaq_prefix)

runtime_seconds = time() - start_time
open("$(snaq_prefix).runtime", "w+") do f
    write(f, "$(runtime_seconds)")
end