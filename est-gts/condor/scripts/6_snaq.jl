start_time = time()

ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])
subset_idx = parse(Int64, ARGS[6])
run_number = parse(Int64, ARGS[7])


snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq/"
snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-run$(run_number)"

if isfile("$(snaq_prefix).runtime")
    @info "$(snaq_prefix).runtime already exists, quitting!"
    exit(0)
end


@info "Loading packages and helpers"
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


pruned_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/pruned_gts$(subset_idx).tre"

tre0_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/ASTRAL$(subset_idx).tre"


# Find the true number of reticulations
@info "Loading true net"
truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
@info "Loading est gts"
est_gts = readMultiTopology(pruned_gt_file)

leaf_names = [leaf.name for leaf in est_gts[1].leaf]
n_retic = simple_prune(truenet, leaf_names).numHybrids

@info "Counting quartets"
q, t = countquartetsintrees(est_gts)

@info "Converting quartets to table"
df = readTableCF(writeTableCF(q, t))

@info "Loading astral tree"
tre0 = readTopology(tre0_file)
if run_number > 3
    tre0 = est_gts[run_number]
end


@info "Starting SNaQ"
seed = set_seed!(ntaxa, rep, ils, ngt, m, subset_idx, run_number)
snaq!(tre0, df, hmax = n_retic, seed = seed, filename = snaq_prefix, runs = 1)

runtime_seconds = time() - start_time
open("$(snaq_prefix).runtime", "w+") do f
    write(f, "$(runtime_seconds)")
end