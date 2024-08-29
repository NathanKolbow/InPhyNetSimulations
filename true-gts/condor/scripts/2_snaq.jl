start_time = time()

@everywhere ntaxa = $(parse(Int64, ARGS[1]))
@everywhere rep = $(parse(Int64, ARGS[2]))
@everywhere ils = $(ARGS[3])
@everywhere ngt = $(parse(Int64, ARGS[4]))
@everywhere m = $(parse(Int64, ARGS[5]))
@everywhere subset_idx = $(parse(Int64, ARGS[6]))
@everywhere run_number = $(parse(Int64, ARGS[7]))


@info "Loading packages and helpers"
@everywhere begin
    using Pkg
    Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")
    
    using PhyloNetworks
    include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
end


astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/astral/"
astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

true_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/true-gts/"
true_gt_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/snaq/"
snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-r$(run_number)"


# Load true gene trees
truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
set_seed!(ntaxa, rep, ils, ngt, m, 0)
gts = simulatecoalescent(truenet, ngt, 1)


# Load astral tree and subset
tre0 = readTopology(astral_file)
subset = sateIdecomp(tre0, m)[subset_idx]


# Prune all trees to fit the subset, find n_retic
for j = 1:length(gts)
    gts[j] = pruneTruthFromDecomp(gts[j], subset)
end
tr0 = pruneTruthFromDecomp(tre0, subset)
n_retic = pruneTruthFromDecomp(truenet, subset).numHybrids


# Run SNaQ
q, t = countquartetsintrees(gts)
df = readTableCF(writeTableCF(q, t))

seed = set_seed!(ntaxa, rep, ils, ngt, m, subset_idx, run_number)
snaq!(tre0, df, hmax = n_retic, seed = seed, filename = snaq_prefix, runs = 1)

runtime_seconds = time() - start_time
open("$(snaq_prefix).runtime", "w+") do f
    write(f, "$(runtime_seconds)")
end