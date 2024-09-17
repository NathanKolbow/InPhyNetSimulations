# Condor runs are taking a LONG TIME
#
# This script is run on our server in a `tmux` session and just
# iterates, running every est-gt sim


@info "Loading packages and helpers"
using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

using PhyloNetworks
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


function run_snaq(ntaxa, rep, ils, ngt, m, subset_idx, run_number)

    start_time = time()
    print_prefix = "n$(ntaxa)#$(rep), $(ils) w/ $(ngt)gt, m=$(m) subset $(subset_idx), $(run_number)/10"

    snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq/"
    snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-run$(run_number)"

    if isfile("$(snaq_prefix).runtime")
        printstyled("[ALREADY COMPLETE] ", color = :cyan)
        println("$(print_prefix)")
        return
    end

    pruned_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)/pruned_gts.treefile"
    tre0_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)/tre0.treefile"

    if !isfile(pruned_gt_file) || !isfile(tre0_file)
        printstyled("[TREEFILES MISSING] ", color = :red)
        println("$(print_prefix)")
        return 
    end

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

    printstyled("[COMPLETE] ", color = :green)
    println("$(print_prefix)")

end


const N_RUNS = 10
const N_SUBSETS = Dict((500, 10) => 77, (500, 20) => 40)

@info "Entering loop"
for run_number in 1:N_RUNS
    for ntaxa in [500]
        for rep in 1:2
            for ils in ["low", "high"]
                for m in [10, 20]
                    for ngt in [100, 1000]
                        for subset_idx in 1:N_SUBSETS[(ntaxa, m)]
                            run_snaq(ntaxa, rep, ils, ngt, m, subset_idx, run_number)
                        end
                    end
                end
            end
        end
    end
end