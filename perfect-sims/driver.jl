include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/subset_comparison_helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/subset_comparison_helpers.jl")

# mu-representation fxns
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/perfect-sims/mu-representation/mu-representation.jl")


function run_sim(ntaxa::Int, rep::Int, ils::String, m::Int)

    DAT_FILE = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/perfect-sims/out.csv"
    @info "n$(ntaxa) r$(rep) $(ils) m$(m)"

    # Check if already exists
    df = CSV.read(DAT_FILE, DataFrame)
    if nrow(filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.m == m, df)) > 0
        @info "\talready completed"
        return
    end

    @info "\tloading data"
    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
    D, namelist = internodedistance(majorTree(truenet))

    @info "\tsubset decomp"
    subsets = cheating_SC_decomp(truenet, m)
    constraints = simple_prune(truenet, subsets)

    @info "\tinphynet"
    inphynet_runtime, mnet = -1, -1
    try
        inphynet_runtime = @elapsed mnet = inphynet(D, constraints, namelist)
    catch e
        @error "\tCaught error: $(typeof(e)); quitting and continuing"
        return -1
    end
    @info "\t\t|truenet.H| = $(truenet.numHybrids)"
    @info "\t\t|mnet.H| = $(mnet.numHybrids)"

    # Metrics
    unrooted_hwcd = -1
    if ntaxa < 1000
        @info "\tunrooted HWCD"
        unrooted_hwcd = hardwiredClusterDistance(truenet, mnet, false)
    end

    @info "\t$(ntaxa >= 1000 ? "" : "un")rooted greedy min hwcd"
    min_greedy_hwcd = (truenet.numHybrids == mnet.numHybrids ? unrooted_hwcd : find_minimum_retic_subset_hwcd_greedy(truenet, mnet, verbose=true, swaponerror=true, rooted=(ntaxa >= 1000))[1])
    
    @info "\tmin NJ hwcd"
    nj_tre = inphynet(D, Vector{HybridNetwork}([]), namelist)
    min_unrooted_nj_hwcd = edge_μ_dist(majorTree(truenet), majorTree(mnet)) # hardwiredClusterDistance(majorTree(truenet), majorTree(mnet), true)
    if ntaxa < 200
        for disp_tre in displayedTrees(truenet, 0.0)
            disp_hwcd = edge_μ_dist(disp_tre, nj_tre)
            if disp_hwcd < min_unrooted_nj_hwcd
                min_unrooted_nj_hwcd = disp_hwcd
                if min_unrooted_nj_hwcd == 0 break end
            end
        end
    end

    results = DataFrame(
        ntaxa=ntaxa, rep=rep, ils=ils, m=m,
        inphynet_runtime=inphynet_runtime,
        nretic_true=truenet.numHybrids,
        nretic_est=mnet.numHybrids,
        unrooted_major_hwcd=edge_μ_dist(majorTree(truenet), majorTree(mnet)) / 2,
        unrooted_hwcd=((ntaxa >= 1000) ? -1 : hardwiredClusterDistance(truenet, mnet, false)),
        min_greedy_hwcd=min_greedy_hwcd,
        is_greedy_rooted=(ntaxa >= 1000),
        min_unrooted_nj_hwcd=min_unrooted_nj_hwcd
    )
    CSV.write(DAT_FILE, results, append=true)

end


if abspath(PROGRAM_FILE) == @__FILE__
    for ntaxa in [50, 100, 200]
        for rep = 1:100
            for ils in ["low", "high"]
                for m in [10, 20, 30, 40]
                    run_sim(ntaxa, rep, ils, m)
                end
            end
        end
    end

    for ntaxa in [1000, 2500]
        for rep = 1:10
            for ils in ["low"]
                for m in [10, 20, 30, 40]
                    run_sim(ntaxa, rep, ils, m)
                end
            end
        end
    end
end
