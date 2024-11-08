include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/subset_comparison_helpers.jl")


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
    inphynet_runtime = @elapsed mnet = inphynet(D, constraints, namelist)

    # Metrics
    @info "\trooted greedy min hwcd"
    rooted_min_greedy_hwcd = find_minimum_retic_subset_hwcd_greedy(truenet, mnet, verbose=true, swaponerror=true, rooted=true)
    @info "\tunrooted greedy min hwcd"
    unrooted_min_greedy_hwcd = (rooted_min_greedy_hwcd == 0) ?
        0 : ((ntaxa >= 1000) ? -1 : find_minimum_retic_subset_hwcd_greedy(truenet, mnet, verbose=true, swaponerror=true, rooted=false))
    
    @info "\tmin NJ hwcd"
    nj_tre = nj(DataFrame(D, namelist))
    min_unrooted_nj_hwcd = edge_μ_dist(majorTree(truenet), majorTree(mnet)) # hardwiredClusterDistance(majorTree(truenet), majorTree(mnet), true)
    if ntaxa < 1000
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
        unrooted_major_hwcd=((ntaxa >= 1000) ? -1 : hardwiredClusterDistance(majorTree(truenet), majorTree(mnet), false)),
        unrooted_hwcd=((ntaxa >= 1000) ? -1 : hardwiredClusterDistance(truenet, mnet, false)),
        min_unrooted_nj_hwcd=min_unrooted_nj_hwcd
    )
    CSV.write(DAT_FILE, results, append=true)

end


if abspath(PROGRAM_FILE) == @__FILE__
    for ntaxa in [50, 100, 200]
        for rep = 1:10
            for ils in ["low", "high"]
                for m in [10, 20, 30]
                    run_sim(ntaxa, rep, ils, m)
                end
            end
        end
    end
end
