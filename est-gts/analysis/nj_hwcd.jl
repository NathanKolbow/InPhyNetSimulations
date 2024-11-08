using CSV, DataFrames, PhyloNetworks, StatsBase, InPhyNet
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
results_df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", DataFrame)


get_est_gts(ntaxa, rep, ils, ngt, m) = return readMultiTopology("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
get_true_gts(ntaxa, rep, ils, ngt, m) = return [readTopology("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile_$(j)") for j=1:ngt]

# ignore this one - the one below is better


args = [];
sum_constraint_hwcd = [];
unrooted_hwcd = [];
unrooted_min_greedy_hwcd = [];
for ntaxa in [50, 100]
    for rep = 1:5
        for ils in ["low", "high"]
            for ngt in [100, 1000]
                for m in [10, 20]

                    row_match = filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, results_df)
                    if nrow(row_match) > 0
                        @info "n$(ntaxa) #$(rep), $(ils), $(ngt)gt, m=$(m)"

                        row = row_match[1,:]
                        push!(args, (ntaxa, rep, ils, ngt, m))
                        push!(sum_constraint_hwcd, row["sum_constraint_hwcd"])
                        push!(unrooted_hwcd, row["unrooted_hwcd"])
                        push!(unrooted_min_greedy_hwcd, row["unrooted_min_greedy_hwcd"])

                    end

                end
            end
        end
    end
end



nj_hwcd = [];
for ntaxa in [50, 100]
    for rep = 1:5
        for ils in ["low", "high"]
            for ngt in [100, 1000]
                for m in [10, 20]

                    row_match = filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, results_df)
                    if nrow(row_match) > 0
                        @info "n$(ntaxa) #$(rep), $(ils), $(ngt)gt, m=$(m)"

                        truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
                        estgts = get_est_gts(ntaxa, rep, ils, ngt, m)
                        D_est, namelist_est = calculateAGID(estgts)
                        nj_tre = nj(DataFrame(D_est, namelist_est))

                        min_hwcd = Inf
                        for disp_tre in displayedTrees(truenet, 0.0)
                            disp_hwcd = hardwiredClusterDistance(disp_tre, nj_tre, false)
                            if disp_hwcd < min_hwcd
                                min_hwcd = disp_hwcd
                                if min_hwcd == 0 break end
                            end
                        end
                        push!(nj_hwcd, min_hwcd)

                    end

                end
            end
        end
    end
end

using Plots, StatsPlots

# dot plot
plot(sum_constraint_hwcd, unrooted_hwcd, label = "vs. sum constraint HWCD", seriestype=:scatter)
Plots.abline!(1, 0, color="red")
plot!(sum_constraint_hwcd .+ nj_hwcd, unrooted_hwcd, label = "vs. sum constraint HWCD", seriestype=:scatter)

# histogram of differences
hwcd_diffs = unrooted_hwcd .- (sum_constraint_hwcd .+ nj_hwcd)
histogram(hwcd_diffs)
qqnorm(hwcd_diffs)  # looks very close to normal!