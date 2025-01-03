include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


approx_df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/analysis/approx_normalized_errors.csv", DataFrame)

out_df = DataFrame(
    ntaxa = Int[], rep = Int[], ils = String[], 
    ngt = Int[], m = Int[], mnet_error = Float64[],
    #consensus_error = Float64[],
    filtered_super_error = Float64[],
    super_error = Float64[]
)
for rep = 1:10
    @info rep
    truenet = load_true_net_ils_adjusted(100, rep, "low")
    filt_sup_net, sup_net = get_splits_tree_networks(100, rep, "low", 100, 20)
    mnet = get_mnet(100, rep, "low", 100, 20)
    common_root!([truenet, filt_sup_net, sup_net, mnet], "OUTGROUP")

    push!(out_df, [
        100, rep, "low", 100, 20,
        hardwiredClusterDistance(truenet, mnet, true),
        #hardwiredClusterDistance(cons_net, mnet, true),
        hardwiredClusterDistance(filt_sup_net, mnet, true),
        hardwiredClusterDistance(sup_net, mnet, true)
    ])
end

CSV.write("compare_vs_other_methods.csv", out_df)




using Plots, StatsPlots
using CSV, DataFrames

df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/analysis/compare_vs_other_methods.csv", DataFrame)

boxplot(df[!,"mnet_error"], label = "mnet")
boxplot!(df[!,"filtered_super_error"], label = "filtered super")
boxplot!(df[!,"super_error"], label = "super")