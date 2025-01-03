include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/perfect-sims/mu-representation/mu-representation.jl")


out_df = DataFrame(
    ntaxa = Int[], rep = Int[], ils = String[], 
    ngt = Int[], m = Int[], avg_gtee = Float64[], 
    norm_output_error = Float64[],
    output_hwcd = Float64[]
)
approx_norm_df = CSV.read("approx_normalized_errors.csv", DataFrame)

for (j, row) in enumerate(eachrow(approx_norm_df))
    try
        print("\r$(j) / $(nrow(approx_norm_df))")
        ntaxa, rep, ils, ngt, m = (row["ntaxa"], row["rep"], String(row["ils"]), row["ngt"], row["m"])
        avg_gtee = calc_avg_gtee(get_true_gts(ntaxa, rep, ils, ngt, m), get_est_gts(ntaxa, rep, ils, ngt, m))
        push!(out_df, [
            ntaxa, rep, ils, ngt, m, avg_gtee, row["min_greedy_approx_norm_hwcd"], row["unrooted_min_greedy_hwcd"]
        ])
        CSV.write("gtee.csv", out_df)
    catch e
    end
end

