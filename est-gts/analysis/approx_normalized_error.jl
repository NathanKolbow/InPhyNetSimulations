include(joinpath(@__DIR__, "..", "..", "helpers/helpers.jl"))
include(joinpath(@__DIR__, "..", "..", "perfect-sims/mu-representation/mu-representation.jl"))


complete_df = CSV.read(joinpath(@__DIR__, "approx_normalized_errors.csv"), DataFrame)

df = CSV.read(joinpath(@__DIR__, "..", "data", "out.csv"), DataFrame)
errors = zeros(nrow(df)) .- 1
min_errors = zeros(nrow(df)) .- 1
input_errors = zeros(nrow(df)) .- 1
nj_errors = zeros(nrow(df)) .- 1
for (j, row) in enumerate(eachrow(df))
    printstyled("\r$(j)/$(nrow(df))", color=:cyan)
    ntaxa, rep, ils, ngt, m = row[["ntaxa", "rep", "ils", "ngt", "m"]]
    truenet = load_true_net_ils_adjusted(ntaxa, rep, String(ils))

    complete_row = filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, complete_df)
    if nrow(complete_row) > 0
        errors[j] = complete_row[1, "approx_norm_hwcd"]
        min_errors[j] = complete_row[1, "min_greedy_approx_norm_hwcd"]
        input_errors[j] = complete_row[1, "approx_norm_input_hwcd"]
        nj_errors[j] = complete_row[1, "approx_norm_nj_hwcd"]
        continue
    end

    # Normalized output error
    errors[j] = row["unrooted_hwcd"] / (2 * (truenet.numEdges - truenet.numTaxa))
    min_errors[j] = row["unrooted_min_greedy_hwcd"] / (2 * (truenet.numEdges - truenet.numTaxa))

    # Normalized input error
    nj_errors[j] = row["min_unrooted_nj_hwcd"] / (2 * truenet.numTaxa - 6)
    
    input_error_num = row["min_unrooted_nj_hwcd"] + row["sum_constraint_hwcd"]
    input_error_div = (2 * truenet.numTaxa - 6)
    cs = get_constraints(ntaxa, rep, String(ils), ngt, m)
    if cs === nothing @info (ntaxa, rep, String(ils), ngt, m); continue; end
    for c in cs
        input_error_div += 2 * (c.numEdges - c.numTaxa)
    end
    input_errors[j] = input_error_num / input_error_div
end
println()
df[!,"approx_norm_hwcd"] = errors
df[!,"min_greedy_approx_norm_hwcd"] = min_errors
df[!,"approx_norm_input_hwcd"] = input_errors
df[!,"approx_norm_nj_hwcd"] = nj_errors
CSV.write("approx_normalized_errors.csv", df)
