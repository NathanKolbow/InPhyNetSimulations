
function get_splits_tree_networks(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    net_file = joinpath(@__DIR__, "..", "/other-methods/splits-tree/output/")
    net_file = joinpath(net_file, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
    if !isfile(net_file) return nothing end
    # Consensus network has <1: ... so is broken at the moment
    return [readTopology(line) for (j, line) in enumerate(readlines(net_file)) if j > 1]
    #return readMultiTopology(net_file)
end


function get_est_gts(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    return readMultiTopology(
        joinpath(@__DIR__, "..", "est-gts", "data", "est-gts", "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
    )
end


function get_true_gts(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    return [readTopology(
        joinpath(@__DIR__, "..", "est-gts", "data", "true-gts", "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)", "$(i).treefile")
    ) for i = 1:ngt]
end


function get_astral_tree(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    return readTopology(
        joinpath(@__DIR__, "..", "est-gts", "data", "astral", "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
    )
end


function get_number_of_constraints(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    subset_file = joinpath(@__DIR__, "..", "/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab")
    if !isfile(subset_file) || length(readlines(subset_file)) == 0 return nothing end
    return Int64(length(readlines(subset_file)) / 10)
end


function get_subsets(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    n_subsets = get_number_of_constraints(ntaxa, rep, ils, ngt, m)
    if n_subsets === nothing return nothing end

    subsets = Vector{String}()
    for j = 1:n_subsets
        tre0_path = joinpath(@__DIR__, "..", "est-gts", "data", "subsets", "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)", "ASTRAL$(j).tre")
        if !isfile(tre0_path) return nothing end
        tre0 = readTopology(tre0_path)
        push!(subsets, tipLabels(tre0))
    end
    return subsets
end


function get_constraints(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    subset_file = joinpath(@__DIR__, "..",  "est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab")
    if !isfile(subset_file) || length(readlines(subset_file)) == 0 return nothing end
    n_subsets = Int64(length(readlines(subset_file)) / 10)

    snaq_newicks = Array{String}(undef, n_subsets, 10)
    runtimes = Array{Float64}(undef, n_subsets, 10)
    nlls = Array{Float64}(undef, n_subsets, 10)

    runtimes .= 0.0
    nlls .= Inf

    for subset_idx = 1:n_subsets
        for run_number = 1:10
            s = time()
            snaq_prefix = joinpath(@__DIR__, "..", "est-gts/data/snaq/")
            snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/subset$(subset_idx)/run$(run_number)"

            # isfile("$(snaq_prefix).runtime") || error("File $(snaq_prefix).runtime does not exist!")
            isfile("$(snaq_prefix).runtime") || return nothing

            out_info = split(readlines("$(snaq_prefix).out")[1], " -Ploglik = ")
            snaq_newicks[subset_idx, run_number] = out_info[1]
            nlls[subset_idx, run_number] = parse(Float64, out_info[2])
            runtimes[subset_idx, run_number] = parse(Float64, readlines("$(snaq_prefix).runtime")[1])
        end
    end

    constraints = Vector{HybridNetwork}([
        readTopology(newick) for newick in snaq_newicks[findmin(nlls, dims=2)[2]]
    ][:,1])
    for c in constraints
        while length(c.node[c.root].edge) != 2
            try
                rootatnode!(c, getchildren(c.node[c.root])[1])
            catch e
                rootatnode!(c, getchildren(c.node[c.root])[2])
            end
        end
    end
    return constraints
end


function get_pct_complete_constraints(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    subset_file = joinpath(@__DIR__, "..", "est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab")
    if !isfile(subset_file) || length(readlines(subset_file)) == 0 return 0.0 end
    n_subsets = Int64(length(readlines(subset_file)) / 10)

    n_complete = 0
    for subset_idx = 1:n_subsets
        for run_number = 1:10
            snaq_prefix = joinpath(@__DIR__, "..", "est-gts/data/snaq/")
            snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/subset$(subset_idx)/run$(run_number)"

            isfile("$(snaq_prefix).runtime") || continue
            n_complete += 1
        end
    end
    return n_complete / (n_subsets * 10)
end


function get_runtimes(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    subset_file = joinpath(@__DIR__, "..", "est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab")
    if !isfile(subset_file) || length(readlines(subset_file)) == 0 return nothing end
    n_subsets = Int64(length(readlines(subset_file)) / 10)

    runtimes = Array{Float64}(undef, n_subsets, 10)
    nlls = Array{Float64}(undef, n_subsets, 10)

    runtimes .= 0.0

    for subset_idx = 1:n_subsets
        for run_number = 1:10
            snaq_prefix = joinpath(@__DIR__, "..", "est-gts/data/snaq/")
            snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/subset$(subset_idx)/run$(run_number)"

            # isfile("$(snaq_prefix).runtime") || error("File $(snaq_prefix).runtime does not exist!")
            isfile("$(snaq_prefix).runtime") || return nothing
            
            runtimes[subset_idx, run_number] = parse(Float64, readlines("$(snaq_prefix).runtime")[1])
        end
    end

    return runtimes
end


function calculate_parallel_runtime(rts::Array{Float64}, ncores::Int)
    cores = [0.0 for _=1:ncores]    # value is time remaining in the current core
    time = 0.0
    for rt in rts
        _next_core = findmin(c -> c, cores)[2]
        elapsed = cores[_next_core]
        time += elapsed
        cores = [c - elapsed for c in cores]
        cores[_next_core] = rt
    end
    return time + maximum(cores)
end
calculate_parallel_runtime(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int, ncores::Int) =
    calculate_parallel_runtime(get_runtimes(ntaxa, rep, ils, ngt, m), ncores)


function get_mnet(ntaxa::Int, rep::Int, ils::String, ngt::Int, m::Int)
    
    est_gts = get_est_gts(ntaxa, rep, ils, ngt, m)
    cs = get_constraints(ntaxa, rep, ils, ngt, m)
    cs !== nothing || return nothing

    for c in cs
        while length(c.node[c.root].edge) != 2
            try
                rootatnode!(c, getchildren(c.node[c.root])[1])
            catch e
                rootatnode!(c, getchildren(c.node[c.root])[2])
            end
        end
    end

    D, namelist = calculateAGID(est_gts)
    return inphynet(D, cs, namelist)
end


function get_completed_params()
    df = CSV.read(joinpath(@__DIR__, "..", "est-gts", "data", "out.csv"), DataFrame)
    params = []
    for row in eachrow(df)
        push!(params, (row["ntaxa"], row["rep"], String(row["ils"]), row["ngt"], row["m"]))
    end
    return params
end