include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


const MAX_CORES = 32

df = DataFrame(ntaxa=Int[], rep=Int[], ils=String[], ngt=Int[], m=Int[], ncores=Int[], runtime=Float64[])
all_params = get_completed_params()
all_rts = zeros(length(all_params), MAX_CORES)

Threads.@threads for j = 1:length(all_params)
    ntaxa, rep, ils, ngt, m = all_params[j]

    rts = get_runtimes(ntaxa, rep, ils, ngt, m)
    if rts === nothing continue end
    @info (ntaxa, rep, ils, ngt, m)
    
    for ncores = 1:MAX_CORES
        all_rts[j, ncores] = calculate_parallel_runtime(rts, ncores)
    end
end

for j = 1:length(all_params)
    ntaxa, rep, ils, ngt, m = all_params[j]
    for ncores = 1:MAX_CORES
        push!(df, [ntaxa, rep, ils, ngt, m, ncores, all_rts[j, ncores]])
    end
end

CSV.write("runtime.csv", df)
