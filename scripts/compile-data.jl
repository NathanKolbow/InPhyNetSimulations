const FORCE_RESAMPLE_ALL = false



using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate();
Pkg.update();

ENV["JULIA_DEPOT_PATH"] = joinpath(pwd(), "..")


using PhyloNetworks, DataFrames, CSV, InPhyNet, Base.Threads
include(joinpath(@__DIR__, "subscripts", "gtee.jl"))

global df = DataFrame(
    ntaxa=Int[], ngt=Int[], ils=String[], nbp=Int[], m=Int[], r=Int[], imethod=String[],
    gtee=Float64[], hwcd=Int[], input_error=Float64[],
    runtime_parallel=Float64[], runtime_serial=Float64[]
)
output_path = joinpath(@__DIR__, "..", "data", "all.csv")


if !FORCE_RESAMPLE_ALL && isfile(output_path)
    @warn "NOT resampling data."
    df = CSV.read(output_path, DataFrame)
    @info "$(nrow(df)) lines in df already."
else
    @warn "Resampling all data - this will take a while."
end


rows = Array{Any}(undef, 4*2*2*2*2*10*2)
is_valid = falses(4*2*2*2*2*10*2)

global j = Atomic{Int}(0)
for ntaxa in [30, 50, 100, 200]
for ngt in [100, 1000]
for ils in ["low", "high"]
for nbp in [100, 1000]
for m in [10, 20]
for r = 1:10
for imethod in ["snaq", "squirrel"]
    global j
    atomic_add!(j, 1)
    if Threads.threadid() == 1
        print("\r$(j[]) / $(length(rows)) (total $(sum(is_valid) + nrow(df)))")
    end

    if nrow(filter(row -> row.ntaxa == ntaxa && row.ngt == ngt && row.ils == ils && row.nbp == nbp && row.m == m && row.r == r && row.imethod == imethod, df)) > 0 continue end

    basedir = joinpath(@__DIR__, "..", "data", string(ntaxa), string(ngt), ils, string(nbp), string(m), string(r))
    if !isdir(basedir) continue end
    if !isfile(joinpath(basedir, "inphynet-$(imethod).net")) || !isfile(joinpath(basedir, "inphynet-$(imethod).runtime"))
        continue
    end
    is_valid[j[]] = true

    tnet = readnewick(joinpath(basedir, "true.net"))
    enet = readnewick(joinpath(basedir, "inphynet-$(imethod).net"))
    enet_runtime = parse(Float64, readlines(joinpath(basedir, "inphynet-$(imethod).runtime"))[1])

    constraint_runtimes = [parse(Float64, line) for line in readlines(joinpath(basedir, "$(imethod).runtime"))]

    egts = readmultinewick(joinpath(basedir, "estgts.tre"))
    tgts = readmultinewick(joinpath(basedir, "truegts.tre"))
    gtee_val = gtee(egts, tgts)

    # Input error calculation
    constraints = readmultinewick(joinpath(basedir, "$(imethod).net"))
    input_cerror = sum(hardwiredclusterdistance(c, prune_network(tnet, tiplabels(c)), false) for c in constraints)
    D, namelist = calculateAGID(egts)
    nj_tre = inphynet(D, namelist)
    D_error = Inf
    for dtre in displayedtrees(tnet, 0.0)
        D_error = min(D_error, hardwiredclusterdistance(nj_tre, dtre, false))
        if D_error == 0 break end
    end

    rows[j[]] = [
        ntaxa, ngt, ils, nbp, m, r, imethod,
        gtee_val, hardwiredclusterdistance(tnet, enet, false), input_cerror + D_error,
        enet_runtime + maximum(constraint_runtimes), enet_runtime + sum(constraint_runtimes)
    ]
end
end
end
end
end
end
end

println("\n$(sum(is_valid)) / $(length(is_valid))")
for idx in findall(is_valid)
    push!(df, rows[idx])
end
CSV.write(output_path, df)
