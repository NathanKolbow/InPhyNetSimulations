using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate();
Pkg.update();


using PhyloNetworks, DataFrames, CSV, InPhyNet
include(joinpath(@__DIR__, "subscripts", "gtee.jl"))

global df = DataFrame(
    ntaxa=Int[], ngt=Int[], ils=String[], nbp=Int[], m=Int[], r=Int[], imethod=String[],
    gtee=Float64[], hwcd=Int[], input_error=Float64[],
    runtime_serial=Float64[], runtime_parallel=Float64[]
)
output_path = joinpath(@__DIR__, "..", "data", "all.csv")

global j = 0
for ntaxa in [50, 100]
for ngt in [100, 1000]
for ils in ["low", "high"]
for nbp in [100, 1000]
for m in [10, 20]
for r = 1:10
for imethod in ["snaq", "squirrel"]
    global j
    j += 1
    print("\r$(j) / 640")

    basedir = joinpath(@__DIR__, "..", "data", string(ntaxa), string(ngt), ils, string(nbp), string(m), string(r))
    if !isdir(basedir) continue end
    if !isfile(joinpath(basedir, "inphynet-$(imethod).net")) || !isfile(joinpath(basedir, "inphynet-$(imethod).runtime"))
        continue
    end

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

    push!(df,
        [
            ntaxa, ngt, ils, nbp, m, r, imethod,
            gtee_val, hardwiredclusterdistance(tnet, enet, false), input_cerror + D_error,
            enet_runtime + maximum(constraint_runtimes), enet_runtime + sum(constraint_runtimes)
        ]
    )
    CSV.write(output_path, df)
end
end
end
end
end
end
end

CSV.write(output_path, df)