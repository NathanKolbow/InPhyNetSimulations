using Sockets
using Pkg

empty!(DEPOT_PATH)
host = gethostname()
if host == "submit1" host = "solislemus-001" end
depot_dir = joinpath(@__DIR__, "..", "depot_paths", split(gethostname(), ".")[1])

if !isdir(depot_dir)
    mkdir(depot_dir)
end
if !isfile(joinpath(depot_dir, "Project.toml"))
    cp(joinpath(@__DIR__, "..", "Project.toml"), joinpath(depot_dir, "Project.toml"))
end

push!(DEPOT_PATH, depot_dir)
ENV["JULIA_PKG_CACHE"] = depot_dir
# ENV["JULIA_DEPOT_PATH"] = depot_dir

using Pkg
Pkg.activate(depot_dir)
Pkg.instantiate()
Pkg.update()
Pkg.precompile()