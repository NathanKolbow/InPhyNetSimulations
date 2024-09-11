using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

using InPhyNet
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


const DAT_FILE = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/pruned-sims/dat.csv"

ntaxas = [parse(Int64, ARGS[1])]


# ntaxas = [500, 1000]
reps = collect(1:10)
ilss = ["low", "med", "high"]
min_ms = [5, 7]
max_ms = [10, 20, 30, 40, 50]
ngts = [10, 100, 1000, 5000]
repititionss = [1, 2, 5]
d_metrics = [calculateAGID, calculateAGIC]
nsim = 10


run_pruned_simulation_suite(
    ntaxas, reps, ilss, min_ms, max_ms, ngts, repititionss, d_metrics, nsim, DAT_FILE
)

