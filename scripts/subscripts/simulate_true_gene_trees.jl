using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using PhyloNetworks, PhyloCoalSimulations, StatsBase, Random

netfile = ARGS[1]
ngt = parse(Int, ARGS[2])
ils = ARGS[3]
seed = parse(Int, ARGS[4])
output = ARGS[5]


# Adjust edges so that the average length matches
# the desired level of ILS
# Low ILS:  2.0
# High ILS: 0.5
truenet = readnewick(netfile)
mean_length = mean(e.length for e in truenet.edge)
desired_mean = ils == "low" ? 1.5 : 0.5
for e in truenet.edge
    e.length = e.length * desired_mean / mean_length
end


# Simulate gene trees
Random.seed!(seed)
gts = simulatecoalescent(truenet, ngt, 1);
writemultinewick(gts, output)