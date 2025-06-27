constraintfile = ARGS[1]
estgtfile = ARGS[2]
output = ARGS[3]

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhyloNetworks, InPhyNet
estgts = readmultinewick(estgtfile)
constraints = readmultinewick(constraintfile)

D, namelist = calculateAGID(estgts);
fullnet = inphynet(D, constraints, namelist)

writenewick(fullnet, output)