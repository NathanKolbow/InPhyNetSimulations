constraintfile = ARGS[1]
estgtfile = ARGS[2]
output = ARGS[3]
rt_output = ARGS[4]

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
Pkg.update()

using PhyloNetworks, InPhyNet
estgts = readmultinewick(estgtfile)
constraints = readmultinewick(constraintfile)

rt = time()
D, namelist = calculateAGID(estgts);
fullnet = inphynet(D, constraints, namelist)
rt = time() - rt

writenewick(fullnet, output)
open(rt_output, "w+") do f
    write(f, string(rt))
end
