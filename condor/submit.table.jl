using Random
combos = []

for n in [25, 50, 100, 200]
for ngt in [100, 1000]
for ils in ["low", "high"]
for nbp in [100, 1000]
for m in [10, 20]
for r in 1:10
for imethod in ["snaq", "phylonet", "phylonet-ml"]
    imethod == "phylonet-ml" && m == 20 && continue
    push!(combos, [n, ngt, ils, nbp, m, r, imethod])
end
end
end
end
end
end
end

combos = combos[randperm(length(combos))]
open(joinpath(@__DIR__, "submit.table"), "w+") do f
    for (n, ngt, ils, nbp, m, r, imethod) in combos
        write(f, "$n,$ngt,$ils,$nbp,$m,$r,$imethod\n")
    end
end
