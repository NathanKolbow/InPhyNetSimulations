using Random, DataFrames, CSV
combos = []


output_path = joinpath(@__DIR__, "..", "data", "all.csv")
df = CSV.read(output_path, DataFrame)


for n in [25, 50, 100, 200]
for ngt in [100, 1000]
for ils in ["low", "high"]
for nbp in [100, 1000]
for m in [10, 20]
for r in 1:10
for imethod in ["snaq", "phylonet", "phylonet-ml"]
    imethod == "phylonet-ml" && m != 10 && continue
    row = filter(ro -> ro.ntaxa == n .&& ro.ngt == ngt .&& ro.ils == ils .&& ro.nbp == nbp .&& ro.m == m .&& ro.r == r .&& ro.imethod == imethod, df)
    if nrow(row) > 0 continue end
    push!(combos, [n, ngt, ils, nbp, m, r, imethod])
end
end
end
end
end
end
end

# Write them in a random order so that any potential server
# slow downs are distributed rather than only being visible
# in some subset of the data
combos = combos[randperm(length(combos))]
open(joinpath(@__DIR__, "submit.table"), "w+") do f
    for (n, ngt, ils, nbp, m, r, imethod) in combos
        write(f, "$n,$ngt,$ils,$nbp,$m,$r,$imethod\n")
    end
end
