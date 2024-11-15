include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/precompile-setup.jl")
using PhyloNetworks


iqtree_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/iqtree/"
consolidated_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/"

function consolidate(ntaxa, rep, ils, ngt, m)
    global iqtree_path
    global consolidated_path

    out_path = joinpath(consolidated_path, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
    if isfile(out_path) && length(readlines(out_path)) > 0
        printstyled("[ALREADY EXISTS] ", color = :green)
        println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")
        return
    end

    trees = Array{HybridNetwork}(undef, ngt)
    for j = 1:ngt
        j_treefile = joinpath(iqtree_path, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)_$(j).treefile")
        if !isfile(j_treefile) || length(readlines(j_treefile)) == 0
            printstyled("[IQTREE INCOMPLETE] ", color = :red)
            println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")
            return
        end

        # Remove non-tree iqtree files
        prefix = joinpath(iqtree_path, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)_$(j)")
        if isfile("$(prefix).bionj") rm("$(prefix).bionj") end
        if isfile("$(prefix).bionj") rm("$(prefix).ckp.gz") end
        if isfile("$(prefix).bionj") rm("$(prefix).iqtree") end
        if isfile("$(prefix).bionj") rm("$(prefix).log") end
        if isfile("$(prefix).bionj") rm("$(prefix).mldist") end
        if isfile("$(prefix).bionj") rm("$(prefix).model.gz") end

        # Save tree
        trees[j] = readTopology(j_treefile)
    end

    missing_idxs = findall(i -> !isassigned(trees, i), 1:ngt)
    if length(missing_idxs) != 0
        throw(ErrorException("Missing $(length(missing_idxs)) estimated gene trees!"))
    end

    writeMultiTopology(trees, out_path)

    printstyled("[COMPLETED] ", color = :yellow)
    println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")
end


if length(ARGS) == 5
    ntaxa = parse(Int64, ARGS[1])
    rep = parse(Int64, ARGS[2])
    ils = ARGS[3]
    ngt = parse(Int64, ARGS[4])
    m = parse(Int64, ARGS[5])

    consolidate(ntaxa, rep, ils, ngt, m)
elseif length(ARGS) != 0
    @error "Usage: julia 3b_consolidate-tree.jl [<ntaxa> <rep> <ils> <ngt> <m>]"
    exit()
else
    for ntaxa in [500, 1000]
        for rep in 1:1
            for ils in ["low", "high"]
                for ngt in [100, 1000]
                    for m in [10, 20]
                        consolidate(ntaxa, rep, ils, ngt, m)
                    end
                end
            end
        end
    end
end