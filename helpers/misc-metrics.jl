include("rf.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/perfect-sims/mu-representation/mu-representation.jl")



function calc_avg_gtee(ts_true::Vector{HybridNetwork}, ts_est::Vector{HybridNetwork})
    return mean(
        gtee(true_gt, est_gt) for (true_gt, est_gt) in zip(ts_true, ts_est)
    )
end

function gtee(t_true::HybridNetwork, t_est::HybridNetwork)
    return (edge_μ_dist(t_true, t_est) / 2) / (2 * t_true.numTaxa - 6)
end


"""
Calculates the `nretics_inside`, `nretics_outside`, and `nretics_duplicated` metrics.
"""
function calculateReticData(truenet::HybridNetwork, constraints::Vector{HybridNetwork})
    true_retic_names = [retic.name for retic in truenet.hybrid]
    constraint_retic_names = Set{String}()

    nretics_duplicated = 0
    for c in constraints
        # Log retics, counting duplicates
        c_reticnames = [retic.name for retic in c.hybrid]
        for retic in c_reticnames
            if retic in constraint_retic_names
                nretics_duplicated += 1
            else
                push!(constraint_retic_names, retic)
            end
        end
    end

    nretics_inside = length(constraint_retic_names)
    nretics_outside = truenet.numHybrids - nretics_inside

    return nretics_inside, nretics_outside, nretics_duplicated
end


"""
Calculates the network `net`'s pseudo-likelihood given gene trees `gts`.
Code taken from PhyloNetworks `pseudolik.jl`
"""
function calculate_net_logpseudolik(net::HybridNetwork, df::PhyloNetworks.DataCF)
    net0 = deepcopy(net)
    PhyloNetworks.parameters!(net0)
    Threads.@threads for q in df.quartet
        PhyloNetworks.extractQuartet!(net0, q)
        PhyloNetworks.calculateExpCFAll!(q.qnet)
    end

    return PhyloNetworks.logPseudoLik(df)
end


"""
Gets the HWCD error between `true_net` and `mnet` after removing reticulations in
`true_net` that do not appear in `mnet`. Hybrids must match in name for this
function to work properly.
"""
function get_error_without_missing_retics(true_net::HybridNetwork, mnet::HybridNetwork; try_root_at_node::String="OUTGROUP")
    true_net_copy = readTopology(writeTopology(true_net))
    true_hybs = true_net_copy.hybrid
    mnet_hybs = mnet.hybrid
    retics_to_remove = setdiff([h.name for h in true_hybs], [h.name for h in mnet_hybs])
    retics_to_remove = intersect(retics_to_remove, [h.name for h in true_hybs])
    
    for retic_name in retics_to_remove
        # Sometimes removing 1 hybrid will also result in another being removed,
        # so we need to make sure the findfirst result is valid
        retic_idx = findfirst([h.name for h in true_hybs] .== retic_name)
        if retic_idx === nothing continue end
    
        hybnode = true_hybs[retic_idx]
        PhyloNetworks.deletehybridedge!(true_net_copy, getparentedgeminor(hybnode))
    end

    try
        rootatnode!(true_net, "OUTGROUP")
    catch
    end
    try
        rootatnode!(mnet, "OUTGROUP")
    catch
    end
    
    return hardwiredClusterDistance(true_net_copy, mnet, true)
end


"""
Calculated hardwiredClusterDistance without multiplicity.
"""
function hwcd_no_multiplicity(truenet::HybridNetwork, mnet::HybridNetwork)
    try_outgroup_root(truenet)
    try_outgroup_root(mnet)

    true_clusters = [r[2:(truenet.numTaxa+1)] for r in eachrow(hardwiredClusters(truenet, tipLabels(truenet)))]
    m_clusters = [r[2:(truenet.numTaxa+1)] for r in eachrow(hardwiredClusters(mnet, tipLabels(truenet)))]

    return 2 * length(symdiff(
        true_clusters,
        m_clusters
    ))
end


"""
Counts the number of minor reticulate edges in `net` that have the same
hardwired clusters. Returns the total number of retics minus this number.
"""
function count_redundant_retics(net::HybridNetwork)
    retic_clusters = [hardwiredCluster(getparentedgeminor(H), tipLabels(net)) for H in net.hybrid]
    return net.numHybrids - length(unique(retic_clusters))
end


"""
Given that `truenet` and `mnet` have the same major tree and the same number
of reticulations, returns how many reticulations are misplaced in `mnet`.
"""
function n_retics_off(truenet::HybridNetwork, mnet::HybridNetwork)
    true_retic_edges = [hardwiredCluster(E, tipLabels(truenet)) for E in truenet.edge if getchild(E).hybrid]
    m_retic_edges = [hardwiredCluster(E, tipLabels(truenet)) for E in mnet.edge if getchild(E).hybrid]

    return length(symdiff(
        true_retic_edges,
        m_retic_edges
    ))
end


"""
Same as `find_minimum_retic_subset_hwcd` but does a greedy search instead of an exhaustive search.
"""
function find_minimum_retic_subset_hwcd_greedy(true_net::HybridNetwork, est_net::HybridNetwork; verbose::Bool=false, swaponerror::Bool=false)
    
    # Make sure the nets are properly rooted
    try_outgroup_root(true_net)
    try_outgroup_root(est_net)

    if est_net.numHybrids == true_net.numHybrids
        return hwcd_no_multiplicity(est_net, true_net), deepcopy(true_net)
    elseif est_net.numHybrids == 0
        return hwcd_no_multiplicity(est_net, majorTree(true_net)), deepcopy(true_net)
    elseif est_net.numHybrids > true_net.numHybrids
        if swaponerror
            temp = est_net
            est_net = true_net
            true_net = temp
        else
            throw(ArgumentError("est_net has > retics as true_net, maybe you put the arguments in the wrong order?"))
        end
    end

    true_copy = readTopology(writeTopology(true_net))
    while true_copy.numHybrids > est_net.numHybrids
        hyb_edge_numbers = [getparentedgeminor(h).number for h in true_copy.hybrid]
        hyb_edge_hwcd = zeros(length(hyb_edge_numbers))
        
        # Gather HWCDs after removing each retic
        for (j, hyb_number) in enumerate(hyb_edge_numbers)
            iter_copy = deepcopy(true_copy)
            PhyloNetworks.deletehybridedge!(iter_copy, iter_copy.edge[findfirst(e -> e.number == hyb_number, iter_copy.edge)])
            try_outgroup_root(iter_copy)
            hyb_edge_hwcd[j] = hwcd_no_multiplicity(iter_copy, est_net)
        end

        # Get rid of the retic w/ worst HWCD
        worst_edge = findmin(hyb_edge_hwcd)
        if verbose @info "HWCDs: $(hyb_edge_hwcd)" end
        if verbose @info "Removing $(worst_edge[2])" end
        PhyloNetworks.deletehybridedge!(true_copy, true_copy.edge[findfirst(e -> e.number == hyb_edge_numbers[worst_edge[2]], true_copy.edge)])
    end

    try_outgroup_root(true_copy)
    try_outgroup_root(est_net)
    return hwcd_no_multiplicity(est_net, true_copy), true_copy

end
# n1, n2 = readMultiTopology("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/networks/n500.netfile")[1:2]
# find_minimum_retic_subset_hwcd_greedy(n2, n1, verbose=true)

"""
`true_net` is a baseline network, `est_net` is a network estimated all the way up from
sequence data (typically). This fxn finds the subset of reticulations in `true_net` such
that `true_net` has only as many retics as `est_net` and the HWCD is minimized between
the two networks.
"""
function find_minimum_retic_subset_hwcd(true_net::HybridNetwork, est_net::HybridNetwork; verbose::Bool=false)
    # Make sure the nets are properly rooted
    try_outgroup_root(true_net)
    try_outgroup_root(est_net)

    # If they have the same number of retics, or `est_net` somehow has more retics
    # than `true_net`, return their HWCD
    if est_net.numHybrids == true_net.numHybrids
        return hardwiredClusterDistance(est_net, true_net, true)
    elseif est_net.numHybrids > true_net.numHybrids
        throw(ArgumentError("est_net has > retics as true_net, maybe you put the arguments in the wrong order?"))
    end

    # Find the subset!
    true_hyb_names = [h.name for h in true_net.hybrid]
    hyb_combinations = combinations(true_hyb_names, true_net.numHybrids - est_net.numHybrids)
    hwcds = Array{Float64}(undef, length(hyb_combinations))
    min_hwcd = Inf
    min_net = nothing
    n_combos = length(hyb_combinations)
    ac = AtomicCounter(0)

    Threads.@threads for (i, hyb_subset_names) in collect(enumerate(hyb_combinations))
        if i == 1 && verbose print("\rBeginning.") end
        @atomic :sequentially_consistent ac.iterspassed += 1
        # 1. Copy the true net
        true_net_copy = readTopology(writeTopology(true_net))
        
        # 2. Find the set of hybrids in `hyb_subset_names`
        remove_hybnodes = []
        for to_remove_name in hyb_subset_names
            push!(remove_hybnodes, true_net_copy.hybrid[findfirst([h.name == to_remove_name for h in true_net_copy.hybrid])])
        end

        # 3. Remove the set of hybrids from `true_net_copy`
        for hyb in remove_hybnodes
            PhyloNetworks.deletehybridedge!(true_net_copy, getparentedgeminor(hyb))
        end

        # 4. Compare w/ hardwiredClusterDistance, save if new minimum
        try
            rootatnode!(true_net_copy, "OUTGROUP")
        catch
        end
        hwcds[i] = hardwiredClusterDistance(true_net_copy, est_net, true)
        if hwcds[i] < min_hwcd
            min_hwcd = hwcds[i]
            min_net = true_net_copy
        end
        if verbose && Threads.threadid() == 1 print("\r$(ac.iterspassed)/$(n_combos): $(hwcds[i]), min = $(min_hwcd)             ") end
    end
    if verbose println("") end
    return minimum(hwcds), min_net
end


function try_outgroup_root(net::HybridNetwork)
    try
        rootatnode!(net, "OUTGROUP")
    catch
    end
end

function hardwiredClusterDistance_unrooted_mine(net1::HybridNetwork, net2::HybridNetwork)
    return hardwiredClusterDistance_unrooted_mine!(deepcopy(net1), deepcopy(net2))
end
function hardwiredClusterDistance_unrooted_mine!(net1::HybridNetwork, net2::HybridNetwork)
    #= fixit: inefficient function, because r1 * r2 "M" matrices of
      hardwiredClusters() are calculated, where ri = # root positions in neti.
      Rewrite to calculate only r1 + r2 M's.
    =#
    removedegree2nodes!(net1) # because re-rooting would remove them in an
    removedegree2nodes!(net2) # unpredictable order
    # find all permissible positions for the root
    net1roots = [n.number for n in net1.node if !n.leaf]
    #= disallow the root at a leaf: adding a degree-2 node adds a cluster
       that could be artificially matched to a cluster from a degree-3 node
       sister to a hybrid edge, when a the leaf edge is the donor. =#
    for i in length(net1roots):-1:1 # reverse order, to delete some of them
        try
            rootatnode!(net1, net1roots[i])
            # tricky: rootatnode adds a degree-2 node if i is a leaf,
            #         and delete former root node if it's of degree 2.
        catch e
            isa(e, PhyloNetworks.RootMismatch) || rethrow(e)
            deleteat!(net1roots, i)
        end
    end
    net2roots = [n.number for n in net2.node if !n.leaf]
    for i in length(net2roots):-1:1
        try
            rootatnode!(net2, net2roots[i])
        catch e
            isa(e, PhyloNetworks.RootMismatch) || rethrow(e)
            deleteat!(net2roots, i)
        end
    end

    bestdissimilarity = typemax(Int)
    bestns = missing
    n_combos = length(net1roots)*length(net2roots)
    ac = AtomicCounter(0)


    println("Iterating over $(n_combos) combinations")
    Threads.@threads for i = 1:length(net1roots)
        # Add this thread to the dict map if it isn't there already
        iter_net1 = deepcopy(net1)
        n1 = net1roots[i]
        rootatnode!(iter_net1, n1)
        vals = Array{Int64}(undef, length(net2roots))

        for j = 1:length(net2roots)
            @atomic :sequentially_consistent ac.iterspassed += 1
            n2 = net2roots[j]
            rootatnode!(net2, n2)

            idx = (i-1)*length(net2roots) + j
            try
                vals[j] = RobinsonFoulds(iter_net1, net2) # rooted = true now
            catch e
                # Something `hardwiredClusterDistance` errors out...
                rethrow(e)
                vals[j] = typemax(Int)
            end
            if vals[j] == 0 return (0, (n1, n2)) end
            if vals[j] < bestdissimilarity
                bestdissimilarity = vals[j]
            end

            print("\r$(ac.iterspassed)/$(n_combos): $(vals[j]), best = $(bestdissimilarity)                              ")
        end

        iter_best = minimum(vals)
        if iter_best <= bestdissimilarity
            bestdissimilarity = iter_best
            bestns = (net1roots[i], net2roots[findmin(vals)[2]])
        end
    end
    # @info "best root nodes: $bestns"
    # warning: original roots (and edge directions) NOT restored
    return bestdissimilarity, bestns
end








function collect_retry_data(netfile_name::String)
    data_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/simulation-study/simulation-scripts/data/"
    netid, replicatenum, ngt, seq_len, ils_level, maxsubsetsize, dmethod = split(netfile_name, "_")[2:8]
    netid = String(netid)
    replicatenum = parse(Int64, replicatenum)
    ngt = parse(Int64, ngt)
    seq_len = parse(Int64, seq_len)
    maxsubsetsize = parse(Int64, maxsubsetsize)
    ils_level = String(ils_level)

    est_constraints = readMultiTopology(joinpath(data_dir, netfile_name))
    est_gts = readMultiTopology(joinpath(data_dir, "estgt_$(netid)_$(replicatenum)_$(ngt)_$(seq_len)_$(ils_level).treefile"))
    est_D, est_namelist = calculateAGIC(est_gts)

    return est_constraints, est_D, est_namelist
end