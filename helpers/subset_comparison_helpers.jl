using PhyloNetworks, InPhyNet, DataStructures, RCall

truegt_folder = joinpath(@__DIR__, "..", "est-gts/data/true-gts")
"""
To accomodate different styles of subset decomp, `subset_fxn` must
take the following arguments, in order:
    - truenet
    - true_gts
    - max_size
"""
function get_hwcd(subset_fxn::Function, ntaxa, rep, ils, ngt, m; verbose=true)
    p(msg::AbstractString) = println(msg)
    if !verbose p(msg::AbstractString) = msg end
    suffix = "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)";

    p("Loading true net...");
    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils);
    if (m == 10 || m == 20) && isfile("$(truegt_folder)/$(suffix).treefile_1")
        p("Loading true gts...");
        true_gts = [readTopology("$(truegt_folder)/$(suffix).treefile_$(j)") for j = 1:ngt]
    else
        p("Simulating true gts...");
        Random.seed!(m);
        true_gts = simulatecoalescent(truenet, ngt, 1)
    end

    p("Conducting subset decomposition...");
    subsets = subset_fxn(truenet, true_gts, m)

    p("Pruning true network into subsets...");
    cs = simple_prune(truenet, subsets)

    p("Calculating D...");
    D, namelist = calculateAGID(true_gts)

    p("Inferring full net...");
    mnet = nothing
    try
        mnet = inphynet(D, cs, namelist)
    catch
        @warn "Using pruneTruthFromDecomp"
        cs = pruneTruthFromDecomp(truenet, subsets)
        mnet = inphynet(D, cs, namelist)
    end

    p("Computing HWCD...");
    rootatnode!(mnet, "OUTGROUP")
    rootatnode!(truenet, "OUTGROUP")

    mnet_major = majorTree(mnet)
    true_major = majorTree(truenet)
    rootatnode!(mnet_major, "OUTGROUP")
    rootatnode!(true_major, "OUTGROUP")

    return hardwiredClusterDistance(mnet, truenet, true), 
        hardwiredClusterDistance(mnet_major, true_major, true),
        truenet.numHybrids,
        mnet.numHybrids
end


function wrapper_sateIdecomp(truenet::HybridNetwork, gts::AbstractVector{HybridNetwork}, m::Real)
    D, namelist = calculateAGID(gts)
    nj_tre = nj(DataFrame(D, namelist))
    return sateIdecomp(nj_tre, m)
end


function cheating_SC_decomp(truenet::HybridNetwork, gts::AbstractVector{HybridNetwork}, m::Real)

    Random.seed!(42)
    # 1. take 1 taxa from H, each sister clade, and each sister clade's parent
    nH = truenet.numHybrids;
    SC_subsets = [];
    for H in truenet.hybrid
        H_cycle = gather_cycle(truenet, H)
        s1 = H_cycle[2]
        s2 = H_cycle[length(H_cycle)-1]
        s3, s4 = (nothing, nothing)
        if length(H_cycle) >= 5
            s3 = H_cycle[3]
        end
        if length(H_cycle) >= 6
            s4 = H_cycle[length(H_cycle)-2]
        end

        all_descs = [];
        subset_options = [];
        for node in [H, s1, s2, s3, s4]
            if node === nothing continue end
            node_descs = get_descendant_leaves(node)
            if node == H
                all_descs = node_descs
                push!(subset_options, all_descs)
                continue
            end

            if all(child in [H, s1, s2, s3, s4] for child in getchildren(node))
                # All children of node are in the cycle, so we need to go UP, not down
                node_descs = setdiff(truenet.leaf, node_descs)
            end
            node_descs = setdiff(union(node_descs, all_descs), all_descs)
            if length(node_descs) > 0
                push!(subset_options, node_descs)
            end
            all_descs = union(node_descs, all_descs)
        end
        # length(unique([H, s1, s2, s3, s4])) == 5 || error("Unique entries has length $(length(unique([H, s1, s2, s3, s4])))")

        # For each set of descendants, pick the one from each set that is closest to H
        H_subset = [];
        for subset in subset_options
            node_dists = [get_node_dist(subset_node, H) for subset_node in subset]
            push!(H_subset, subset[findmin(node_dists)[2]])
        end

        push!(SC_subsets, [n.name for n in H_subset])
    end
    
    # 2. if any subsets have overlapping taxa, collapse them
    #    this could obviously be made more robust by not using
    #    random sampling, but we ignore that here
    # @info "n subsets before collapse: $(length(SC_subsets))"
    @info "B"
    none_overlapping = false;
    j = 0;
    while !none_overlapping
        j += 1;
        if j >= 1000 error("Reached 1000 iterations of while loop, quitting...") end
        none_overlapping = true
        for i = 1:(length(SC_subsets)-1)
            for j = (i+1):length(SC_subsets)
                sub1 = SC_subsets[i]; sub2 = SC_subsets[j];
                if length(intersect(sub1, sub2)) != 0
                    none_overlapping = false
                    SC_subsets[i] = union(SC_subsets[i], SC_subsets[j])
                    deleteat!(SC_subsets, j)
                    break
                end
                if !none_overlapping break end
            end
            if !none_overlapping break end
        end
    end
    # @info all(length(s) == length(unique(s)) for s in SC_subsets)
    # @info "n subsets after collapse: $(length(SC_subsets))"
    # @info "sum(length(ss)) = $(sum([length(s) for s in SC_subsets]))"
    # @info "length(unique(ss)) = $(length(unique(reduce(vcat, SC_subsets))))"
    # @show SC_subsets

    # 3. prune the taxa in these subsets off of a copied truenet
    #    so that we can guarantee these subsets appear together
    #    after sateI subset selection
    trimmed_taxa_set = tipLabels(truenet)
    for S in SC_subsets
        # Delete all but the first taxon in the set
        for taxon in S[2:length(S)]
            deleteat!(trimmed_taxa_set, findfirst(t -> t == taxon, trimmed_taxa_set))
        end
    end

    D, namelist = calculateAGID(gts)
    idx_filter = [j for j = 1:length(namelist) if namelist[j] in trimmed_taxa_set]
    D = D[idx_filter, idx_filter]
    namelist = namelist[idx_filter]
    tre0 = nj(DataFrame(D, namelist))
    sateI_subsets = sateIdecomp(tre0, m)

    for S in SC_subsets
        matching_sateI = sateI_subsets[findfirst(sat -> S[1] in sat, sateI_subsets)]
        for taxon in S[2:length(S)]
            push!(matching_sateI, taxon)
        end
    end
    # @info "Max subset size: $(maximum(length(s) for s in sateI_subsets))"
    return sateI_subsets

end


function cheating_SC_decomp(truenet::HybridNetwork, m::Real)

    Random.seed!(42)

    # 1. take 1 taxa from H, each sister clade, and each sister clade's parent
    nH = truenet.numHybrids;
    SC_subsets = [];
    for H in truenet.hybrid
        H_cycle = gather_cycle(truenet, H)
        s1 = H_cycle[2]
        s2 = H_cycle[length(H_cycle)-1]
        s3, s4 = (nothing, nothing)
        if length(H_cycle) >= 5
            s3 = H_cycle[3]
        end
        if length(H_cycle) >= 6
            s4 = H_cycle[length(H_cycle)-2]
        end

        all_descs = [];
        subset_options = [];
        for node in [H, s1, s2, s3, s4]
            if node === nothing continue end
            node_descs = get_descendant_leaves(node)
            if node == H
                all_descs = node_descs
                push!(subset_options, all_descs)
                continue
            end

            if all(child in [H, s1, s2, s3, s4] for child in getchildren(node))
                # All children of node are in the cycle, so we need to go UP, not down
                node_descs = setdiff(truenet.leaf, node_descs)
            end
            node_descs = setdiff(union(node_descs, all_descs), all_descs)
            if length(node_descs) > 0
                push!(subset_options, node_descs)
            end
            all_descs = union(node_descs, all_descs)
        end
        # length(unique([H, s1, s2, s3, s4])) == 5 || error("Unique entries has length $(length(unique([H, s1, s2, s3, s4])))")

        # For each set of descendants, pick the one from each set that is closest to H
        H_subset = [];
        for subset in subset_options
            node_dists = [get_node_dist(subset_node, H) for subset_node in subset]
            push!(H_subset, subset[findmin(node_dists)[2]])
        end

        push!(SC_subsets, [n.name for n in H_subset])
    end
    
    # 2. if any subsets have overlapping taxa, collapse them
    #    this could obviously be made more robust by not using
    #    random sampling, but we ignore that here
    # @info "n subsets before collapse: $(length(SC_subsets))"
    none_overlapping = false;
    j = 0;
    while !none_overlapping
        j += 1;
        if j >= 1000 error("Reached 1000 iterations of while loop, quitting...") end
        none_overlapping = true
        for i = 1:(length(SC_subsets)-1)
            for j = (i+1):length(SC_subsets)
                sub1 = SC_subsets[i]; sub2 = SC_subsets[j];
                if length(intersect(sub1, sub2)) != 0
                    none_overlapping = false
                    SC_subsets[i] = union(SC_subsets[i], SC_subsets[j])
                    deleteat!(SC_subsets, j)
                    break
                end
                if !none_overlapping break end
            end
            if !none_overlapping break end
        end
    end
    # @info all(length(s) == length(unique(s)) for s in SC_subsets)
    # @info "n subsets after collapse: $(length(SC_subsets))"
    # @info "sum(length(ss)) = $(sum([length(s) for s in SC_subsets]))"
    # @info "length(unique(ss)) = $(length(unique(reduce(vcat, SC_subsets))))"
    # @show SC_subsets

    # 3. prune the taxa in these subsets off of a copied truenet
    #    so that we can guarantee these subsets appear together
    #    after sateI subset selection
    trimmed_taxa_set = tipLabels(truenet)
    for S in SC_subsets
        # Delete all but the first taxon in the set
        for taxon in S[2:length(S)]
            deleteat!(trimmed_taxa_set, findfirst(t -> t == taxon, trimmed_taxa_set))
        end
    end
    
    D, namelist = internodedistance(majorTree(truenet))
    idx_filter = [j for j = 1:length(namelist) if namelist[j] in trimmed_taxa_set]
    D = Matrix{Float64}(D[idx_filter, idx_filter])
    namelist = namelist[idx_filter]
    # tre0 = PhyloNetworks.nj!(D, namelist) UNBELIEVABLY SLOW!!!!
    tre0 = inphynet(D, Vector{HybridNetwork}([]), namelist)
    sateI_subsets = sateIdecomp(tre0, m)

    for S in SC_subsets
        matching_sateI = sateI_subsets[findfirst(sat -> S[1] in sat, sateI_subsets)]
        for taxon in S[2:length(S)]
            push!(matching_sateI, taxon)
        end
    end
    # @info "Max subset size: $(maximum(length(s) for s in sateI_subsets))"
    return sateI_subsets

end


function simple_prune(truenet::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    nets = Array{HybridNetwork}(undef, length(subsets))
    for (i, set) in enumerate(subsets)
        tempnet = deepcopy(truenet)
        namelist = [leaf.name for leaf in tempnet.leaf]
        for name in namelist
            if !(name in set)
                PhyloNetworks.deleteleaf!(tempnet, name)
            end
        end
        nets[i] = readTopology(writeTopology(tempnet))
    end
    return nets
end


const Node = PhyloNetworks.Node;
"""
Finds and returns the nodes and edges that form the cycle associated with `H`.
"""
function gather_cycle(net::HybridNetwork, H::Node)
    H.hybrid || error("H must be a hybrid node")

    # 2-cycles
    if getparent(H) == getparentminor(H) return [H, getparent(H)] end

    Q = Queue{Vector{Node}}(); enqueue!(Q, [H, getparent(H)]);
    valid_paths = Vector{Vector{Node}}();
    
    while !isempty(Q)
        curr_path = dequeue!(Q)
        path_tip = curr_path[length(curr_path)]
        for adj_node in InPhyNet.getnodes(path_tip)
            if length(curr_path) > 2 && adj_node == H
                return [curr_path; adj_node]
            end
            if adj_node in curr_path continue end
            enqueue!(Q, [curr_path; adj_node])
        end
    end
end


"""
Gets the distance between two nodes
"""
function get_node_dist(node1::Node, node2::Node)
    Q = Queue{Vector{Node}}(); enqueue!(Q, [node1]);
    valid_path = nothing;

    while !isempty(Q) && valid_path === nothing
        curr_path = dequeue!(Q)
        tip_node = curr_path[length(curr_path)]

        for neighbor in InPhyNet.getnodes(tip_node)
            if neighbor == node2
                valid_path = [curr_path; node2]
                break
            elseif neighbor in curr_path
                continue
            else
                enqueue!(Q, [curr_path; neighbor])
            end
        end
    end
    if valid_path === nothing error("????") end

    path_length = 0.0
    for i = 1:(length(valid_path)-1)
        from = valid_path[i]
        to = valid_path[i+1]
        connecting_edge = from.edge[findfirst(e -> from in e.node && to in e.node, from.edge)]
        path_length += connecting_edge.length
    end
    return path_length
end




empirical_SC_k10(t, gt, m) = empirical_SC_decomp(t, gt, m, k=10)
empirical_SC_k15(t, gt, m) = empirical_SC_decomp(t, gt, m, k=15)
empirical_SC_k25(t, gt, m) = empirical_SC_decomp(t, gt, m, k=25)
empirical_SC_k50(t, gt, m) = empirical_SC_decomp(t, gt, m, k=50)
function empirical_SC_decomp(tnet::HybridNetwork, gts::Vector{HybridNetwork}, m::Int; k::Int=25)

    # 1. NJ tree
    @info "Inferring NJ tree"
    D, namelist = calculateAGID(gts)
    nj_tre = nj(DataFrame(D, namelist))
    rootatnode!(nj_tre, "OUTGROUP")
    cladewiseorder!(nj_tre)
    ordered_leaves = [node.name for node in nj_tre.node[nj_tre.cladewiseorder_nodeIndex] if node.leaf && node.name != "OUTGROUP"]

    # 2. Subset decomposition
    @info "Performing broad ToB subset decomp"
    Random.seed!(42)
    # k = 25
    subsets_clade = [];
    n_subsets = length(ordered_leaves) ÷ k
    for j = 1:n_subsets
        push!(subsets_clade, ordered_leaves[((j-1)*k+1):(j*k)])
        # Shuffle here so that `subsets_ancestral` is randomly sampled
        subsets_clade[j] = sample(subsets_clade[j], k, replace=false)
        push!(subsets_clade[j], "OUTGROUP")
    end

    subsets_ancestral = [];
    for j = 1:n_subsets
        push!(subsets_ancestral, [])
        for i = 1:k
            push!(subsets_ancestral[j], ordered_leaves[(i-1)*n_subsets+j])
        end
        push!(subsets_ancestral[j], "OUTGROUP")
    end


    # 3. Write the pruned gts for these subsets to files
    @info "Pruning and writing clade subsets"
    Threads.@threads for j=1:length(subsets_clade)
        set = subsets_clade[j]
        writeMultiTopology([pruneTruthFromDecomp(gt, Vector{String}(set)) for gt in gts], "temp/clade_subset$(j).tre")
    end
    @info "Pruning and writing ancestral subsets"
    Threads.@threads for j=1:length(subsets_ancestral)
        set = subsets_ancestral[j]
        writeMultiTopology([pruneTruthFromDecomp(gt, Vector{String}(set)) for gt in gts], "temp/ancestral_subset$(j).tre")
    end


    # 4. Infer the trees of blobs
    @info "Inferring trees of blobs 1/2"
    tobs = [];
    blob_centers = [];
    R"library(ape); library(MSCquartets);"
    for j = 1:n_subsets
        @info "1/2: $(j)/$(n_subsets)"
        fname = "temp/clade_subset$(j).tre"
        tobname = "temp/clade_tob$(j).tre"
        @rput fname
        @rput tobname
        if isfile(tobname) rm(tobname) end
        R"tt <- TINNIK(fname, alpha=0.01, plot=FALSE)"; # MAYBE TRY alpha=0.05??
        R"write.tree(tt$ToB, file=tobname)";
        tob = readTopology(tobname);
        push!(tobs, tob)
        push!(blob_centers, [node for node in tob.node if length(node.edge) > 3])
    end
    @info "Inferring trees of blobs 2/2"
    for j = 1:n_subsets
        @info "2/2: $(j)/$(n_subsets)"
        fname = "temp/ancestral_subset$(j).tre"
        tobname = "temp/ancestral_tob$(j).tre"
        @rput fname
        @rput tobname
        if isfile(tobname) rm(tobname) end
        R"tt <- TINNIK(fname, alpha=0.01, plot=FALSE)"; # MAYBE TRY alpha=0.05??
        R"write.tree(tt$ToB, file=tobname)";
        tob = readTopology(tobname);
        push!(tobs, tob)
        push!(blob_centers, [node for node in tob.node if length(node.edge) > 3])
    end

    temp_tobs = tobs
    temp_centers = blob_centers
    tobs = [];
    blob_centers = [];
    for (tob, cents) in zip(temp_tobs, temp_centers)
        if length(cents) == 0 continue end
        for c in cents
            push!(blob_centers, c)
            push!(tobs, tob)
        end
    end
    @info "Found $(length(blob_centers)) blobs ($(tnet.numHybrids) hybrids in true net)"


    # 5. Gather cycle subsets
    @info "Gathering cycle subsets"
    SC_subsets = [];
    for (tob, center_node) in zip(tobs, blob_centers)
        cycle_nodes = [ifelse(e.node[1] == center_node, e.node[2], e.node[1]) for e in center_node.edge]
        cycle_node_descs = [get_descendant_leaves(node) for node in cycle_nodes];
        idxs = sortperm([length(s) for s in cycle_node_descs])

        all_descs = [];
        subset_options = [];
        for (node, node_descs) in zip(cycle_nodes[idxs], cycle_node_descs[idxs])
            if node.leaf push!(all_descs, node); push!(subset_options, [node]); continue end

            if all(child in cycle_nodes for child in getchildren(node))
                node_descs = setdiff(tob.leaf, node_descs)
            end
            node_descs = setdiff(union(node_descs, all_descs), all_descs)
            if length(node_descs) > 0
                push!(subset_options, node_descs)
            end
            all_descs = union(node_descs, all_descs)
        end

        H_subset = [];
        for subset in subset_options
            node_dists = [get_node_dist(subset_node, center_node) for subset_node in subset]
            push!(H_subset, subset[findmin(node_dists)[2]])
        end
        push!(SC_subsets, [n.name for n in H_subset])
    end

    #### remove "OUTGROUP" from each subset
    for s in SC_subsets
        out_idx = findfirst(tax -> tax == "OUTGROUP", s)
        if out_idx !== nothing deleteat!(s, out_idx) end
    end


    # 6. collapse overlapping taxa
    @info "Collapsing overlapping taxa"
    none_overlapping = false;
    j = 0;
    while !none_overlapping
        j += 1;
        if j >= 1000 error("Reached 1000 iterations of while loop, quitting...") end
        none_overlapping = true
        for i = 1:(length(SC_subsets)-1)
            for j = (i+1):length(SC_subsets)
                sub1 = SC_subsets[i]; sub2 = SC_subsets[j];
                if length(intersect(sub1, sub2)) != 0
                    none_overlapping = false
                    SC_subsets[i] = union(SC_subsets[i], SC_subsets[j])
                    deleteat!(SC_subsets, j)
                    break
                end
                if !none_overlapping break end
            end
            if !none_overlapping break end
        end
    end


    # 7. prune these taxa off of `nj_tre` so that we can guarantee
    #    these taxa appear together after subset decomp
    @info "Pruning NJ tree and performing sateI decomp"
    trimmed_taxa_set = tipLabels(nj_tre)
    for S in SC_subsets
        for taxon in S[2:length(S)]
            deleteat!(trimmed_taxa_set, findfirst(t -> t == taxon, trimmed_taxa_set))
        end
    end

    nj_copy = pruneTruthFromDecomp(nj_tre, trimmed_taxa_set)
    sateI_subsets = sateIdecomp(nj_copy, m)

    @info "Imputing required subsets onto sateI decomp"
    for S in SC_subsets
        matching_sateI = sateI_subsets[findfirst(sat -> S[1] in sat, sateI_subsets)]
        for taxon in S[2:length(S)]
            push!(matching_sateI, taxon)
        end
    end

    return sateI_subsets   # FINAL SUBSET CHOICES

end



function smart_prune(net::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    if length(subsets) == 1
        return [pruneTruthFromDecomp(net, subsets[1])]
    elseif length(subsets) <= 4
        return pruneTruthFromDecomp(net, subsets)
    else
        i = length(subsets) ÷ 2
        split1 = subsets[1:i]
        split2 = subsets[(i+1):length(subsets)]
        n1 = pruneTruthFromDecomp(net, reduce(vcat, split1))
        n2 = pruneTruthFromDecomp(net, reduce(vcat, split2))
        return [smart_prune(n1, split1); smart_prune(n2, split2)]
    end
end
# Trying to beat: median=1.119s, min-max=1-1.4s