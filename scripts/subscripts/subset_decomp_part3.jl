using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

# Parse ARGS
est_gt_file = ARGS[1]
subset_dir = ARGS[2]
astral_file = ARGS[3]
output = ARGS[4]
m = parse(Int, ARGS[5])

# ntaxa = parse(Int64, ARGS[1])
# rep = parse(Int64, ARGS[2])
# ils = ARGS[3]
# ngt = parse(Int64, ARGS[4])
# m = parse(Int64, ARGS[5])

# Load packages
@info "Loading packages"
using PhyloNetworks, InPhyNet, DataStructures, Random, StatsBase
Random.seed!(42)

# Helper function
function smart_prune(net::HybridNetwork, subsets::AbstractVector{<:AbstractVector{<:AbstractString}})
    if length(subsets) == 1
        return [prune_network(net, subsets[1])]
    elseif length(subsets) <= 4
        return prune_network(net, subsets)
    else
        i = length(subsets) ÷ 2
        split1 = subsets[1:i]
        split2 = subsets[(i+1):length(subsets)]
        n1 = prune_network(net, reduce(vcat, split1))
        n2 = prune_network(net, reduce(vcat, split2))
        return [smart_prune(n1, split1); smart_prune(n2, split2)]
    end
end

# Helper function
function get_descendant_leaves(node::PhyloNetworks.Node)
    q = Vector{PhyloNetworks.Node}([node])
    leaves = Vector{PhyloNetworks.Node}([])
    while length(q) > 0
        curr = pop!(q)
        if curr.leaf
            push!(leaves, curr)
        else
            for child in getchildren(curr)
                push!(q, child)
            end
        end
    end
    return leaves
end

# Helper function
const Node = PhyloNetworks.Node
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

# Read trees of blobs
@info "Reading trees of blobs"
tobs = [readTopology(joinpath(subset_dir, file)) for file in readdir(subset_dir) if length(file) > 4 && file[(length(file)-3):length(file)] == ".tob"]

# Gather blob centers
blob_centers = []
DEBUG_tobs = []
for tob in tobs
    for node in tob.node
        if length(node.edge) > 3 push!(blob_centers, node); push!(DEBUG_tobs, tob) end
    end
end
const N_BLOBS = length(blob_centers)

# Generate subsets from these blob centers
@info "Generating subsets from blob centers"
blob_subsets = [];
for (j, center_node) in enumerate(blob_centers)
    cycle_nodes = [ifelse(e.node[1] == center_node, e.node[2], e.node[1]) for e in center_node.edge]
    cycle_node_descs = [get_descendant_leaves(node) for node in cycle_nodes]
    idxs = sortperm([length(s) for s in cycle_node_descs])

    all_descs = [];
    subset_options = [];
    for (node, node_descs) in zip(cycle_nodes[idxs], cycle_node_descs[idxs])
        if node.leaf push!(all_descs, node); push!(subset_options, [node]); continue end

        if any(child == center_node for child in getchildren(node))
            continue
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
    push!(blob_subsets, [n.name for n in H_subset])
end

# Remove OUTGROUP from each subset (it shouldn't be any, but just in case)
for s in blob_subsets
    out_idx = findfirst(tax -> tax == "OUTGROUP", s)
    if out_idx !== nothing deleteat!(s, out_idx) end
end

# Combine any subsets with overlapping taxa
@info "Collapsing subsets with overlapping taxa"
none_overlapping = false;
j = 0;
n_skipped_due_to_m = 0;
while !none_overlapping
    global n_skipped_due_to_m
    global none_overlapping
    global blob_subsets
    global j
    j += 1
    if j >= 1000 error("Reached 1000 iterations of while loop, quitting...") end
    none_overlapping = true
    for i = 1:(length(blob_subsets)-1)
        for j = (i+1):length(blob_subsets)
            sub1 = blob_subsets[i]; sub2 = blob_subsets[j];
            if length(intersect(sub1, sub2)) != 0
                if length(union(sub1, sub2)) <= m
                    none_overlapping = false
                    blob_subsets[i] = union(blob_subsets[i], blob_subsets[j])
                    deleteat!(blob_subsets, j)
                    break
                else
                    # If we skip the collapse, we can't have duplicate entries so we have to remove those taxa from the largest subset
                    matching_taxa = intersect(sub1, sub2)
                    for taxa in matching_taxa
                        if length(sub1) > length(sub2)
                            deleteat!(sub1, findfirst(t -> t == taxa, sub1))
                        else
                            deleteat!(sub2, findfirst(t -> t == taxa, sub2))
                        end
                    end
                    none_overlapping = false
                    n_skipped_due_to_m += 1
                end
            end
            if !none_overlapping break end
        end
        if !none_overlapping break end
    end
end
if n_skipped_due_to_m > 0 printstyled("[ WARN ] ", color=:yellow); println("Skipped $(n_skipped_due_to_m) collapses due to m-value") end
printstyled("[ RESULT ]", color=:red); printstyled(" $(N_BLOBS) originally inferred blobs\n")
printstyled("[ RESULT ]", color=:red); printstyled(" $(length(blob_subsets)) restricted sets:\n", color=:black)
for (i, s) in enumerate(blob_subsets)
    printstyled("\tS_i, n=$(length(s)): ", color=:black); printstyled("$(s)\n", color=:red)
end

# Read ASTRAL tree
@info "Reading ASTRAL tree"
tre0 = readTopology(astral_file)
D0, namelist0 = internodedistance(tre0)

# Get trimmed taxa set that only includes that first taxa in each blob subset
trimmed_taxa_set = tipLabels(tre0)
for s in blob_subsets
    global trimmed_taxa_set
    for taxon in s[2:length(s)]
        deleteat!(trimmed_taxa_set, findfirst(t -> t == taxon, trimmed_taxa_set))
    end
end

# Compute SATe-I decomposition w/o blob taxa
@info "Computing SATe-I decomposition w/o blob taxa"
tre0_trimmed = prune_network(tre0, trimmed_taxa_set)
final_subsets = []
while final_subsets == []
    global min_size, final_subsets
    try
        final_subsets = sateIdecomp(tre0_trimmed, m)
    catch
        min_size -= 1
        if min_size < 0
            error("Could not find valid subset decomposition")
        end
    end
end

# Impute blob taxa back into SATe-I decomposition
@info "Imputing blob taxa into SATe-I decomposition"
for s in blob_subsets
    global final_subsets
    matching_subset = final_subsets[findfirst(sat -> s[1] in sat, final_subsets)]
    if length(matching_subset) + length(s) - 1 > m
        deleteat!(matching_subset, findfirst(tax -> tax == s[1], matching_subset))
        push!(final_subsets, s)
    else
        for taxon in s[2:length(s)]
            push!(matching_subset, taxon)
        end
    end
end
final_subsets = final_subsets[[length(s) > 0 for s in final_subsets]]

# If any subsets are larger than `m` at this point, split them in half
L = [length(s) for s in final_subsets];
while maximum(L) > m + 2    # give a little leeway, b/c we'd rather not break up cycles...
    global L
    _max_subset = findmax(L)[2]
    max_subset = final_subsets[_max_subset]
    _split_idx = length(max_subset) ÷ 2

    push!(final_subsets, max_subset[(_split_idx+1):length(max_subset)])
    final_subsets[_max_subset] = max_subset[1:_split_idx]
    L = [length(s) for s in final_subsets];
end

# If any subsets have <5 taxa, move their members to other subsets
name_map = Dict([taxon_name => i for (i, taxon_name) in enumerate(namelist0)])
L = [length(s) for s in final_subsets];
while minimum(L) < 5
    global D0
    global final_subsets
    global name_map
    global L

    L = [length(s) for s in final_subsets];
    min_subset_idx = findmin(L)[2]
    min_subset = final_subsets[min_subset_idx]
    deleteat!(final_subsets, min_subset_idx)
    L = [length(s) for s in final_subsets];

    valid_subsets = final_subsets[findall(l -> l + length(min_subset) < 5 || l == minimum(L), L)]

    avg_dists = [mean(D0[name_map[taxa], name_map[valid_member]] for valid_member in v_subset for taxa in min_subset) for v_subset in valid_subsets]
    min_avg_dist_idx = findmin(avg_dists)[2]

    for tax in min_subset
        push!(valid_subsets[min_avg_dist_idx], tax)
    end
end

L = [length(s) for s in final_subsets];
@info "Subset decomposition: min=$(minimum(L)), max=$(maximum(L)), median=$(median(L)), mean=$(mean(L))"

open(output, "w+") do f
    for (si, s) in enumerate(final_subsets)
        for (iitem, item) in enumerate(s)
            if iitem == length(s)
                write(f, "$(item)")
            else
                write(f, "$(item),")
            end
        end
        if si != length(final_subsets)
            write(f, "\n")
        end
    end
end
