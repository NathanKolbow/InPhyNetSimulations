# Helper functions for simulating real data
# Uses software in directory `network-merging/software/`
using PhyloNetworks, PhyloCoalSimulations, StaticArraysCore, StaticArrays

function simulate_sequence_data(gts::Vector{HybridNetwork}, output_file_prefix::String, estgt_file::String="", data_dir::String="")
    seq_file_paths = ["$(output_file_prefix)_$(i)" for i=1:length(gts)]
    if all(isfile(f) for f in seq_file_paths)
        @debug "Sequence data files already exist"
        return seq_file_paths
    elseif isfile(estgt_file)
        @debug "Estimated gene trees already exist, skipping sequence generation."
        return seq_file_paths
    end

    # find the correct `-s` flag to use w/ seq-gen
    # s = find_seqgen_s(gts)
    s = 0.036
    

    # simulate sequences w/ seq-gen
    seq_file_paths = run_seqgen_multi(s, gts, output_file_prefix, data_dir = data_dir)
    return seq_file_paths
end


function estimate_gene_trees(seq_files::Vector{String}, estgt_output_file::String)
    if isfile(estgt_output_file)
        @debug "Estimated gene tree file already exists"
        return
    end
    
    estgt_newicks = Array{String}(undef, length(seq_files))
    Threads.@threads for i=1:length(seq_files)
        seq_file = seq_files[i]

        # run iqtree
        run_iqtree(seq_file)

        # read iqtree results
        estgt_newicks[i] = readlines("$(seq_file).treefile")[1]
        clean_est_gt_files("fake_placeholder", seq_file)
    end

    # save iqtree results to file
    open(estgt_output_file, "w+") do f
        for newick in estgt_newicks
            write(f, "$(newick)\n")
        end
    end
end


function silently(f)
    redirect_stdout(Base.DevNull()) do
        redirect_stderr(Base.DevNull()) do 
            return f()
        end
    end
end


"""
Finds how many reticulations are in the subnetwork of `true_net` composed of the taxa in `subset_taxa`.
Used to determine how many reticulations should be inferred w/ SNaQ.
"""
function retics_in_subnet(true_net, subset_taxa)
    subnet = pruneTruthFromDecomp(true_net, subset_taxa)
    return subnet.numHybrids
end


function load_true_net_ils_adjusted(ntaxa::Int64, replicatenum::Int64, ils::String)
    truenet = readMultiTopology(get_network_filepath(ntaxa))[replicatenum]
    newick = writeTopology(truenet)
    avg_bl = get_avg_bl(truenet)
    newick = "($(newick[1:(length(newick)-1)]):$(avg_bl),OUTGROUP:1.0);"
    truenet = readTopology(newick)

    # ils level : desired average branch length
    # - low       : 2.0
    # - med       : 1.0
    # - high      : 0.5
    # - very high : 0.1

    # minimum branch length in low ILS case BEFORE overall adjustment: 0.5
    # maximum branch length in low ILS case BEFORE overall adjustment: 5.0

    # 1. make the network have low ILS branch lengths. after adjusting for ILS, we extend branch
    #    lengths such that root to tip distance is unchanged, so this needs to be our baseline
    desired_avg = 2.
    avg_bl = get_avg_bl(truenet)
    for e in truenet.edge
        e.length *= desired_avg / avg_bl
    end

    # 1b. increase branch length min to 0.5 and decrease max to 5.0
    for e in truenet.edge
        if !e.hybrid || !e.isMajor
            e.length = max(e.length, 0.5)
            e.length = min(e.length, 5.0)
        end
    end

    # 1c. reset to desired average again, now with increased minimum lengths
    desired_avg = 2.
    avg_bl = get_avg_bl(truenet)
    for e in truenet.edge
        e.length *= desired_avg / avg_bl
    end

    # 2. record desired final depths
    leaf_depths = get_leaf_depths(truenet)

    # 3. adjust network to desired branch length avg for desired level of ILS
    desired_avg = ils == "low" ? 2. : (ils == "med" ? 1. : (ils == "high" ? 0.5 : 0.1))
    avg_bl = get_avg_bl(truenet)
    for e in truenet.edge
        e.length *= desired_avg / avg_bl
    end
    
    # 4. extend leaves to desired length
    extend_leaves!(truenet, leaf_depths)

    # 5. adjust gammas to be randomly drawn from exp(1 / 7)
    Random.seed!(Int64(ils[1]))     # make the gamma parameter constant for a network given ILS
    for hyb in truenet.hybrid
        g = rand(Exponential(1. / 7.))
        g = 1 - min(g, 0.5) # this is major parent gamma, so it is at least 0.5
        
        getparentedge(hyb).gamma = g
        getparentedgeminor(hyb).gamma = 1 - g
    end

    return truenet
end



### Helper functions for the main functions above ###

# Extends the branches above leaves in `net` such that the root to tip distance
# matches that in `desired_depths`.
function extend_leaves!(net::HybridNetwork, desired_depths::AbstractVector{Float64})
    curr_depths = get_leaf_depths(net)
    for (i, adder) in enumerate(desired_depths .- curr_depths)
        edge = getparentedge(net.leaf[i])
        edge.length += adder
    end
end


# Finds the root to tip distance for each tip in net's major tree.
function get_leaf_depths(net::HybridNetwork)
    tre = majorTree(net)
    depths = zeros(net.numTaxa)

    for (i, leaf) in enumerate(net.leaf)
        curr = leaf
        while curr != net.node[net.root]
            depths[i] += getparentedge(curr).length
            curr = getparent(curr)
        end
    end

    return depths
end


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


# Finds the `-s` parameter for seq-gen such that gene
# tree estimation error is in a reasonable range
#
# Returns: path to sequence file
function find_seqgen_s(gts::Vector{HybridNetwork}; data_dir::AbstractString="")
    error("DEPRECATED")

    min_s = 0.00001
    max_s = 0.1
    curr_s = 0.003

    desired_gtee = 0.25
    tolerance = 0.01

    ntrees = length(gts)
    gtees = zeros(ntrees)

    while true
        Threads.@threads for i=1:ntrees
            tree = gts[i]
            true_newick = writeTopology(tree)

            temp_seqfile = joinpath(data_dir, "seqgen_$(i).phy")
            temp_gtfile = joinpath(data_dir, "truegt_$(i).treefile")
            writeTopology(tree, temp_gtfile)

            # seq-gen
            run_seqgen(curr_s, temp_gtfile, temp_seqfile)

            # iqtree
            run_iqtree(temp_seqfile)

            # save the result
            est_newick = readlines("$(temp_seqfile).treefile")[1]

            # calculate gtee
            gtee_nrf = calc_gtee(true_newick, est_newick)
            gtees[i] = parse(Float64, gtee_nrf)

            # clean up
            clean_est_gt_files(temp_gtfile, temp_seqfile)
            rm_suppress(temp_seqfile)
        end

        avg_gtee = mean(gtees)
        if avg_gtee - desired_gtee > tolerance
            min_s = curr_s
            curr_s = mean([max_s, min_s])
        elseif avg_gtee - desired_gtee < -tolerance
            max_s = curr_s
            curr_s = mean([max_s, min_s])
        else
            return curr_s
        end
        curr_s = round(curr_s, sigdigits = 4)
    end
end


function run_seqgen(seqgen_s::AbstractFloat, temp_gtfile::AbstractString, temp_seqfile::AbstractString; seq_length::Int64=1000)
    error("DEPRECATED")
    #software_path = "/m nt/dv/wid/pro jects4/SolisLemus-network- merging/software/"
    if Sys.isapple()
        run(pipeline(`$software_path/seq-gen-macos -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    elseif Sys.islinux()
        run(pipeline(`$software_path/seq-gen-linux -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    else
        run(pipeline(`$software_path/seq-gen-windows.exe -s$seqgen_s -n1 -f0.3,0.2,0.2,0.3 -t3.0 -mHKY -op -l$(seq_length) $temp_gtfile`, stdout=temp_seqfile, stderr=devnull))
    end
end


# `truegt_file` contains several topologies, and we need to run seq-gen for each of them
function run_seqgen_multi(seqgen_s::AbstractFloat, gts::Vector{HybridNetwork}, output_file_prefix::AbstractString; seq_length::Int64=1000, data_dir::AbstractString="")
    error("DEPRECATED")
    ntrees = length(gts)
    seq_file_paths = Array{String}(undef, ntrees)
    Threads.@threads for i=1:ntrees
        temp_gt_file = joinpath(data_dir, "tempgt_$(i).treefile")
        writeTopology(gts[i], temp_gt_file)

        temp_seq_file = joinpath(data_dir, "tempseq_out_$(i).phy")

        run_seqgen(seqgen_s, temp_gt_file, temp_seq_file, seq_length = seq_length)
        mv(temp_seq_file, "$(output_file_prefix)_$(i)")

        seq_file_paths[i] = "$(output_file_prefix)_$(i)"
        clean_est_gt_files(temp_gt_file, temp_seq_file)
    end
    return seq_file_paths
end


function run_iqtree(temp_seqfile::AbstractString)
    error("DEPRECATED")
    software_path = ""
    if Sys.isapple()
        try
            run(pipeline(`$software_path/iqtree-1.6.12-macos -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-macos -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    elseif Sys.islinux()
        try
            run(pipeline(`$software_path/iqtree-1.6.12-linux -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-linux -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    else
        try
            run(pipeline(`$software_path/iqtree-1.6.12-windows.exe -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        catch e
            run(pipeline(`$software_path/iqtree-1.6.12-windows.exe -redo -quiet -s $temp_seqfile`, stdout=devnull, stderr=devnull))
        end
    end
end


function calc_gtee(true_newick::AbstractString, est_newick::AbstractString)
    scriptpath = joinpath(@__DIR__, "..", "software", "compare_two_trees.py")
    gtee_nrf = Pipe()
    run(pipeline(`python3 $scriptpath -t1 $true_newick -t2 $est_newick`, stdout=gtee_nrf))
    close(gtee_nrf.in)
    return String(read(gtee_nrf))
end


function clean_iqtree_files(input_seqfile::AbstractString)
    rm_suppress(input_seqfile*".bionj")
    rm_suppress(input_seqfile*".ckp.gz")
    rm_suppress(input_seqfile*".iqtree")
    rm_suppress(input_seqfile*".log")
    rm_suppress(input_seqfile*".mldist")
    rm_suppress(input_seqfile*".model.gz")
    rm_suppress(input_seqfile*".treefile")
end


function clean_est_gt_files(temp_gtfile::AbstractString, temp_seqfile::AbstractString)
    rm_suppress(temp_gtfile)
    rm_suppress(temp_seqfile*".bionj")
    rm_suppress(temp_seqfile*".ckp.gz")
    rm_suppress(temp_seqfile*".iqtree")
    rm_suppress(temp_seqfile*".log")
    rm_suppress(temp_seqfile*".mldist")
    rm_suppress(temp_seqfile*".model.gz")
    rm_suppress(temp_seqfile*".treefile")
end

rm_suppress(file::AbstractString) = try rm(file) catch e end
rm_suppress(files::AbstractArray) = [rmsuppress(file) for file in files]