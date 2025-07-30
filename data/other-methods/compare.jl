include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/other-methods/newick-conversion.jl")

out_df = DataFrame(
    ntaxa = Int[], rep = Int[], ils = String[], ngt = Int[],
    mnet10_hwcd = Float64[], mnet20_hwcd = Float64[], consensus_hwcd = Float64[], neighbor_hwcd = Float64[]
)
mnet_df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/analysis/approx_normalized_errors.csv", DataFrame)

for ngt in [1000]
    @info ngt
    for ils in ["low", "high"]
        @info "\t$ils"
        consensus_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/other-methods/consensus-nets/n200-$(ils)-$(ngt)gt-m20.netfile"
        nnet_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/other-methods/neighbor-net/n200-$(ils)-$(ngt)gt-m20.netfile"

        consensus_nets = [readTopology(remove_nested_markers(line)) for line in readlines(consensus_file)]
        nnet_nets = [readTopology(remove_nested_markers(line)) for line in readlines(nnet_file)]

        for rep = 1:15
            cs10 = get_constraints(200, rep, ils, ngt, 10)
            if cs10 === nothing continue end
            cs20 = get_constraints(200, rep, ils, ngt, 20)
            if cs20 === nothing continue end

            @info "\t\t$(rep)"

            mnet_rows = filter(r -> r.ntaxa == 200 .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt, mnet_df)
            mnet10_hwcd = -1
            try mnet10_hwcd = filter(r -> r.m == 10, mnet_rows)[1, "approx_norm_hwcd"] catch end
            mnet20_hwcd = -1
            try mnet20_hwcd = filter(r -> r.m == 20, mnet_rows)[1, "approx_norm_hwcd"] catch end

            truenet = load_true_net_ils_adjusted(200, rep, ils)
            
            consensus_hwcd = -1
            if length(consensus_nets) >= rep
                common_root!([truenet, consensus_nets[rep]], "OUTGROUP")
                consensus_hwcd = hardwiredClusterDistance(truenet, consensus_nets[rep], true)
                consensus_hwcd /= (truenet.numEdges - truenet.numTaxa + consensus_nets[rep].numEdges - consensus_nets[rep].numTaxa)
            end

            nnet_hwcd = -1
            if length(nnet_nets) >= rep
                common_root!([truenet, nnet_nets[rep]], "OUTGROUP")
                nnet_hwcd = hardwiredClusterDistance(truenet, nnet_nets[rep], true)
                nnet_hwcd /= (truenet.numEdges - truenet.numTaxa + nnet_nets[rep].numEdges - nnet_nets[rep].numTaxa)
            end

            push!(out_df, [
                200, rep, ils, ngt,
                mnet10_hwcd, mnet20_hwcd, consensus_hwcd, nnet_hwcd
            ])
            CSV.write("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/other-methods/comparison.csv", out_df)

        end
    end
end

