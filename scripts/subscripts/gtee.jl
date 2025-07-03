
using PhyloNetworks, StatsBase


function gtee(gts1::Vector{HybridNetwork}, gts2::Vector{HybridNetwork})
    gtees::Vector{Float64} = []
    for (gt1, gt2) in zip(gts1, gts2)
        push!(
            gtees,
            hardwiredclusterdistance(gt1, gt2, false) / (2 * gt1.numtaxa - 6)
        )
    end
    return mean(gtees)
end
