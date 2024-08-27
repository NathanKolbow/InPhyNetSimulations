using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")



function sim_true_gts(ntaxa::Int64, rep::Int64, ils::String, ngt::Int64, m::Int64)
    unsplit_output_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/"
    unsplit_output_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"
    split_output_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/split/"
    split_output_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

    set_seed!(ntaxa, rep, ils, ngt, m, 0)       # 0 for the seed b/c that isn't used in any of the true gt sims
    gts = simulatecoalescent(truenet, ngt, 1)
    writeMultiTopology(gts, unsplit_output_file)
    for i = 1:ngt
        writeTopology(gts[i], split_output_file*"_$(i)")
    end
end

# Parameters from `est-gts/sim-outline.md` #

for ntaxa_param in [500, 1000]
    for ils_param in ["low", "high"]
        for rep_param in [1]
            for ngt_param in [100, 1000]
                for m_param in [22]

                    sim_true_gts(ntaxa_param, rep_param, ils_param, ngt_param, m_param)

                end
            end
        end
    end
end


