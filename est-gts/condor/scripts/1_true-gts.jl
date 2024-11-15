ENV["LD_LIBRARY_PATH"] = "/mnt/dv/wid/projects1/WID-Software/rocky8/x86_64/R-4.4.0/bin/"


ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])

out_pref = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/"
if all(isfile(out_pref * "$(i).treefile") for i=1:ngt)
    @info "Already complete - quitting."
    exit(0)
end

include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")



function sim_true_gts(ntaxa::Int64, rep::Int64, ils::String, ngt::Int64, m::Int64)
    @info "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)"
    
    output_file_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/true-gts/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/"
    if !isdir(output_file_prefix) mkdir(output_file_prefix) end

    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

    set_seed!(ntaxa, rep, ils, ngt, m, 0)       # 0 for the seed b/c that isn't used in any of the true gt sims
    gts = simulatecoalescent(truenet, ngt, 1)

    for i = 1:ngt
        writeTopology(gts[i], output_file_prefix*"$(i).treefile")
    end
end

# Parameters from `est-gts/sim-outline.md` #

sim_true_gts(ntaxa, rep, ils, ngt, m)


# for ntaxa_param in [500, 1000]
#     for ils_param in ["low", "high"]
#         for rep_param in 1:5
#             for ngt_param in [100, 1000]
#                 for m_param in [10, 20]

#                     sim_true_gts(ntaxa_param, rep_param, ils_param, ngt_param, m_param)

#                 end
#             end
#         end
#     end
# end


