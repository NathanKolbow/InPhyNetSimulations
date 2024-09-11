@info "Loading packages"
using CSV, DataFrames, InPhyNet, PhyloNetworks


results = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/out.csv", DataFrame)
full_df = DataFrame(
    ntaxa=Int64[], rep=Int64[], ils=String[], ngt=Int64[], m=Int64[], subset_idx=Int64[], run_number=Int64[]
)
n500r12_df = DataFrame(
    ntaxa=Int64[], rep=Int64[], ils=String[], ngt=Int64[], m=Int64[], subset_idx=Int64[], run_number=Int64[]
)
n1000r1_df = DataFrame(
    ntaxa=Int64[], rep=Int64[], ils=String[], ngt=Int64[], m=Int64[], subset_idx=Int64[], run_number=Int64[]
)


@info "Checking parameters"
astral_folder = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/astral/"
snaq_folder = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/snaq/"
for ntaxa in [500, 1000]
    for rep in ifelse(ntaxa == 500, [1, 2], [1])
        for ils in ["low", "high"]
            for ngt in [100, 1000]
                for m in [10, 20]
                    iter_file = joinpath(astral_folder, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile")
                    
                    if isfile(iter_file) && length(readlines(iter_file)) > 0
                        if nrow(filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, results)) > 0
                            printstyled("[COMPLETE] ", color = :green)
                            println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")
                            continue
                        end

                        tre0 = readTopology(iter_file)
                        n_subsets = length(sateIdecomp(tre0, m))

                        partial = false
                        for s = 1:n_subsets
                            for j = 1:10
                                run_file = joinpath(snaq_folder, "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(s)-r$(j).runtime")
                                if isfile(run_file)
                                    partial = true
                                    continue
                                end

                                push!(full_df, [ntaxa, rep, ils, ngt, m, s, j])
                                if ntaxa == 500 && rep in 1:2
                                    push!(n500r12_df, [ntaxa, rep, ils, ngt, m, s, j])
                                elseif ntaxa == 1000 && rep == 1
                                    push!(n1000r1_df, [ntaxa, rep, ils, ngt, m, s, j])
                                end
                            end
                        end

                        if partial
                            printstyled("[PARTIAL] ", color = :yellow)
                        else
                            printstyled("[INCOMPLETE] ", color = :red)
                        end
                        println("n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)")
                    end
                end
            end
        end
    end
end


n500_r12_out_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/submit/snaq_n500r12.tab"
n1000r1_out_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/submit/snaq_n1000r1.tab"
full_out_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/submit/snaq.tab"

@info "Writing $(nrow(full_df)) lines to snaq.tab"
@info "Writing $(nrow(n500r12_df)) lines to snaq_n500r12.tab"
@info "Writing $(nrow(n1000r1_df)) lines to snaq_n1000r1.tab"

CSV.write(full_out_file, full_df, writeheader = false)
CSV.write(n500_r12_out_file, n500r12_df, writeheader = false)
CSV.write(n1000r1_out_file, n1000r1_df, writeheader = false)