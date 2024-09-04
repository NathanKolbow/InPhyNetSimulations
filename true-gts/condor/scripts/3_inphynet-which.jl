using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


for ntaxa in [500, 1000]
    for ngt in [100, 1000]
        for ils in ["low", "high"]
            for rep in [1, 2, 3]
                for m in [10, 20]

                    astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/astral/"
                    astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

                    if !isfile(astral_file) || length(readlines(astral_file)) != 1 continue end
                    tre0 = readTopology(astral_file)
                    subsets = sateIdecomp(tre0, m)

                    n_subsets = length(subsets)

                    all_present = true
                    some_present = true
                    runs_done = 0
                    
                    for subset_idx = 1:n_subsets
                        runs = 0
                        for run_number = 1:10
                            snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/snaq/"
                            snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-r$(run_number)"

                            if !isfile("$(snaq_prefix).runtime")
                                all_present = false
                                continue
                            end
                            runs += 1
                        end
                        if runs == 0 some_present = false end
                        runs_done += runs
                    end

                    if all_present printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m)\n", color = :green)
                    elseif some_present printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :orange)
                    elseif runs_done > 0 printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :cyan) end
                        # else printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :red) end
                end
            end
        end
    end
end
