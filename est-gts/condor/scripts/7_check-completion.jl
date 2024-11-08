include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")

using CSV, DataFrames
df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", DataFrame)



all_string = ""
all_dag_string = ""
some_string = ""
none_string = ""
less_than_none_string = ""
n_complete = 0
n_total = 0
for ntaxa in [50, 100] #, 200]#, 500, 1000]
    for rep in 1:10
        for ils in ["low", "high"]
            for ngt in [100, 1000]
                for m in [10, 20]
                    global all_string, all_dag_string, some_string, none_string, less_than_none_string, n_complete, n_total

                    n_total += 1

                    astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/astral/"
                    astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"
                    snaqtab_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab"

                    if !isfile(astral_file) || length(readlines(astral_file)) != 1 || !isfile(snaqtab_file)
                        less_than_none_string = "$(less_than_none_string)julia create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                        continue
                    end
                    n_subsets = Int(length(readlines(snaqtab_file)) / 10)

                    all_present = true
                    some_present = true
                    runs_done = 0
                    inphynet_inferred = nrow(filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, df)) > 0
                    
                    for subset_idx = 1:n_subsets
                        runs = 0
                        for run_number = 1:10
                            snaq_prefix = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq/"
                            snaq_prefix *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)-subset$(subset_idx)-run$(run_number)"

                            if !isfile("$(snaq_prefix).runtime")
                                all_present = false
                                continue
                            end
                            runs += 1
                        end
                        if runs == 0 some_present = false end
                        runs_done += runs
                    end
                    
                    if all_present printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) $(ifelse(inphynet_inferred, "[✓]", "[X]"))\n", color = :green)
                    elseif some_present printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :cyan)
                    elseif runs_done > 0 printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :yellow)
                    else printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*runs_done / (n_subsets*10), digits=2))%)\n", color = :red) end

                    if all_present
                        n_complete += 1
                        if !inphynet_inferred
                            all_string = "$(all_string)julia 8_inphynet.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m);"
                            all_dag_string = "$(all_dag_string)julia create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                        end
                    elseif some_present
                        if inphynet_inferred @error "[SEMI] Present in out.csv but not all SNaQ networks inferred: $(ntaxa) $(rep) $(ils) $(ngt) $(m)" end
                        some_string = "$(some_string)julia create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                    else
                        if inphynet_inferred @error "[NONE] Present in out.csv but not all SNaQ networks inferred: $(ntaxa) $(rep) $(ils) $(ngt) $(m)" end
                        none_string = "$(none_string)julia create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                    end
                end
            end
        end
    end
end


@info "$(n_complete)/$(n_total) settings finished"
@info "Command to run complete analyses missing from out.csv in DAGs:"
printstyled("$(all_dag_string)\n\n", color=:green)
@info "Command to run complete analyses missing from out.csv in tmux:"
printstyled("$(all_string)\n\n", color=:green)
@info "Command to submit DAGs for semi-complete analyses:"
printstyled("$(some_string)\n\n", color=:cyan)
@info "Command to submit DAGs for incomplete analyses:"
printstyled("$(none_string)\n\n", color=:yellow)
@info "Command to submit DAGs for analyses that are YET TO BE INITIALIZED:"
printstyled("$(less_than_none_string)\n\n", color=:red)

