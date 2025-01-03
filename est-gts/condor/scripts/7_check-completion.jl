include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")

using CSV, DataFrames
df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", DataFrame)



all_string = ""
all_dag_string = ""
some_string = ""
none_string = ""
has_no_dag_string = ""
n_complete = 0
n_total = 0
for ntaxa in [50, 100, 200]#, 500, 1000]
    a_df = filter(r -> r.ntaxa == ntaxa, df)
    for rep in 1:25
        b_df = filter(r -> r.rep == rep, a_df)
        for ils in ["low", "high"]
            c_df = filter(r -> r.ils == ils, b_df)
            for ngt in [100, 1000, 3000]
                d_df = filter(r -> r.ngt == ngt, c_df)
                for m in [10, 20]
                    iter_df = filter(r -> r.m == m, d_df)
                    n_total += 1
                    global all_string, all_dag_string, some_string, none_string, has_no_dag_string, n_complete, n_total

                    # Totally complete
                    if nrow(iter_df) > 0
                        printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) [✓]\n", color = :green)
                        n_complete += 1
                        continue
                    end

                    pct_snaq = get_pct_complete_constraints(ntaxa, rep, ils, ngt, m)
                    # SNaQ is complete, but InPhyNet not run
                    if pct_snaq == 1.0
                        printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) [X]\n", color = :green)
                        n_complete += 1
                        all_string = "$(all_string)j 8_inphynet.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m);"
                        all_dag_string = "$(all_dag_string)j create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                        continue
                    end

                    # Partially complete
                    if pct_snaq > 0
                        some_string = "$(some_string)j create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                        printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) ($(round(100*pct_snaq, digits=2))%)\n", color = :cyan)
                        continue
                    end

                    none_string = "$(none_string)j create_submit_dag.jl $(ntaxa) $(rep) $(ils) $(ngt) $(m) true;"
                    printstyled("$(ntaxa) $(rep) $(ils) $(ngt) $(m) (0%)\n", color = :red)
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
@info "Command to submit DAGs for 0% analyses:"
printstyled("$(none_string)\n\n", color=:red)
