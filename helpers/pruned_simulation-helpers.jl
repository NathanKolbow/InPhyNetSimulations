using PhyloNetworks, PhyloCoalSimulations




function run_pruned_simulation(truenet::HybridNetwork, min_m::Int64, max_m::Int64, ngt::Int64, d_metric::Function, repititions::Int64; verbose::Bool = false)

    rootatnode!(truenet, "OUTGROUP")

    # 1. Generate gene trees
    if verbose @info "Simulating gene trees" end
    gts = simulatecoalescent(truenet, ngt, 1)

    # 2. Calculate distance matrix
    if verbose @info "Calculating distance matrix" end
    D, namelist = d_metric(gts)

    # 3. NJ tree
    if verbose @info "Computing NJ tree" end
    tre0 = nj(DataFrame(D, namelist))

    mnet = nothing
    total_inphynet_time = 0.
    for rep = 1:repititions
        # 4. Subset decomposition
        if verbose @info "[$(rep)/$(repititions)] Computing subset decomposition" end
        rootatnode!(tre0, "OUTGROUP")
        subsets = sateIdecomp(tre0, min_m, max_m)

        # 5. Prune networks
        if verbose @info "[$(rep)/$(repititions)] Pruning constraint networks" end
        cs = pruneTruthFromDecomp(truenet, subsets)

        # 6. InPhyNet
        if verbose @info "[$(rep)/$(repititions)] Running InPhyNet" end
        total_inphynet_time += @elapsed mnet = inphynet(D, cs, namelist)
        if rep < repititions tre0 = majorTree(mnet) end
    end

    # 7. Calculate metrics
    common_root!([mnet, truenet, tre0], "OUTGROUP")
    maj_mnet = majorTree(mnet)
    maj_truenet = majorTree(truenet)

    rootatnode!(maj_mnet, "OUTGROUP")
    rootatnode!(maj_truenet, "OUTGROUP")

    hwcd = hardwiredClusterDistance(mnet, truenet, true)
    maj_hwcd = hardwiredClusterDistance(maj_mnet, maj_truenet, false)
    nj_hwcd = hardwiredClusterDistance(tre0, maj_truenet, false)

    nretics_est = mnet.numHybrids
    nretics_true = truenet.numHybrids

    return total_inphynet_time, hwcd, maj_hwcd, nj_hwcd, nretics_est

end



function run_pruned_simulation_suite(ntaxas, reps, ilss, min_ms, max_ms, ngts, repititionss, d_metrics, nsim, dat_file)
    n_combos = length(min_ms) * length(max_ms) * length(ngts) * length(repititionss) * length(d_metrics)

    for ntaxa in ntaxas
        for rep in reps
            for ils in ilss
                truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)
                combos_run = 0

                for min_m in min_ms
                    for max_m in max_ms
                        for ngt in ngts
                            for repititions in repititionss
                                for d_metric in d_metrics
                                    d_metric_string = metric_string(d_metric)
                                    combos_run += 1

                                    nretics_true = truenet.numHybrids
                                    runtime = Array{Float64}(undef, nsim)
                                    hwcd = Array{Float64}(undef, nsim)
                                    maj_hwcd = Array{Float64}(undef, nsim)
                                    nj_hwcd = Array{Float64}(undef, nsim)
                                    nretics_est = Array{Int64}(undef, nsim)

                                    for j = 1:nsim
                                        if Threads.threadid() == 1
                                            print("\rn$(ntaxa) #$(rep) $(ils) - $(combos_run) / $(n_combos) [$(j-1) / $(nsim)]")
                                        end

                                        try
                                            runtime[j], hwcd[j], maj_hwcd[j], nj_hwcd[j], nretics_est[j] = 
                                                run_pruned_simulation(truenet, min_m, max_m, ngt, d_metric, repititions, verbose = false)
                                        catch
                                            continue
                                        end

                                        print("\rn$(ntaxa) #$(rep) $(ils) - $(combos_run) / $(n_combos) [$(j) / $(nsim)]")
                                    end
                                        
                                    CSV.write(dat_file, DataFrame(
                                        ntaxa=ntaxa, rep=rep, ils=ils, min_m=min_m, max_m=max_m, ngt=ngt, d_metric=d_metric_string,
                                        hwcd=hwcd, major_hwcd=maj_hwcd, nj_hwcd=nj_hwcd, nretics_est=nretics_est
                                    ), append = true)
                                end
                            end
                        end
                    end
                end

                println("")
            end
        end
    end

end
metric_string(fxn::Function) = ifelse(fxn == calculateAGID, "D", "C")


function common_root!(nets::Vector{HybridNetwork}, preferred_root::String)

    try
        rootatnode!(mnet, preferred_root)
        rootatnode!(truenet, preferred_root)
        rootatnode!(tre0, preferred_root)
    catch
    end

    for t_name in tipLabels(nets[1])
        try
            rootatnode!(mnet, t_name)
            rootatnode!(truenet, t_name)
            rootatnode!(tre0, t_name)
            return
        catch
        end
    end

    error("Could not find common root")
end