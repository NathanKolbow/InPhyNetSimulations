@info "Loading packages and helpers"
using InPhyNet
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")

if length(ARGS) < 6
    throw(ErrorException("Usage: julia driver.jl <ntaxa> <replicate> <ils> <ngt> <m> <nsim>"))
end

ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])
nsim = parse(Int64, ARGS[6])

out_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/data/true-gt-output/n$(ntaxa).csv"

################
# Saved values #
################
nretics_true = Array{Int64}(undef, nsim)
nretics_est = Array{Int64}(undef, nsim)
HWCDs = Array{Int64}(undef, nsim)
majortreeRFs = Array{Int64}(undef, nsim)
constraint_sums = Array{Int64}(undef, nsim)
pairwises = BitVector(undef, nsim)
runtimes = Array{Float64}(undef, nsim)
seeds = Array{Int64}(undef, nsim)
################

#################################################################################################
# We do a quick, fake run of `true_gt_simulation` to force all relevant functions to precompile #
#################################################################################################
@info "Forcing precompilation"
timetaken = @elapsed true_gt_simulation(load_true_net_ils_adjusted_level1(50, 1, "med"), 2, 30, 1)
@info "Took $(round(timetaken, digits=0)) seconds"
#################################################################################################


@info "Loading true network"
truenet = load_true_net_ils_adjusted_level1(ntaxa, rep, ils)


@info "Entering simulation loop"
start_time = time()
display_progress(0, nsim, start_time)

for i = 1:nsim
    # We can't actually do multi-threading here! Seeds are SHARED between threads,
    # so we cannot get reproducible results if we multithread this

    seeds[i] = set_seed!(ntaxa, rep, ils, ngt, m, i)
    HWCDs[i], majortreeRFs[i], constraint_sums[i], nretics_true[i], nretics_est[i], pairwises[i], runtimes[i] = true_gt_simulation(truenet, ngt, m, seeds[i])

    display_progress(i, nsim, start_time)
end


@info "Acquiring lockfile"
lk = lock_file(out_file)

try
    @info "Writing output to $(out_file)"
    CSV.write(
        out_file,
        DataFrame(
            ntaxa = repeat([ntaxa], nsim),
            nretics_true = nretics_true,
            nretics_est = nretics_est,
            max_subset_size = repeat([m], nsim),
            HWCD = HWCDs,
            majortreeRF = majortreeRFs,
            constraint_error_sum = constraint_sums,
            replicate_num = repeat([rep], nsim),
            estimated_pairwise = pairwises,
            runtime_seconds = runtimes,
            ngt = repeat([ngt], nsim),
            ils = repeat([ils], nsim),
            seeds = seeds
        ),
        append=true
    )
catch
finally
    @info "Releasing lockfile"
    close(lk)
end
