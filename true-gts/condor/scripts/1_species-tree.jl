@info "Loading project"
using Pkg
Pkg.activate("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/")

@info "Loading helpers"
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])

@info "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)"


astral_jar = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/software/Astral/astral.5.7.1.jar"

astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/astral/"
astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

true_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/true-gts/"
true_gt_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"


truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

if !isfile(astral_file) || length(readlines(astral_file)) == 0
    set_seed!(ntaxa, rep, ils, ngt, m, 0)
    gts = simulatecoalescent(truenet, ngt, 1)

    # Save gts to file so that Astral can use them
    @info "Writing true gts to file"
    writeMultiTopology(gts, true_gt_file)

    # Run Astral
    @info "Running Astral"
    silently() do 
        run(`java -jar $(astral_jar) -i $(true_gt_file) -o $(astral_file) -t 0 -s 0`)
    end
    rm(true_gt_file)
end

tre0 = readTopology(astral_file)
maj = majorTree(truenet)
rootatnode!(tre0, "OUTGROUP")
rootatnode!(maj, "OUTGROUP")

hwcd = hardwiredClusterDistance(maj, tre0, true)
@info "HWCD: $(hwcd)"


# Subset decomposition
subsets = sateIdecomp(tre0, m)

small_idxs = findall(set -> length(set) < 5, subsets)
if length(small_idxs) > 0
    throw(ErrorException("Subset sizes: $(sort([length(s) for s in subsets]))"))
end

tab_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/submit/snaq.tab"
tab_lines = readlines(tab_file)
if "$(ntaxa),$(rep),$(ils),$(ngt),$(m),$(length(subsets))" in tab_lines
    @info "Lines already written to snaq.tab"
    exit()
end

open(tab_file, "a+") do file
    for subset_idx = 1:length(subsets)
        for run_number = 1:10
            write(file, "$(ntaxa),$(rep),$(ils),$(ngt),$(m),$(subset_idx),$(run_number)\n")
        end
    end
end
@info "Wrote $(10 * length(subsets)) lines to snaq.tab"