@info "Loading helpers"
using PhyloNetworks, InPhyNet
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


astral_jar = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/software/Astral/astral.5.7.1.jar"
function generate_tre0(ntaxa, rep, ils, ngt, m)
    global astral_jar


    astral_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/astral/"
    astral_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

    true_gt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/data/true-gts/"
    true_gt_file *= "n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"


    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

    if !isfile(astral_file) || length(readlines(astral_file)) == 0
        set_seed!(ntaxa, rep, ils, ngt, m, 0)
        gts = simulatecoalescent(truenet, ngt, 1)

        # Save gts to file so that Astral can use them
        writeMultiTopology(gts, true_gt_file)

        # Run Astral
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
    

    # Subset decomposition
    subsets = sateIdecomp(tre0, 5, m)

    small_idxs = findall(set -> length(set) < 5, subsets)
    if length(small_idxs) > 0
        @error "Minimum subset size: $(minimum([length(s) for s in subsets]))"
        return 0
    end

    tab_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/true-gts/condor/submit/snaq.tab"
    tab_lines = readlines(tab_file)
    if "$(ntaxa),$(rep),$(ils),$(ngt),$(m),$(length(subsets))" in tab_lines
        return 0
    end

    open(tab_file, "a+") do file
        for subset_idx = 1:length(subsets)
            for run_number = 1:10
                write(file, "$(ntaxa),$(rep),$(ils),$(ngt),$(m),$(subset_idx),$(run_number)\n")
            end
        end
    end

    return 10 * length(subsets)

end


if length(ARGS) == 5
    ntaxa = parse(Int64, ARGS[1])
    rep = parse(Int64, ARGS[2])
    ils = ARGS[3]
    ngt = parse(Int64, ARGS[4])
    m = parse(Int64, ARGS[5])
    generate_tre0(ntaxa, rep, ils, ngt, m)
elseif length(ARGS) > 0
    @error "Usage: julia 1_species-tree.jl [<ntaxa> <rep> <ils> <ngt> <m>]"
else
    lines_written = 0
    n_written = 0
    n_past = 0

    for ntaxa in [500, 1000]
        for rep = 1:5
            for ils in ["low", "high"]
                for ngt in [100, 1000]
                    for m in [10, 20]
                        global lines_written
                        global n_written
                        global n_past

                        n_past += 1
                        print("\rn$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)        [$(n_past) / 80]")

                        n_lines = generate_tre0(ntaxa, rep, ils, ngt, m)
                        lines_written += n_lines
                        n_written += length(n_lines) > 0
                    end
                end
            end
        end
    end
    @info "Wrote lines for $(n_written) parameter combinations, $(lines_written) total lines."
end