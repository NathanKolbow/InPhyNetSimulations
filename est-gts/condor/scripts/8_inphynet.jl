start_time = time()

ntaxa = parse(Int64, ARGS[1])
rep = parse(Int64, ARGS[2])
ils = ARGS[3]
ngt = parse(Int64, ARGS[4])
m = parse(Int64, ARGS[5])

# 0. Check if the DF already has an entry for these params
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/precompile-setup.jl")
using CSV, DataFrames

@info "Checking for existing record..."
df = CSV.read("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", DataFrame)
if nrow(filter(r -> r.ntaxa == ntaxa .&& r.rep == rep .&& r.ils == ils .&& r.ngt == ngt .&& r.m == m, df)) > 0
    @info "Entry already exists, quitting."
    exit()
end


@info (ntaxa, rep, ils, ngt, m)
@info "Loading helpers..."
include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")


# 1. Figure out how many subsets this network has
@info "Checking number of subsets..."
truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

subset_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m)/snaqtab"
estgt_file = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts/n$(ntaxa)-r$(rep)-$(ils)-$(ngt)gt-m$(m).treefile"

n_subsets = Int64(length(readlines(subset_file)) / 10)
@info "\t$(n_subsets) subsets"


# 2. Load data, all the while making sure everything actually exists
constraints = get_constraints(ntaxa, rep, ils, ngt ,m)
runtimes = get_runtimes(ntaxa, rep, ils, ngt, m)

cs_taxa = [c.numTaxa for c in constraints]
@info "Subset sizes: min=$(minimum(cs_taxa)), max=$(maximum(cs_taxa)), median=$(median(cs_taxa)), mean=$(mean(cs_taxa))"


# 4. Calculate D
@info "Calculating D..."
gts = readMultiTopology(estgt_file)
D, namelist = calculateAGID(gts)


# 5. Root all the constraints
for c in constraints
    while length(c.node[c.root].edge) != 2
        try
            rootatnode!(c, getchildren(c.node[c.root])[1])
        catch e
            rootatnode!(c, getchildren(c.node[c.root])[2])
        end
    end
end


# 5. Construct network with InPhyNet
@info "Running InPhyNet..."
inphynet_runtime = @elapsed mnet = inphynet(D, constraints, namelist)


# 6. Check and save results
common_root!([mnet, truenet], "OUTGROUP")

@info tipLabels(majorTree(mnet))
@info tipLabels(majorTree(truenet))

major_hwcd = hardwiredClusterDistance(majorTree(mnet), majorTree(truenet), true)
unrooted_major_hwcd = major_hwcd
if major_hwcd != 0
    global unrooted_major_hwcd
    @info "Calculating major unrooted HWCD:"
    if ntaxa <= 200
        unrooted_major_hwcd = hardwiredClusterDistance(majorTree(mnet), majorTree(truenet), false)
    else
        @error "SKIPPING UNROOTED CALCS FOR NOW"
    end
end


rooted_hwcd = hardwiredClusterDistance(mnet, truenet, true)
unrooted_hwcd = rooted_hwcd

if rooted_hwcd != 0
    global unrooted_hwcd
    @info "Calculating unrooted HWCD:"
    if ntaxa <= 200
        unrooted_hwcd = hardwiredClusterDistance(mnet, truenet, false)
    else
        @error "SKIPPING UNROOTED CALCS FOR NOW"
    end
end


@info "Retics: true=$(truenet.numHybrids), est=$(mnet.numHybrids), sum_constraints=$(sum([c.numHybrids for c in constraints]))"
@info "Calculating rooted min retic subset HWCD:"

rooted_min_greedy_hwcd, min_rooted_truenet = find_minimum_retic_subset_hwcd_greedy(truenet, mnet, verbose=true, swaponerror=true)
unrooted_min_greedy_hwcd = rooted_min_greedy_hwcd

if rooted_min_greedy_hwcd != 0
    global unrooted_min_greedy_hwcd
    @info "Calculating unrooted min retic subset HWCD:"
    if ntaxa <= 200
        unrooted_min_greedy_hwcd, _ = find_minimum_retic_subset_hwcd_greedy(truenet, mnet, verbose=true, swaponerror=true)
    else
        @error "SKIPPING UNROOTED CALCS FOR NOW"
    end
end

@info "Calculating minimum displayed tree NJ difference"
nj_tre = inphynet(D, Vector{HybridNetwork}([]), namelist)
min_unrooted_nj_hwcd = unrooted_major_hwcd
for disp_tre in displayedTrees(truenet, (ntaxa > 200 ? 0.48 : 0.0))
    global min_unrooted_nj_hwcd
    disp_hwcd = hardwiredClusterDistance(disp_tre, nj_tre, false)
    if disp_hwcd < min_unrooted_nj_hwcd
        min_unrooted_nj_hwcd = disp_hwcd
        if min_unrooted_nj_hwcd == 0 break end
    end
end
@info "Done w/ displayed trees"



nretic_true = truenet.numHybrids
nretic_est = mnet.numHybrids

snaq_runtime_sum = sum(runtimes)
snaq_runtime_serial = maximum(runtimes)

@info "Pruning constraints"
pruned_cs = simple_prune(truenet, [tipLabels(c) for c in constraints])
@info "Calculating constraint HWCDs"
sum_constraint_hwcd = sum([hardwiredClusterDistance(pruned_cs[i], constraints[i], false) for i = 1:length(constraints)])
sum_constraint_retics = sum(c.numHybrids for c in constraints)

perfect_infer_unrooted_hwcd = -100
try
    global perfect_infer_unrooted_hwcd
    if ntaxa <= 200
        perfect_infer_mnet = inphynet(D, pruned_cs, namelist)
        perfect_infer_unrooted_hwcd = hardwiredClusterDistance(truenet, perfect_infer_mnet, false)
    end
catch e
end


@show sum_constraint_hwcd
@show min_unrooted_nj_hwcd
@show min_unrooted_nj_hwcd + sum_constraint_hwcd

@show unrooted_hwcd
@show unrooted_min_greedy_hwcd

@show nretic_true
@show nretic_est
@show snaq_runtime_sum
@show snaq_runtime_serial

end_time = time()

println("Took $(round(end_time - start_time, digits=2)) seconds to calculate metrics.")


CSV.write("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv", 
    DataFrame(
        ntaxa=ntaxa, rep=rep, ils=ils, ngt=ngt, m=m,

        snaq_runtime_sum=snaq_runtime_sum,
        snaq_runtime_serial=snaq_runtime_serial,
        inphynet_runtime=inphynet_runtime,

        nretic_true=nretic_true, nretic_est=nretic_est,
        
        rooted_major_hwcd=major_hwcd,
        unrooted_major_hwcd=unrooted_major_hwcd,
        rooted_hwcd=rooted_hwcd, unrooted_hwcd=unrooted_hwcd,
        rooted_min_greedy_hwcd=rooted_min_greedy_hwcd, unrooted_min_greedy_hwcd=unrooted_min_greedy_hwcd,
        sum_constraint_hwcd=sum_constraint_hwcd,
        sum_constraint_retics=sum_constraint_retics,
        perfect_constraint_hwcd=perfect_infer_unrooted_hwcd,
        min_unrooted_nj_hwcd=min_unrooted_nj_hwcd,
        n_constraints=length(constraints)
    ), append=true
)



