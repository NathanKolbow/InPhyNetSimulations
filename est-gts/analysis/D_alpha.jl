include("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/helpers/helpers.jl")

comp_params = get_completed_params()
train_set = [
    (ntaxa, rep, ils, ngt, m) for ntaxa in [200] for rep in 1:5 for ils in ["low", "high"] for ngt in [100, 1000, 3000] for m in [10, 20] if (ntaxa, rep, ils, ngt, m) in comp_params
]
test_set = [
    (ntaxa, rep, ils, ngt, m) for ntaxa in [200] for rep = 6:10 for ils in ["low", "high"] for ngt in [100, 1000, 3000] for m in [10, 20] if (ntaxa, rep, ils, ngt, m) in comp_params
]


alphas = 0.0:0.05:1.0
alpha_hwcd_sums = zeros(length(alphas))

for (ntaxa, rep, ils, ngt, m) in train_set
    print("\rLoading true net")
    truenet = load_true_net_ils_adjusted(ntaxa, rep, ils)

    print("\rLoading est gts ")
    gts = get_est_gts(ntaxa, rep, ils, ngt, m)
    print("\rCalculating D from gts")
    gt_D, gt_namelist = calculateAGID(gts)
    print("\rCalculating AGIC from gts")
    gt_AGIC, _ = calculateAGIC(gts)
    # print("\rCalculating D from truenet")
    # tre0_D, _ = internodedistance(majorTree(truenet), namelist = gt_namelist)
    print("\rGetting constraint networks    ")
    constraints = get_constraints(ntaxa, rep, ils, ngt, m)

    for (j, alpha) in enumerate(alphas)
        print("\rα = $(alpha) ($(j) / $(length(alphas)))                         ")
        D = gt_D .* alpha + gt_AGIC .* (1 - alpha)
        # D = gt_D .* alpha + tre0_D .* (1 - alpha)
        mnet = inphynet(D, constraints, gt_namelist)
        common_root!([truenet, mnet], "OUTGROUP")
        alpha_hwcd_sums[j] += hardwiredClusterDistance(truenet, mnet, true)
    end
    println()

    @info alpha_hwcd_sums
end

