# Subset Decomposition Schema

The automatic subset decomposition schema works as follows:
1. Take rough "clade" and "ancestral" subsets (excluding the outgroup) of size 50 informed by NJ. For a network w/ 500 taxa + 1 outgroup there will be 20 subsets (10 for "clade", 10 for "ancestral").
2. Infer a tree of blobs on each of the subsets in (1).
3. For each blob in each of the inferred trees of blobs, gather the "cycle nodes" that form that blob's cycle, then create a set of taxa that includes 1 descendant of each cycle node
4. Perform normal SATe-I decomposition, but ensure that the sets generated in (3) always stay together. How we do this:
    - SATe-I will be performed on a species tree inferred by ASTRAL, call it $T$.
    - For each "blob-subset" $s_i^B$ with size $n_{s_i^B}=|s_i^B|$, remove $S[2:n_{s_i^B}]$ from $T$
    - Perform SATe-I subset decomposition to obtain subsets $S^*_i$
    - For each "blob-subset" $s_i^B$, find the SATe-I subset that contains $s_i^B[1]$ and add $s_i^B[2:n_{s_i^B}]$ to it
    - This is a very rough procedure that may cause subsets to be larger than the specified $m$, but it works

# Implementation

The schema above is implemented as follows:
1. `5a_subset_decomp.jl` takes typical CMD args and performs step (1) and writes the pruned gene trees to `data/subsets/<sim ID>/clade<index>.tre` and `data/subsets/<sim ID>/ancestral<index>.tre`
2. `5b_subset_decomp.R` takes typical CMD args **PLUS** two additional args: `<clade/ancestral>` and `<index>`, then reads the respective gene trees from `data/subsets/<sim ID>/<clade/ancestral><index>.tre`, performs step (2) and writes the inferred tree of blobs for that subset to `data/subsets/<sim ID>/<clade/ancestral><index>.tob`
3. Once all ToBs have been inferred, `5c_subset_decomp.jl` performs steps (3) and (4), then, for each subset, writes the pruned ASTRAL tree and pruned estimated gene trees to `data/subsets/<sim ID>/ASTRAL<index>.tre` and `.../pruned_gts<index>.tre`, respectively