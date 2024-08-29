## Pipeline outline

1. Starting tree with Astral (`manually-julia`)
    1. Generate true gene trees
    2. Infer species tree with Astral
    3. Save Astral tree
    4. Perform subset decomposition and, `for i=1:n_subsets`, append `<network-info>,i` to `true-gts/condor/submit/snaq.tab`

2. Infer networks with SNaQ (`condor-julia`)
    1. Load Astral tree
    2. Do subset decomposition
    3. Load true gene trees, pruning them appropriately as they come in
    4. Infer constraint with SNaQ
    5. **Note:** Runs are done individually instead of 10 at a time so that things finish in Condor more efficiently.

3. InPhyNet (`manually-julia`)


## Simulation parameters

- ntaxa:        500[, 1000, 2500]   **(only running 500 right now, queue 1000 & 2500 later)**
- ils:          low, high           (potentially `med` later)
- replicates:   1-10                (potentially expand to more replicates later)
- ngt:          100, 1000           (potentially more/fewer later depending on initial results)
- m:            10, 20              (potentially more params later)


## Running the sims

Start to finish for `n=500`, `rep=1`, `ils=low`, `ngt=100`, `m=10`

```bash
# 1.
julia scripts/1_species-tree.jl 500 1 low 100 10

# 2.
condor_submit submits/snaq.submit

# 3.
# not implemented yet
```