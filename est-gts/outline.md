## Pipeline outline

1. Simulate true gene trees (`manually-julia`)
2. Simulate sequences from gene trees (`manually-bash`)
3. Estimate gene trees (`condor-bash`)
4. Infer the starting tree with Astral (`manually-bash`)
5. Subset decomposition (`manually-julia`)
6. Infer networks with SNaQ (`condor-julia`)
7. InPhyNet (`manually-julia`)


## Simulation parameters

- ntaxa:        500, 1000
- ils:          low, high   (potentially `med` later)
- replicates:   1           (expand to 10 reps per ntaxa after initial results)
- ngt:          100, 1000   (potentially more/fewer later depending on results)
- m:            22


## Example

Here is an example of how to get all the data for `n=500`, replicate #1 w/ 100 gts, low ILS and `m=22`:

```bash
# 1.
est-gts/condor/scripts/julia 1_true-gts.jl

# 2.
est-gts/condor/scripts/2_seqs.sh 500 1 low 100 22

# 3.
python3 est-gts/condor/submits/iqtree.tab.py        # Setup this file accordingly first
condor_submit est-gts/condor/submits/iqtree.submit

# 4.
est-gts/condor/scripts/4_astral.sh 500 1 low 100 22

# 5.
julia est-gts/condor/scripts/5_subset_decomp.jl 500 1 low 100 22

# 6.
python3 est-gts/condor/submits/snaq.tab.py          # Setup this file accordingly first
condor_submit est-gts/condor/submits/snaq.submit

# 7.
# not written yet
```