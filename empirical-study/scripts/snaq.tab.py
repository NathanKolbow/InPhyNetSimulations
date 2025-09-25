import os
from os.path import isfile

subset_dir = f"/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/empirical-study/subset_data/"
n_written = 0

with open("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/empirical-study/scripts/snaq.tab", "w") as tab:
    for which_group in ["gymnosperms"]:
        for subset_number in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
            for h in [0, 1, 2, 3]:
                for run_number in range(20):
                    run_number += 1
                    iter_dir = f"{subset_dir}/{which_group}/subset{subset_number}/"
                    
                    # If was completed when we were doing individual runs
                    filepath = f"{iter_dir}/snaq_s{subset_number}_h{h}_run{run_number}.out"
                    if isfile(filepath) and os.stat(filepath).st_size > 0:
                        continue
                    
                    n_written += 1
                    tab.write(f"{which_group},{subset_number},{h},{run_number}\n")
                    
print(f"Wrote {n_written} lines to snaq.tab")
