import os.path


n_written = 0
with open("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/submits/snaq.tab", "w") as tab:
    for ntaxa in [500, 1000]:
        for rep in [1]:
            for ils in ["low", "high"]:
                for ngt in [100, 1000]:
                    for m in [10, 20]:
                        for subset_idx in range(1, ngt+1, 1):
                            snaq_iter_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/snaq"
                            iter_dir = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/subsets"
                            iter_dir = f"{iter_dir}/n{ntaxa}-r{rep}-{ils}-{ngt}gt-m{m}-subset{subset_idx}"
                            
                            if not os.path.isdir(iter_dir):
                                break
                            
                            if not os.path.isfile(f"{iter_dir}/pruned_gts.treefile"):
                                break
                            
                            if not os.path.isfile(f"{iter_dir}/tre0.treefile"):
                                break
                            
                            for run_number in range(1, 10+1, 1):
                                if os.path.isfile(f"{snaq_iter_dir}/n{ntaxa}-r{rep}-{ils}-{ngt}gt-m{m}-subset{subset_idx}-run{run_number}.runtime"):
                                    continue
                                
                                with open(f"{iter_dir}/pruned_gts.treefile", "r") as f:
                                    if len(f.readlines()) != ngt:
                                        continue
                                
                                with open(f"{iter_dir}/tre0.treefile", "r") as f:
                                    if len(f.readlines()) != 1:
                                        continue
                                
                                tab.write(f"{ntaxa},{rep},{ils},{ngt},{m},{subset_idx},{run_number}\n")
                                n_written += 1


print(f"Wrote {n_written} lines to snaq.tab.py")