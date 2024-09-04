import os.path

consolidated_path = "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/est-gts"
with open("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/submits/iqtree.tab", "w") as file:
    for ntaxa in [500, 1000]:
        for rep in [1]:
            for ils in ["low", "high"]:
                for ngt in [100, 1000]:
                    for m in [10, 20]:
                        
                        if os.path.isfile(f"{consolidated_path}/n{ntaxa}-r{rep}-{ils}-{ngt}gt-m{m}.treefile"):
                                continue
                        
                        for gt_idx in range(ngt):
                            gt_idx = gt_idx + 1
                            prio = 1000 / ngt
                            file.write(f"{ntaxa},{rep},{ils},{ngt},{m},{gt_idx},{prio}\n")