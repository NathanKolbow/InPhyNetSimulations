import sys
import physquirrel as psq

def infer_constraint(subset_taxa, msa):
    # Trim the MSA
    trimmed_msa = msa.deep
    
    Q = trimmed_msa.delta_heuristic()
    return Q.squirrel()

if __name__ == "__main__":
    subsetfile = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    
    full_msa = psq.MSA(fasta)
    networks = []
    
    with open(subsetfile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            N = infer_constraint(line.split(","), full_msa)
            networks.append(N)
    
    with open(output, "w+") as f:
        for N in networks:
            f.write(N.create_enewick())
            f.write("\n")