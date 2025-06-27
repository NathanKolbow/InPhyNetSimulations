import time
import sys
import physquirrel as psq
from Bio import AlignIO
import Bio

def infer_constraint(subset_taxa, msa, tempdir):
    # Trim the MSA
    sub = [rec for rec in msa if rec.id in subset_taxa]
    submsa = Bio.Align.MultipleSeqAlignment(sub)
    AlignIO.write(submsa, "%s/trimmed.fasta" % tempdir, "fasta")
    
    trimmed_msa = psq.MSA("%s/trimmed.fasta" % tempdir)
    Q = trimmed_msa.delta_heuristic()
    return Q.squirrel()

if __name__ == "__main__":
    subsetfile = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    rtoutput = sys.argv[4]
    tempdir = sys.argv[5]
    
    full_msa = AlignIO.read(fasta, "fasta")
    networks = []
    runtimes = []
    
    with open(subsetfile, "r") as f:
        for line in f:
            start_time = time.perf_counter()
            line = line.rstrip("\n")
            N = infer_constraint(line.split(","), full_msa, tempdir)
            networks.append(N)
            runtimes.append(time.perf_counter() - start_time)
    
    with open(output, "w+") as f:
        for N in networks:
            f.write(N.create_enewick())
            f.write("\n")
    with open(rtoutput, "w+") as f:
        for rt in runtimes:
            f.write("%s" % rt)
            f.write("\n")
    