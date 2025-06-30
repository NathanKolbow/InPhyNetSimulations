"""
Script adapted from `compare_tree_lists.py`, which was, quote:
Written by EKM (molloy.erin.k@gmail.com) in October 2016.
"""
import sys
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives

def compare_trees(tr1, tr2):
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)

    do = 0
    if com != lb1 or com != lb2:
        do = 1

    com = list(com)
    tns = dendropy.TaxonNamespace(com)

    if do:
        tr1.retain_taxa_with_labels(com)
        tr2.retain_taxa_with_labels(com)

        tr1.migrate_taxon_namespace(tns)
        tr2.migrate_taxon_namespace(tns)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)

    return(nl, ei1, ei2, fp, fn, rf)

if __name__ == "__main__":
    gtfile1 = sys.argv[1]
    gtfile2 = sys.argv[2]
    taxa = dendropy.TaxonNamespace()
    gtees = []
    
    with open(gtfile1, 'r') as f1:
        with open(gtfile2, 'r') as f2:
            i = 1
            for line1, line2 in zip(f1, f2):
                gt1 = dendropy.Tree.get(string=line1,
                                        schema="newick",
                                        rooting="force-unrooted",
                                        taxon_namespace=taxa)
                gt2 = dendropy.Tree.get(string=line2,
                                        schema="newick",
                                        rooting="force-unrooted",
                                        taxon_namespace=taxa)
                gtees.append(compare_trees(gt1, gt2)[5])
    print(gtees)
    gtee = sum(gtees) / len(gtees)
    print("GTEE: %.2f" % gtee)