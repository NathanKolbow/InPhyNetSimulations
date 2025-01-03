"""
Robust comparison of two trees using dendropy

Written by EKM (molloy.erin.k@gmail.com) in October 2016.
"""
import argparse
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives


def compare_trees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives

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


def main(args):
    taxa = dendropy.TaxonNamespace()
    stre = dendropy.Tree.get(path=args.species,
                             schema="newick",
                             rooting="force-unrooted",
                             taxon_namespace=taxa)
    with open(args.output, 'w') as f:
        with open(args.gene, 'r') as gf:
            with open(args.esti, 'r') as ef:
                i = 1
                for line1, line2 in zip(gf, ef):
                    if (i >= 1) and (i <= 1000):
                        data = "exon"
                    elif (i >= 1001) and (i <= 2000):
                        data = "intron"
                    else:
                        data = "uce"
                    gtre = dendropy.Tree.get(string=line1,
                                             schema="newick",
                                             rooting="force-unrooted",
                                             taxon_namespace=taxa)
                    etre = dendropy.Tree.get(string=line2,
                                             schema="newick",
                                             rooting="force-unrooted",
                                             taxon_namespace=taxa)

                    [ad_nl, ad_e1, ad_e2,
                     ad_fp, ad_fn, ad_rf] = compare_trees(stre, gtre)
                    [ge_nl, ge_e1, ge_e2,
                     ge_fp, ge_fn, ge_rf] = compare_trees(gtre, etre)
                    [td_nl, td_e1, td_e2,
                     td_fp, td_fn, td_rf] = compare_trees(stre, etre)
                    f.write("%s,%s,%d,%d,%d,%d,%d,%d,%f,%d,%d,%d,%d,%d,%f,%d,%d,%d,%d,%d,%f\n" %
                            (args.prefix, data, i,
                             ad_nl, ad_e1, ad_e2, ad_fp, ad_fn, ad_rf,
                             ge_nl, ge_e1, ge_e2, ge_fp, ge_fn, ge_rf,
                             td_nl, td_e1, td_e2, td_fp, td_fn, td_rf))
                    i = i + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare tree lists")

    parser.add_argument("-s", "--species", type=str,
                        help="True species tree", required=True)

    parser.add_argument("-g", "--gene", type=str,
                        help="True gene tree list", required=True)

    parser.add_argument("-e", "--esti", type=str,
                        help="Estimated gene tree list", required=True)

    parser.add_argument("-p", "--prefix", type=str,
                        help="Prefix", required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file", required=True)

    main(parser.parse_args())
