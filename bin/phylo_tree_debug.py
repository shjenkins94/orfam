#!/usr/bin/env python

import argparse

from Bio import Phylo

__version__ = 1.0

INFO = """
phylo_tree - excludes non-OR genes on the basis of a phylogenetic tree.
When a given sequence is located in the outgroup clade, it should be discarded.
version: %s

Zuoyi Jian (jianzuoyi@gmail.com)
lab 518, Conservation Biology on Endangered Wildlife, SiChuan University
""" % __version__


def phylo_tree(args):
    with open(args.nwk, "rU") as in_nwk:
        tree = Phylo.read(in_nwk, "newick")
    outgroup = []
    orgroup = []

    with open(args.outgroup_id, "rU") as in_outgroup:
        non_or = in_outgroup.read().splitlines()

    for non_or_id in non_or:
        outgroup.append(tree.find_clades(non_or_id))

    terminals = tree.get_terminals()
    for terminal in terminals:
        if terminal not in outgroup:
            orgroup.append(terminal)

    # If genes with that in the outgroup share one MRCA, it may be a non-OR
    # gene and will be discarded.
    if not tree.is_monophyletic(outgroup):
        mrca = tree.common_ancestor(non_or)
        terminals = mrca.get_terminals()
        for terminal in terminals:
            if terminal in orgroup:
                print(terminal.name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=INFO)
    parser.add_argument("nwk",
                        help=("newick file"))
    parser.add_argument("outgroup_id",
                        help=("ids of outgroup sequences (txt)"))

    args = parser.parse_args()
    phylo_tree(args)
