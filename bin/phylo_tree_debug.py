#!/usr/bin/env python

import argparse

from Bio import Phylo

__version__ = 2.0

INFO = """
phylo_tree - has 2 modes depending on if optional argument classi_id is
present.
version: %s

Original program by Zuoyi Jian
Modified by Scott Jenkins
""" % __version__


# 2. Convert nwk to a tree and outgroup_id to a list.
def phylo_tree(args):
    with open(args.nwk, "r") as in_nwk:
        tree = Phylo.read(in_nwk, "newick")
    with open(args.outgroup_id, "r") as in_outgroup:
        non_or = in_outgroup.read().splitlines()

    outgroup = []
    orgroup = []

# 3. Add clades that contain outgroups to outgroup.
    for non_or_id in non_or:
        outgroup.append(next(tree.find_clades(non_or_id)))

    # 4. The program diverges based on whether the classi_id parameter was
    #    supplied.
    if not args.classi_id:
        # 6.a if outgroup is not monophyletic, then a non or is in the mrca
        #   clade. once thats discarded, will others be revealed?
        if not tree.is_monophyletic(outgroup):
            omrca = tree.common_ancestor(non_or)
            for terminal in omrca.get_terminals():
                if terminal not in orgroup:
                    print(terminal.name)
    else:
        # 5.b Convert class1_id to a list
        with open(args.classi_id, "r") as in_classi_id:
            c1 = in_classi_id.read().splitlines()

        # 6.b Add clades with class I ORs to classI
        classI = []

        for c1_id in c1:
            classI.append(next(tree.find_clades(c1_id)))

        tree.root_with_outgroup(outgroup)

        if not tree.is_monophyletic(classI):
            cImrca = tree.common_ancestor(classI)
            for terminal in cImrca.get_terminals():
                if terminal not in classI:
                    print(terminal.name)


# 1. Take in arguments nwk and outgroup_id, and send them to phylo_tree
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=INFO)
    parser.add_argument("nwk", help=("phylogenetic tree (nwk)"))
    parser.add_argument("outgroup_id", help=("ids of outgroup (txt)"))
    parser.add_argument("-C", "--classi_id", help=("ids of class I (txt)"))
    args = parser.parse_args()
    phylo_tree(args)
