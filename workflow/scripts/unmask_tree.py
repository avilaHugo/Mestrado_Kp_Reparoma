#!/usr/bin/env bash


import sys
from Bio import Phylo
from Bio.Phylo.NewickIO import parse

def main(mask, tree):
    # read in the tree and parse it into a tree object
    with open(tree) as nj:
        tree = next(parse(nj))

    # read in the mask file and create a dictionary mapping the original names to the new names
    with open(mask) as f:
        mask = dict(map(lambda line: line.strip().split('\t'), f))

    # iterate through the tree and replace the names of the nodes with the new names from the mask
    for node in tree.get_terminals():
        node.name = mask[node.name]

    # write the modified tree to a new file
    # Phylo.write(tree, "modified_tree.nwk", "newick")
    print(tree.format('newick'))



if __name__ == '__main__':
    main(*sys.argv[1:])
