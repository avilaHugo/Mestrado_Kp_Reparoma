#!/usr/bin/env python3
import sys
from itertools import chain
from Bio import Phylo


def main(fn):
     
    nodes = chain.from_iterable(map(lambda tree: tree.find_clades(), Phylo.parse(fn, "newick")))

    node_names = filter(bool, map(lambda node: node.name, nodes))

    print(*set(node_names), sep='\n')


if __name__ == '__main__':
    main(*sys.argv[1:]) 
