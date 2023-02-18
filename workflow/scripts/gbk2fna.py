#!/usr/bin/env python3

import sys
from itertools import chain
from Bio import SeqIO

def main() -> None:
    try:
        with sys.stdin as f:
            for rec in SeqIO.parse(f, 'genbank'):
                for fet in rec.features:
                    if fet.type == 'CDS' and fet.qualifiers.get('locus_tag', None):
                        header = fet.qualifiers['locus_tag'][0]
                        sequence = fet.extract(rec).seq
                        print(f'>{header}\n{sequence}')
    except (KeyboardInterrupt, BrokenPipeError) as e:
        sys.exit(0)

if __name__ == '__main__':
    main()
