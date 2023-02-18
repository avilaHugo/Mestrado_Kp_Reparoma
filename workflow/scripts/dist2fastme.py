#!/usr/bin/env python3

import sys
import random
import string 


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def read_mask(fn):
    with open(fn, 'r') as f:
        return dict(map(lambda line: line.strip().split('\t'), f))

def create_mask(fn):
    with open(fn, 'r') as f:
        genome_ids = {}
        for genome_id in next(f).strip().split('\t')[1:]:
            while (new_id := id_generator(20)) in genome_ids:
                continue
            genome_ids[new_id] = genome_id
    return { v : k for k, v in genome_ids.items() }

def save_mask(fn, created_mask):
    with open(fn, 'w') as f:
        f.write('\n'.join(map('\t'.join, created_mask.items())) + '\n')

def main(mash_dist, mask_file):
    try:
        mask = read_mask(mask_file)
    except FileNotFoundError:
        print(f'Mask not fount, creating new mask at {mask_file} ...', file=sys.stderr)
        mask = create_mask(mash_dist)
        save_mask(mask_file, mask)

    print(len(mask))
    with open(mash_dist, 'r') as f:
        _ = next(f)
        for lines in map(str.strip, f):
            genome_id, *dists = lines.split('\t')
            print(
                mask[genome_id],
                *map(lambda dist: format(float(dist), '.10f'), dists),
                sep='\t'
            )

            

        

        
        

    

if __name__ == '__main__':
    main(*sys.argv[1:])
