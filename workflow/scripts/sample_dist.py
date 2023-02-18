#!/usr/bin/env python3
import sys
import random

def sample_dist(fn, sample_size):
    with open(fn, 'r') as f:
        lines = map(str.strip, f)
        header = int(next(lines))
        random_recs = set(random.sample(range( header ), sample_size))
        print(sample_size)
        for _, line in filter(lambda x: x[0] in random_recs, enumerate(lines)):
            genome_id, *dists = line.split('\t')
            print(genome_id, *(dist for i, dist in enumerate(dists) if i in random_recs), sep='\t')


if __name__ == '__main__':
    fn, sample_size = sys.argv[1:]
    sample_size = int(sample_size)
    sample_dist(fn, sample_size)
