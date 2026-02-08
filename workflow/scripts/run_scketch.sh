#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR=/home/hugo.avila/Projects/reparoma/results/mash_dist/bootstrap

grep -v '#' MASH_SEEDS.txt | while read seed;do
    SEED_DIR="${ROOT_DIR}/${seed}"
    echo
    echo "STARTING: ${seed}"

    time cat FILTERED_kleborate.txt | \
    parallel -j 10 'mash sketch -o '$SEED_DIR'/{/} -S '$seed' {}' && \
    mash paste "${ROOT_DIR}/seed_${seed}.skts" -l <(ls $SEED_DIR/*fna.msh) && \
    echo "FINISHED PASTING: seed ${seed}" && \
    find $SEED_DIR -name '*fna.msh' -delete
    
    echo "FINISHED: ${seed}"
    echo
done 
