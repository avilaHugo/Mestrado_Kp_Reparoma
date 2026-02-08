#!/usr/bin/env bash

set -e

REPAROMA_DB=/home/hugo.avila/Projects/reparoma/data/reparoma/DB
BIN=/home/hugo.avila/Projects/reparoma/submodules/pseudofinder/pseudofinder.py
OUTPUT_PREFIX=/home/hugo.avila/Projects/reparoma/results/pseudofinder
INPUT_FILE=$1
RUN_THREADS=${2:-4}


echo python3 $BIN \
    annotate \
    --genome ${INPUT_FILE} \
    --outprefix "$(mkdir -pv ${OUTPUT_PREFIX}/$(basename $INPUT_FILE .gbk) |& awk '{print $4}' | tr -d \' )/$(basename $INPUT_FILE .gbk)" \
    --database $REPAROMA_DB \
    --threads $RUN_THREADS \
    --diamond
