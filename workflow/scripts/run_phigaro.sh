#!/usr/bin/env bash

INPUT_FILE=$1
MYTHREADS=${2:-2}
SAMPLE_PREFIX=$( basename ${INPUT_FILE} .fna )
OUTPUT_DIR="/home/hugo.avila/Projects/reparoma/results/phigaro/${SAMPLE_PREFIX}"

# mkdir -p $OUTPUT_DIR

echo "STARTED: $SAMPLE_PREFIX"

time phigaro \
    --fasta-file "${INPUT_FILE}" \
    --extension tsv \
    --output "${OUTPUT_DIR}" \
    --threads $MYTHREADS \
    --delete-shorts &> "$( dirname ${OUTPUT_DIR} )/${SAMPLE_PREFIX}.run.log"

echo "FINISHED: $SAMPLE_PREFIX"
