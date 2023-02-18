#!/usr/bin/env bash


WASTRID_BIN=/home/hugo.avila/Projects/reparoma/submodules/wastrid/wastrid
INPUT_NAME=$1
MYTHREADS=${2:-2}
OUTPUT_NAME="${INPUT_NAME//.nwk/.wastrid.nwk}"
LOG_NAME="${OUTPUT_NAME//.nwk/.log}"


time $WASTRID_BIN \
    -t $MYTHREADS \
    -i $INPUT_NAME \
    --preset vanilla \
    -o $OUTPUT_NAME |& tee $LOG_NAME
