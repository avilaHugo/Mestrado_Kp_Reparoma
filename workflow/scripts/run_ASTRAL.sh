#!/usr/bin/env bash


ASTRAL_BIN=/home/hugo.avila/Projects/reparoma/submodules/ASTER/bin/astral-hybrid
GET_LABELS_BIN=/home/hugo.avila/Projects/reparoma/workflow/scripts/extract_labels.py
# MAPPING_FILE=$(mktemp)
ASTRAL_CORES=40
ASTRAL_RAM=40g
INPUT_FILE=$1
OUTPUT_FILE="${INPUT_FILE//.nwk/.astral.nwk}"
LOG_FILE="${INPUT_FILE//.nwk/.astral.log}"
KLEBORATE='/home/hugo.avila/Projects/reparoma/data/kleborate_and_checkm_filtered_genomes.tsv'


# python3 $GET_LABELS_BIN $INPUT_FILE | \
#     sed 's:$:,:' | \
#     grep -Ff - $KLEBORATE | \
#     sed 's: :_:;s:,: :' > $MAPPING_FILE

# awk -F',' '!a[$2]++ {print $2}' $LABELS_TMP | \
#     while read specie;do
#         printf "${specie}: " && grep ",${specie}" $LABELS_TMP | \
#         cut -f 1 -d',' | \
#         tr '\n' ',' | \
#         sed 's:,$::' && echo
#     done > $MAPPING_FILE

time $ASTRAL_BIN \
    -t $ASTRAL_CORES \
    -i $INPUT_FILE \
    -o $OUTPUT_FILE |& tee $LOG_FILE

rm -v $MAPPING_FILE
