#!/usr/bin/env bash

set -euo pipefail

input_file=$1
prefix="${input_file%%.*}"
gfa_file="${prefix}.gfa"
og_file="${prefix}.og"
paths="${prefix}.haplotypes_nodes.tsv"


vg construct -M "$input_file" -m 1000 | vg view - > "$gfa_file" && \
    odgi build -g "$gfa_file" -o "$og_file" && \
    odgi paths -i "$og_file" -H > "$paths"




# && odgi viz -i kp.og -o kp_1D.png

