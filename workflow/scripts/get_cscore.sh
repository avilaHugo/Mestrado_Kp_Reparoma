#!/usr/bin/env bash

msa_script_path=$1
aln_file=$2

python3 \
    "$msa_script_path" \
    "$aln_file" | \
    awk -v OFS=',' '{print $1, $3}' > "${aln_file//.aln/}.cscore.csv"
