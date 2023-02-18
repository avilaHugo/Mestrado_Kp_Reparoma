#!/usr/bin/env bash

set -euo pipefail

input_fasta="$1"
output_aln="$2"
number_of_threads="${3-1}"

exec 2> "${input_fasta//.faa/.log}"

mafft \
    --thread $number_of_threads \
    --auto \
    $input_fasta > $output_aln 
