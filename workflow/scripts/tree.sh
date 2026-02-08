#!/usr/bin/env bash

set -euo pipefail

script_path='/home/hugo.avila/Projects/reparoma/workflow/scripts/get_tree.sh'



awk -F ':' '$2 >= 500 {print $1}' REDUCED  | \
    sed 's:.faa:.aln:;s:^:'$PWD'/results/checkm_tree/fastas/:' | \
    parallel \
    -j 5 \
    'echo "STARTED: {/.}" && \
    cd {=s:fastas.+::;s:$:trees:=}/{/.} && \
    time bash '"${script_path}"' {/.}.aln 10 &> {/.}.main.log && \
    echo "DONE: {/.}"'

# # ORGINAL
#     'echo "STARTED: {/.}" && \
#     mkdir -p {=s:fastas.+::;s:$:trees:=}/{/.} && \
#     cd {=s:fastas.+::;s:$:trees:=}/{/.} && \
#     cp {} {=s:fastas.+::;s:$:trees:=}/{/.} && \
#     time bash '"${script_path}"' {/.}.aln &> {/.}.main.log && echo "DONE: {/.}"'
