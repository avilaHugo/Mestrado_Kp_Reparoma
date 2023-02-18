#!/usr/bin/env bash 


set -e

# snakemake -pc 60 --use-conda $@

# find -L results/prokka/ \
    # -type f \
    # -name '*.gbk.gz' | \
    # parallel -j 150 'any2fasta -nu {} 1>data/{=s:.+/::;s:.gz$::=} 2>data/{=s:.+/::;s:.gz$::=}.log' 

# poppunk_visualise --ref-db /home/hugo.avila/Projects/reparoma/results/poppunk/s_genus_poppunk_db --model-dir /home/hugo.avila/Projects/reparoma/results/poppunk/fit_models/dbscan --output /home/hugo.avila/Projects/reparoma/results/poppunk/viz/dbscan --microreact --threads 40


# CheckM

# checkm taxonomy_wf -t 64 \
#     --nt \
#     -x fna \
#     species 'Klebsiella oxytoca' \
#     /home/hugo.avila/Projects/reparoma/data/Klebsiella_oxytoca \
#     /home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_oxytoca

# checkm taxonomy_wf  -t 60 \
#     --nt \
#     -x fna \
#     genus 'Klebsiella' \
#     /home/hugo.avila/Projects/reparoma/data/Klebsiella_genus \
#     /home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_genus

# SPECIE='/home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_oxytoca'
# checkm qa \
#     -o 9 \
#     -f "${SPECIE}/Klebsiella_oxytoca.faa" \
#     -t 40 \
#     "${SPECIE}/"'Klebsiella oxytoca.ms' \
#     "${SPECIE}"

# SPECIE='/home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_genus'
# checkm qa \
#     -o 9 \
#     -f "${SPECIE}/Klebsiella.faa" \
#     -t 40 \
#     "${SPECIE}/"'Klebsiella.ms' \
#     "${SPECIE}"


# SPECIE='/home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_pneumoniae'
# checkm qa \
#     -o 9 \
#     -f "${SPECIE}/Klebsiella_pneumoniae.faa" \
#     -t 40 \
#     "${SPECIE}/"'Klebsiella pneumoniae.ms' \
#     "${SPECIE}"

# for i in {6..8};do
#     SPECIE='/home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_oxytoca'
#     FMT=$i
#     checkm qa \
#         -o $FMT \
#         -f "${SPECIE}/marker_genes_fmt_${FMT}.tsv" \
#         -t 40 \
#         --tab_table \
#         "${SPECIE}/"'Klebsiella oxytoca.ms' \
#         "${SPECIE}"
# done

# SPECIE='/home/hugo.avila/Projects/reparoma/results/checkm/Klebsiella_pneumoniae'
# FMT=4
# checkm qa \
#     -o $FMT \
#     -f "${SPECIE}/marker_genes_fmt_${FMT}.tsv" \
#     -t 40 \
#     --tab_table \
#     "${SPECIE}/"'Klebsiella pneumoniae.ms' \
#     "${SPECIE}"



# MAFFT

# Create a single script file to use this generic line
mafft \
    --thread $MAFFT_CORES \
    --quiet --auto \
    fasta_file > output_file

