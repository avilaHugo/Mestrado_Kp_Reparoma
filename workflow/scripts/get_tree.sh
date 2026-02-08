#!/usr/bin/env bash

set -euo pipefail

ALN_NAME="$1"
PHY_CORES="${2:-2}"

. ~/miniconda3/etc/profile.d/conda.sh

ALN_NAME="${ALN_NAME}.reduced"

# STEP: Getting a starter tree
conda activate fasttree_env
export OMP_NUM_THREADS=$PHY_CORES
FastTreeMP -lg -gamma -seed 12345 \
    "${ALN_NAME}" > "fasttree.nwk"
unset OMP_NUM_THREADS
conda deactivate

 
# STEP: remove outliers 
conda activate treeshrink_env
run_treeshrink.py -t "fasttree.nwk" -o tree_filt
conda deactivate

[[ "$(tr '\t' '\n' < tree_filt/output.txt | grep -v '^$' | wc -l)" > 0 ]] && {

    echo MAIOR $ALN_NAME

    conda activate seqkit_env

    ALN_NAME_SHRINKED="${ALN_NAME//.reduced/}.shrinked"

    seqkit \
        grep -nf <(sed '1d' "${ALN_NAME}" | awk '{print $1}' | grep -vFf <(tr '\t' '\n' < tree_filt/output.txt | grep -v '^$') -) \
        "${ALN_NAME//.reduced/}" \
        > "${ALN_NAME_SHRINKED}" && \
        ALN_NAME="${ALN_NAME_SHRINKED}"

    conda deactivate

}

# STEP: Getting best tree
conda activate raxml_env
raxmlHPC-PTHREADS-AVX2 -m PROTCATLG -F -f D -D -s "${ALN_NAME}" -T "${PHY_CORES}" -p 42    -n FTCAT -t "tree_filt/output.nwk" && echo FTCAT
conda deactivate

# STEP: OPTMIZE tree branches with GAMMA
conda activate iqtree_env
ls RAxML_result.* | \
    while IFS=. read prefix suffix;do
        [ "${suffix}" == FTCAT ] || continue # for my dataset running only for the FTCAT usually worked better

        iqtree \
            -m LG+G4 \
            -s "${ALN_NAME}" \
            -te "${prefix}.${suffix}" \
            --prefix "${suffix}" \
            -T "${PHY_CORES}" && echo "DONE scoring: ${suffix}"

    done
conda deactivate
