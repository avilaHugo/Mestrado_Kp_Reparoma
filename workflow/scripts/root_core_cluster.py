#!/usr/bin/env python3
"""
Generate multifastas to align.

Usage:
    root_core_cluster.py ( --panaroo_genus_dir=PATH ) ( --specie=STR ) ( --roots_file=PATH )
                         ( --annotation_dir=PATH ) ( --panaroo_dir=PATH ) ( --fasta_suffix=STR )
                         ( --output_dir=PATH ) ( --core_clusters_file=PATH ) [ --threads=N ]

Options:
    -h --help                     Show this screen.
    -s --specie=STR               Specie to use (remember to quote the string).
    -r --roots_file=PATH          Path to the roots file (tsv from mmseqs).
    -a --annotation_dir=PATH      Path to the annotation directory (main prokka annotation dir).
    -p --panaroo_dir=PATH         Path to the panaroo output directory.
    -f --fasta_suffix=STR         Suffix of the prefixed fasta files: 'SPECIE#ALLEL#CLUSTER'.
    -c --core_clusters_file=PATH  Path to the core genes file (one per line).
    -o --output_dir=PATH          Path to the output directory.
    -t --threads=N                Threads number [default: 1].
"""

# Native imports
import os
import sys
import datetime
import datetime
import re
from pathlib import Path
from collections import defaultdict
import uuid
import logging
import subprocess
from multiprocessing import Pool
import itertools
logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))

# 3th party modules
from docopt import docopt

def read_roots(roots_file: str, separator="\t") -> list:
    """
    Read a file with roots for a given species.
    """
    root_clusters = defaultdict(list)
    with open(roots_file, 'r') as roots_file:
        for line in roots_file.readlines():
            pan_cluster_name, pan_sequence = line.strip().replace('"', '').split(separator)
            root_clusters[pan_cluster_name].append(pan_sequence)

    result = {}
    for root_cluster in root_clusters:
        names = defaultdict(set)
        for root_sequence in root_clusters[root_cluster]:
            names[root_sequence.split('#')[0]].add(root_sequence)
        result[root_cluster] = dict(names)

    return result



def list_and_filter_files(annotation_dir: str, genome_ids: set, fasta_suffix: str) -> list:
    """
    List all the files in the annotation directory and filter them by specie and fasta_suffix.
    """
    genomes_fasta = {}
    for fasta in Path(annotation_dir).rglob('*' + fasta_suffix):
        genome_id = os.path.basename(fasta).replace(fasta_suffix, "")
        if genome_id in genome_ids:
            assert not genome_id in genomes_fasta, "Duplicated genome id."
            genomes_fasta[genome_id] = str(fasta.absolute())
    return genomes_fasta

def read_core_clusters_from_file(core_clusters_file: str) -> dict:
    with open(core_clusters_file, 'r') as f:
        core_clusters = set(line.strip() for line in f.readlines())
    return core_clusters

def extract_species_genome_ids_from_panaroo_out(panaroo_dir: str, panaroo_table: str = 'gene_presence_absence.csv', separator: str = ',' ) -> set:
    """
    Extract the species genomes ids from the panaroo output.
    """
    panaroo_out = os.path.join(panaroo_dir, panaroo_table)
    assert os.path.exists(panaroo_out), "Panaroo output does not exist: {}.".format(panaroo_out)
    with open(panaroo_out, 'r') as f:
        genome_ids = set(f.readlines()[0].strip().split(separator)[3:])
    return genome_ids

def list_roots_fasta_files(panaroo_genus_dir: str, panaroo_representatives: str = 'pan_genome_reference.prefixed.fa.gz') -> None:
    """
    List the roots fasta files.
    """
    roots_fastas = {}
    for panaroo_specie in filter(lambda x: x.is_dir(), Path(panaroo_genus_dir).iterdir()):
        for fasta in panaroo_specie.rglob('*' + panaroo_representatives ):
            specie = os.path.basename(os.path.dirname(fasta))
            fasta_abspath = str(fasta.absolute())
            assert not specie in roots_fastas, "Duplicated species."
            assert os.path.exists(fasta_abspath), "FASTA FILE NOT FOUN"
            roots_fastas[specie] = fasta_abspath
    return roots_fastas

def extract_record_from_fasta_faidx(fasta_file: str, record_id: str) -> str:
    """
    Extract a record from a samtools indexed fasta file.
    """
    try:
        record = subprocess.check_output(['samtools', 'faidx', fasta_file, record_id], encoding='utf-8').split('\n')
        header = record[0]
        sequence = ''.join(record[1:-1])
        return f'{header}\n{sequence}'
    except subprocess.CalledProcessError as e:
        logger.info("Extracting record_id ('{}') from fasta file, RECORD might not exists on fasta: {}.".format(record_id, e))

def main(panaroo_genus_dir: str, core_clusters_file: str, specie: str, roots_file: str, annotation_dir: str, panaroo_dir: str, fasta_suffix: str, output_dir: str, threads: int = 1, *args, **kwargs) -> None:

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )

    assert all([os.path.exists(x) for x in [roots_file, annotation_dir, panaroo_dir]]), "One or more paths do not exist."

    logger.info("Extracting species genome_ids from panaroo output ...")
    genome_ids = extract_species_genome_ids_from_panaroo_out(panaroo_dir=panaroo_dir)

    logger.info("Reading species core genes ...")
    core_clusters = read_core_clusters_from_file(core_clusters_file=core_clusters_file)

    logger.info("Listing genomes fastas ...")
    genomes_fastas = list_and_filter_files(annotation_dir=annotation_dir, genome_ids=genome_ids, fasta_suffix=fasta_suffix)

    logger.info("Reading roots ...")
    roots = read_roots(roots_file=roots_file)

    logger.info("Listing roots fasta files ...")
    roots_fastas = list_roots_fasta_files(panaroo_genus_dir=panaroo_genus_dir)

    logger.info("Generating multifastas ...")

    fasta_ids = set()
    fasta_renamed = []

    assert roots_fastas, "Empty dict"

    for core_cluster in core_clusters:
        logger.info("Generating multifasta for core cluster: {}".format(core_cluster))

        # check if id is in fasta_ids, if it is, create another
        cluster_fasta_id = uuid.uuid1().hex
        while cluster_fasta_id in fasta_ids:
            cluster_fasta_id = uuid.uuid1().hex
        fasta_ids.add(cluster_fasta_id)

        cluster_fasta_file = os.path.join(output_dir, cluster_fasta_id + ".fasta")
        logger.info("Writing multifasta to: {}".format(cluster_fasta_file))

        with open(cluster_fasta_file, 'w') as f_out:
            with Pool(processes=threads) as pool:
                logger.info("Initiating pool of {} threads for core_cluster: {}".format(threads, core_cluster))
                specie_sequences = pool.starmap(
                    extract_record_from_fasta_faidx,
                    [(genomes_fastas[genome_id] , f'{genome_id}#{core_cluster}') for genome_id in genome_ids if genomes_fastas.get(genome_id) ]
                )
                specie_sequences = tuple(filter(None, specie_sequences))
                logger.info("Finished pool of {} threads for core_cluster: {}".format(threads, core_cluster))

                # Root sequences
                root_sequences = []

                logger.info("Rooting cluster {} ...".format(core_cluster))
                cluster_roots = roots.get(f'{specie}#{core_cluster}', [])
                if cluster_roots or len(cluster_roots) > 1:
                    cluster_roots = {k : v for k,v in cluster_roots.items() if k != specie}
                    logger.info("Cluster {} has {} roots.".format(core_cluster, sum([len(v) for v in cluster_roots.values()])))

                    root_sequences = pool.starmap(
                        extract_record_from_fasta_faidx,
                        [ (roots_fastas[root_specie], root_id) for root_specie, root_ids in cluster_roots.items() for root_id in root_ids ]
                        )
                    root_sequences = tuple(filter(None, root_sequences))
                    fasta_renamed.append(f'{core_cluster},{cluster_fasta_file},rooted')
                else:
                    fasta_renamed.append(f'{core_cluster},{cluster_fasta_file},unrooted')
                    logger.info("Cluster {} has no roots.".format(core_cluster))

                # merge specie and root sequences and write to file
                logger.info("Merging specie (count: {}) and root (count: {}) sequences for cluster {}".format(len(specie_sequences), len(root_sequences), core_cluster))
                merged_sequences = iter(itertools.chain(specie_sequences, root_sequences))

                # write to file
                logger.info("Writing merged sequences to file {} for cluster {}".format(cluster_fasta_file, core_cluster))
                f_out.write('\n'.join(merged_sequences))


    summary_name = os.path.abspath(os.path.join(output_dir, 'renamed_fasta_names.csv'))
    logger.info("Writing summary file: {} ...".format(summary_name))

    with open(summary_name, 'w') as f:
        f.write("\n".join(fasta_renamed))

    logger.info("Finished.")

if __name__ == "__main__":
    clean_args = lambda args: { k.replace('-', '') : v for k, v in args.items() }
    args = clean_args(docopt(__doc__))
    args['threads'] = int(args['threads'])
    assert args['threads'] > 0, "Threads must be greater than 0."
    main(**args)
