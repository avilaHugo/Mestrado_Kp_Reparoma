#!/usr/bin/env python3
"""
Reads a panaroo csv outputs and gbk files to create concatenated genes reference genomes

Usage:
    get_pangenome_genes.py -h
    get_pangenome_genes.py ( --panaroo_results_preffix=PATH ) ( --coregenome_threshold=FLOAT )
                           ( --annotation_path_list=PATH ) [ --core_cluster_names_to_file=PATH  ]
                           [ --threads=N ] [ --debug ]

Options:
    -i --panaroo_results_preffix=PATH       Dir of the panaroo result dir with prefix of the results tables to use.
    -c --coregenome_threshold=FLOAT         Min percentange (float > 0 and <= 1 ) Ex: 0.60,0.84,0.95, 1.
    -a --annotation_path_list=PATH          Path to the directorie where the prokka are stored.
    -p --core_cluster_names_to_file=PATH    Write the names of coregenome clusters.
    -t --threads=N                          Write the names of coregenome clusters.
    -d --debug                              Run in debug mode (more verbose log).
"""

import os
import sys
import gzip
import logging
import subprocess
from multiprocessing import Pool
from pathlib import Path
from collections import namedtuple
from Bio import SeqIO
from docopt import docopt

logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))

def get_min_persistence(coregenome_threshold: float, genome_count: int) -> int :
    assert 0 < coregenome_threshold <= 1, "'coregenome_threshold' needs to be greater than 0 or less equal than 1."
    return int(genome_count * (1 - coregenome_threshold))

def get_cluster_names(panaroo_Rtab: str, coregenome_threshold: float) -> tuple:
    with open(panaroo_Rtab, 'r') as f:
        file_lines = iter(map(lambda x: x.strip().split('\t'), f.readlines()))
        header = next(file_lines)

        genome_count = len(header) - 1
        logger.info(f'Genome count: {genome_count}')

        min_persistence = get_min_persistence(coregenome_threshold=coregenome_threshold, genome_count=genome_count)
        logger.info(f'Min persistence: {min_persistence}')

        logger.info('Reading panaroo Rtab ...')
        return set(map(lambda x: x[0], filter(lambda x: x.count('0') < min_persistence, file_lines)))

def get_coregenome_data(panaroo_csv: str, cluster_names: set):
    logger.info(f'Extracting locustags from: {panaroo_csv}')
    data = {}
    with open(panaroo_csv, 'r') as f:
        file_lines = iter(map(lambda x: x.strip().split(','), f.readlines()))
        data['columns'] = next(file_lines)
        coregenome_lines = filter(lambda x: x[0] in cluster_names, file_lines)
        for i in coregenome_lines:
            data[i[0]] = i

    clusters_names = [i for i in data.keys() if i != 'columns']
    logger.info('Finished locustags extraction.')

    genome_data = namedtuple('genomeData', 'genome_id locustags')

    for idx, columns in tuple(enumerate(data['columns']))[3:]:
        locustags = dict(map(lambda x: (data[x][idx], x), clusters_names))
        genome_id = columns
        yield genome_data(genome_id, locustags)

def process_genome_list(genome_list_path: str) -> dict:
    with open(genome_list_path, 'r') as f:
        data = {}
        for line in iter(map(str.strip, f.readlines())):
            assert os.path.isfile(line)
            genome_id = os.path.basename(line).replace('.gbk.gz', '')
            assert not data.get(genome_id), f'Genome ID ({genome_id}) is not unique.'
            data[genome_id] = line
    return data

def write_coregenome_cluster_names(clusters_names, output_path):
    with open(output_path, 'w') as f:
        for i in clusters_names:
            f.write(f'{i}\n')

def samtools_compress_and_index(fasta_path: str):
    assert os.path.isfile(fasta_path), f'File does not exist: {fasta_path}'

    logger.debug(f'Compressing : {fasta_path}')
    subprocess.run(['bgzip', fasta_path], check=True)

    logger.debug(f'Indexing : {fasta_path}')
    subprocess.run(['samtools', 'faidx', fasta_path + '.gz'], check=True)

    assert os.path.isfile(fasta_path + '.gz'), f'File does not exist: {fasta_path}.gz'

def fasta_clusters_from_gbk(genome_id, locustags, genome_path, no_clubber = True):
    output_path = os.path.join(os.path.dirname(genome_path), genome_id + '.panclusters.fa')

    if no_clubber and ( Path(output_path + '.gz' ).exists() and Path(output_path + '.gz.fai' ).exists() ):
        logger.info('Skiping file creation: {}'.format(output_path))
        return
        
    logger.info(f'Creating file: {output_path}')
    with open(output_path, 'w') as f_out:
        with gzip.open(genome_path, 'rt') as f:
            for record in SeqIO.parse(f, 'genbank'):
                for feature in filter(lambda x: x.type == 'CDS', record.features):
                    locustag = feature.qualifiers['locus_tag'][0]
                    if locustag in locustags:
                        cluster_name =  locustags[locustag]
                        # f_out.write(f'>{genome_id}#{cluster_name}#{locustag}\n{feature.location.extract(record).seq}\n')
                        f_out.write(f'>{genome_id}#{cluster_name}\n{feature.location.extract(record).seq}\n')
    samtools_compress_and_index(output_path)
    logger.info(f'FINISHED: {output_path}')

def main(panaroo_results_preffix, coregenome_threshold, annotation_path_list, core_cluster_names_to_file, threads=1, debug=False, *args, **kwargs):
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )

    coregenome_threshold = float(coregenome_threshold)

    logger.debug(str(locals()))

    logger.info('Reading genome paths list ...')
    genome_paths = process_genome_list(annotation_path_list)

    cluster_names = get_cluster_names(panaroo_Rtab=panaroo_results_preffix + '.Rtab', coregenome_threshold=coregenome_threshold)

    logger.info(f'Coregenome genes count: {len(cluster_names)}')
    if core_cluster_names_to_file:
        logger.info('Writing core cluster names to file: {}'.format(core_cluster_names_to_file))
        write_coregenome_cluster_names(cluster_names, core_cluster_names_to_file)

    coregenome_data = get_coregenome_data(panaroo_csv=panaroo_results_preffix + '.csv', cluster_names=cluster_names)

    logger.info('Extracting fasta sequences ...')
   
    logger.info('Creating jobs ...') 
    jobs = []
    for genome_data in coregenome_data:
        idx = genome_paths.get(genome_data.genome_id, None)
        if idx:
            jobs.append((genome_data.genome_id, genome_data.locustags, idx))

    logger.info('Starting jobs ...')
    with Pool(processes=int(threads)) as p:
        _ = p.starmap(fasta_clusters_from_gbk, jobs)

    logger.info('All jobs finished.')

if __name__ == '__main__':
    clean_dashes_from_args = lambda x: { i.replace('--', '') : y for i,y in x.items() }
    main(**clean_dashes_from_args(docopt(__doc__)))
