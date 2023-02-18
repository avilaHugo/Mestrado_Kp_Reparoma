#!/usr/bin/env python3
"""
List fasta aligmets (builded with mafft), converts to qiime2 artifact,
run aligment mask and run phylogeny

Usage:
    ./convert_mask_and_run_phylogeny.py ( --alignment_fastas_dir=PATH ) ( --output_dir=PATH )
                                        [ --alignment_fastas_suffix=STR ] [ --threads=INT ] [ --use_fasttree ]
                                        [ --bootstrap=INT ] [ --upper_case_aln ] [ --sample_n_genes=INT ]

Options:
    --alignment_fastas_dir=PATH    Input dir with MAFFT aligments to mask.
    --output_dir=PATH              Dir to store output files.
    --alignment_fastas_suffix=STR  Suffix of fasta aligmnets (for find ) [default: .fasta]
    --bootstrap=INT                Number of bootstrap replicats.
    --use_fasttree                 Use fasttree to generate trees [default: raxml]
    --threads=INT                  Number of threads to run phylogeny [default: 1]
    --upper_case_aln               Enchore fasta aln is uppercase: 'acgt' -> 'ACGT'.
    --sample_n_genes=INT         Only run n trees.
"""

# native modules
import os
import sys
import subprocess
from multiprocessing import Pool
from pathlib import Path
from functools import partial
import random
import logging
logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))

# 3rd party modules
from docopt import docopt

def convert_alignment_to_qiime2_format(input_file, output_dir, alignment_fastas_suffix):
    output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace(alignment_fastas_suffix, '.qza'))

    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    logger.info('Coverting {} to qiime2 format.'.format(input_file))
    cmd = [
        'qiime',
        'tools',
        'import',
        '--input-path', input_file ,
        '--output-path', output_file_name,
        '--type', "FeatureData[AlignedSequence]"
    ]

    subprocess.run(
        cmd,
        check=True,
    )
    logger.info('Finished Covertion: {}.'.format(input_file))

    return output_file_name

def mask_alignments(input_file, output_dir):
    output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace('.qza', '.masked.qza'))

    logger.info('Masking: {}.'.format(input_file))
    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    cmd = [
        'qiime',
        'alignment',
        'mask',
        '--i-alignment', input_file,
        '--o-masked-alignment', output_file_name,
    ]

    subprocess.run(
        cmd,
        check=True,
    )
    logger.info('Finished Mask: {}.'.format(input_file))

    return output_file_name

def run_raxml(input_file, output_dir, threads, bootstrap=None):
    if bootstrap:
        logger.info('Bootstrap is active: {} replicates.'.format(bootstrap))
        output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace('.masked.qza', '.GTRCAT.BS.tree.qza'))
        cmd = [
            'qiime',
            'phylogeny',
            'raxml-rapid-bootstrap',
            '--i-alignment', input_file,
            '--p-seed', '1723',
            '--p-rapid-bootstrap-seed', '9384',
            '--p-bootstrap-replicates', str(bootstrap),
            '--p-substitution-model', 'GTRCAT',
            '--o-tree', output_file_name,
            '--p-n-threads', str(threads),
        ]
    else:
        logger.info('Bootstrap is deactivated.'.format(bootstrap))
        output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace('.masked.qza', '.GTRCAT.tree.qza'))
        cmd = [
            'qiime',
            'phylogeny',
            'raxml',
            '--p-seed', '1723',
            '--i-alignment', input_file,
            '--p-substitution-model', 'GTRCAT',
            '--o-tree', output_file_name,
            '--p-n-threads', str(threads),
        ]

    logger.info('Running tree for: {}.'.format(input_file))
    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    subprocess.run(
        cmd,
        check=True,
    )

    logger.info('Finished tree: {}.'.format(input_file))

    return output_file_name

def extract_tree_from_qiime2(input_file, output_dir):
    output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace('.GTRCAT.tree.qza', '.GTRCAT.tree'))

    logger.info('Extracting tree from: {}.'.format(input_file))
    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    dir_output = os.path.join(output_dir, os.path.basename(input_file).replace('.GTRCAT.tree.qza', ''))

    cmd = [
        'qiime',
        'tools',
        'export',
        '--input-path', input_file,
        '--output-path', dir_output,
    ]

    subprocess.run(
        cmd,
        check=True,
    )

    extracted_tree = os.path.join(dir_output, 'tree.nwk')

    # Moving extracted file
    Path(extracted_tree).rename(output_file_name)

    # Removing empty dir
    os.rmdir(dir_output)

    logger.info('Finished tree extration: {}.'.format(input_file))

    return output_file_name

def uppercase_aln_fasta(input_file, output_dir):
    output_file_name = os.path.join(output_dir, os.path.basename(input_file))

    logger.info('Uppercase file: {}.'.format(input_file))
    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    cmd = "awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {} > {}".format(input_file, output_file_name)
    subprocess.run(
        cmd,
        check=True,
        shell=True,
    )
    logger.info('Finished Uppercase file: {}.'.format(input_file))
    return output_file_name

def run_fasttree(input_file, output_dir, threads):
    output_file_name = os.path.join(output_dir, os.path.basename(input_file).replace('.masked.qza', '.GTRCAT.tree.qza'))

    logger.info('Running tree for: {}.'.format(input_file))
    if os.path.exists(output_file_name):
        logger.info('File alredy exits, skiping creation: {}'.format(output_file_name))
        return output_file_name

    cmd = [
        'qiime',
        'phylogeny',
        'fasttree',
        '--i-alignment', input_file,
        '--o-tree', output_file_name,
        '--p-n-threads', str(threads),
    ]

    subprocess.run(
        cmd,
        check=True,
    )

    logger.info('Tree finished: {}.'.format(input_file))

    return output_file_name

def main(*args, **kwargs):

    # Logging setup
    logging.basicConfig(
        level=logging.DEBUG,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )
    logger.info('ARGS: {}'.format(kwargs))

    threads = int(kwargs['--threads'])
    assert threads >= 1, 'Threads need to be more equal than 1'
    raxml_threads = 4
    if threads < raxml_threads:
        raxml_threads = 1
    raxml_instances = max(1, int(threads/raxml_threads))
    logger.info('Running with {} threads.'.format(threads))
    logger.info('Running with {} raxml_instances.'.format(raxml_instances))
    logger.info('Running with {} raxml_threads.'.format(raxml_threads))

    kwargs['--alignment_fastas_dir'] = os.path.abspath(kwargs['--alignment_fastas_dir'])
    kwargs['--output_dir'] = os.path.abspath(kwargs['--output_dir'])

    assert os.path.exists(kwargs['--alignment_fastas_dir'])
    if not os.path.exists(kwargs['--output_dir']):
        logger.warning('Output dir already exists, maybe some files will be overwriten !!!')
    Path(kwargs['--output_dir']).mkdir(parents=True, exist_ok=True)

    logger.info('Listing fasta aln files from {} ...'.format(kwargs['--alignment_fastas_dir']))

    find_list = subprocess.run(
        [
            'find',
            kwargs['--alignment_fastas_dir'],
            '-mindepth', '1',
            '-maxdepth', '1',
            '-type', 'f',
            '-name', '*' + kwargs['--alignment_fastas_suffix'],
        ],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout.decode('utf-8').splitlines()

    if kwargs['--sample_n_genes']:
        logger.info('SAMPLING GENES IS ACTIVE, only {} will be runned !!'.format(kwargs['--sample_n_genes']))
        find_list = random.sample(find_list, int(kwargs['--sample_n_genes']))
        logger.info('Length of input fastas after sampling: {}'.format(len(find_list)))

    if kwargs['--upper_case_aln']:
        logger.info('MAKING SHORE ALIGNMENTS ARE UPPERCASED')
        with Pool(processes=threads) as p:
            find_list = p.map(
                partial(
                    uppercase_aln_fasta,
                    output_dir=kwargs['--output_dir'],
                ),
                find_list
            )

    assert find_list, 'No files found at {} with suffix {}'.format(kwargs['--alignment_fastas_dir'], kwargs['--alignment_fastas_suffix'])

    logger.info('Converting alignments to qiime2 format ...')
    with Pool(processes=threads) as p:
        converted_files = p.map(
            partial(
                convert_alignment_to_qiime2_format,
                output_dir=kwargs['--output_dir'],
                alignment_fastas_suffix=kwargs['--alignment_fastas_suffix'],
            ),
            find_list,
        )
    logger.info('Finished alignments convertion.')

    logger.info('Masking alignments ...')
    with Pool(processes=threads) as p:
        masked_files = p.map(
            partial(
                mask_alignments,
                output_dir=kwargs['--output_dir'],
            ),
            converted_files,
        )


    if kwargs['--bootstrap']:
        kwargs['--bootstrap'] = int(kwargs['--bootstrap'])
        assert kwargs['--bootstrap'] >= 1, 'Invalid bootstrap value.'
        logger.info("Bootstrap is active: {} replicates".format(kwargs['--bootstrap']))


    if not kwargs['--use_fasttree'] == 'raxml':
        logger.info('Using fasttree with {} threads and {} instances'.format(raxml_threads, raxml_instances))
        with Pool(processes=raxml_instances) as p:
            qiime2_trees = p.map(
                partial(
                    run_fasttree,
                    output_dir=kwargs['--output_dir'],
                    threads=raxml_threads,
                ),
                masked_files,
            )
    else:
        with Pool(processes=raxml_instances) as p:
            qiime2_trees = p.map(
                partial(
                    run_raxml,
                    output_dir=kwargs['--output_dir'],
                    threads=raxml_threads,
                    bootstrap=kwargs['--bootstrap'],
                ),
                masked_files,
            )


    with Pool(processes=threads) as p:
        _ = p.map(
            partial(
                extract_tree_from_qiime2,
                output_dir=kwargs['--output_dir'],
            ),
            qiime2_trees,
        )

    logger.info('FINISHED ALL JOBS !')

if __name__ == '__main__':
    main(**docopt(__doc__))
