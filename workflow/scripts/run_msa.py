#!/usr/bin/env python3
"""
Runs the MSA algorithm on a list of fasta files.

Usage:
    run_msa.py ( --input_dir=PATH ) ( --output_dir=PATH ) [ --threads=INT ]

Options:
    --input_dir=PATH    The directory containing the fasta files.
    --output_dir=PATH   The directory to write the output to.
    --threads=INT       The number of threads to use. [default: 1]
"""

import os
import sys
import subprocess
import multiprocessing
import logging
from functools import partial
from pathlib import Path
logger = logging.getLogger(os.path.basename(__file__).replace('.py', ''))
from docopt import docopt

def run_mafft(fasta_file, output_dir, mafft_threads):
    """
    Runs mafft on a fasta file.
    """
    logger.info('Running mafft on {}'.format(fasta_file))
    output_file = os.path.join(output_dir, os.path.basename(fasta_file).replace('.fasta', '.aln.fasta'))
    cmd = [
        'mafft',
        '--thread', str(mafft_threads),
        '--quiet',
        '--auto',
        fasta_file,
        '>',
        output_file,
    ]

    cmd = ' '.join(cmd)

    logger.info('cmd: {}'.format(cmd))

    try:
        subprocess.run(cmd, check=True, shell=True)
        logger.info('Finished mafft on {}'.format(fasta_file))
    except subprocess.CalledProcessError as e:
        logger.error('Error running mafft on {}'.format(fasta_file))
        logger.error(e)
    except Exception as e:
        logger.error('Error running mafft on {}'.format(fasta_file))
        logger.error(e)

def main(*args, **kwargs):

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(name)s][%(asctime)s][%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
        ]
    )

    
    threads = int(kwargs['--threads'])
    mafft_threads = 4
    mafft_instances = max(1, int(threads/mafft_threads))

    logger.info('Starting run_msa.py with {} threads.'.format(mafft_instances))
    logger.info('Each mafft instance with {} threads.'.format(mafft_threads))

    with multiprocessing.Pool(processes=mafft_instances) as pool:
        pool.map(
            partial(run_mafft, output_dir=kwargs['--output_dir'], mafft_threads=mafft_threads),
            list(map(lambda x: str(x), Path(kwargs['--input_dir']).glob('*.fasta'))),
        )

    logger.info('ALL ALIGNMENTS FINISHED !!!')

if __name__ == '__main__':
    main(**docopt(__doc__))
