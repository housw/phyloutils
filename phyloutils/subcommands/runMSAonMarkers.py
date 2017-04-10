"""
given a folder with marker proteins inside,
run multiple resquence alignment on each marker.
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError

_logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument('inputMarkerProteinsFolder',
                                   help="input folder contains identified marker proteins")
    parser.add_argument('-o', '--outdir',
                                   help='output directory, default is input directory')
    parser.add_argument('-s', '--suffix', default='.mfa',
                                   help='suffix to identify input protein fasta file, default=.mfa')


def command(args):
    """ given an input folder contains marker protein sequences, run multiple sequence alignment
        on each marker
    """
    input_folder = args.inputMarkerProteinsFolder
    suffix = args.suffix
    if args.outdir:
        output_folder = args.outdir
    else:
        output_folder = input_folder

    for f in os.listdir(input_folder):
        if f.endswith(suffix):
            filestem = os.path.splitext(f)[0]
            # run msa
            try:
                align = subprocess.Popen([ 'muscle',
                                        '-in', os.path.join(input_folder, f),
                                        '-out', os.path.join(output_folder, filestem+".msf")
                                        ])
                align.wait()

                trim = subprocess.Popen(['trimal',
                                         '-in', os.path.join(output_folder, filestem+".msf"),
                                         '-out', os.path.join(output_folder, filestem+"_trimAl.msf"),
                                         '-gt', '0.9',
                                         '-cons', '60',
                                         '-w', '3'
                                        ])
                trim.wait()
            except Exception as e:
                err = " Error raised at runMSAonMarkers: %s"%e
                _logger.error(err)
                raise SubCommandError(err)
