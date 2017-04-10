"""
convert multiple sequence alignments in multiple fasta format to phylip-relaxed format.
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError
from Bio import AlignIO

_logger = logging.getLogger(__name__)


def convertFasta2Phylip(input_fasta, output_prefix, format='phylip-sequential'):
    """ convert input_fasta to phylip-relaxed format using Biopython
    """
    AlignIO.convert(input_fasta, 'fasta', output_prefix+".phy", format)


def add_arguments(parser):

    parser.add_argument('inputSuperAlignmentFolder',
                        help="input folder contains supermatrix alignments")
    parser.add_argument('-o', '--outdir',
                        help='output directory, default is input directory')
    parser.add_argument('--phylip_relaxed', action='store_true',
                        help='output phylip-relaxed\
                            format, default is phylip-sequential')
    parser.add_argument('-s', '--suffix', default='_supermatrix.msf',
                        help='suffix to identify input super alignment file, default=_supermatrix.msf')
    return parser


def command(args):
    """convert multiple sequence alignments in multiple fasta format to phylip-sequential format
    """
    input_folder = args.inputSuperAlignmentFolder
    suffix = args.suffix
    phylip_relaxed = args.phylip_relaxed

    if args.outdir:
        output_folder = args.outdir
    else:
        output_folder = input_folder

    outfmt = "phylip-sequential"
    if phylip_relaxed:
        outfmt = 'phylip-relaxed'

    # convert to phylip-relaxed format
    for f in os.listdir(input_folder):
        if f.endswith(suffix):
            input_fasta = os.path.join(input_folder, f)
            basename = os.path.splitext(f)[0]
            output_prefix = os.path.join(output_folder, basename)
            _logger.info(" Converting %s to %s format ..."%(f, outfmt))

            # it's compulsory to use phylip-sequential format for PartitionFinderProtein
            convertFasta2Phylip(input_fasta, output_prefix, format=outfmt)
