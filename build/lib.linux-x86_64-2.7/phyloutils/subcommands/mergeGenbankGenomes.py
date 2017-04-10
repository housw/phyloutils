"""
merge contigs/scaffolds/plasmids for each folder in GenbankGenomes into one fasta.
"""

import argparse
import logging
from Bio import SeqIO
import subprocess
import os
import sys
from . import merge_multiple_fasta_to_single_record


_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """ add arguments for parser
    """
    parser.add_argument('input_GenbankGenomes', help="\
                         input directory path to GenbankGenomes")
    parser.add_argument('-n', '--numberOfNs', type=int, default=100,
                         help="number of Ns used to seperate each fasta \
                             record, default=100")
    parser.add_argument('-o', '--outdir',
                        help="output directory, default='MergedGenbankGenomes'",
                        default="MergedGenbankGenomes")
    return parser


def command(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        merge fasta records for each assembly into one fasta record.
    """
    delim = args.numberOfNs*'N'
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    _logger.info("Creating output folder: "+output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        _logger.info("Merging genome sequence: "+genome_name)

        with open(os.path.join(output_folder, genome_name+".mfa"), "w") as oh:
            for f in os.listdir(os.path.join(input_folder, folder)):
                if f.endswith("_genomic.fna.gz") and \
                   ("_cds_from_genomic" not in f) and \
                   ("_rna_from_genomic" not in f):
                    gz_fna = os.path.join(input_folder, folder, f)
                    genome_seq = merge_multiple_fasta_to_single_record(gz_fna, delim)

                    # write out fasta record for current folder
                    if genome_seq:
                        oh.write(">"+genome_name+"\n")
                        oh.write(str(genome_seq)+"\n")
                    else:
                        _logger.info("No genome fna was found for %s"%genome_name)
