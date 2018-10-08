"""
merge contigs/scaffolds/plasmids for one mfa
"""

import argparse
import logging
from Bio import SeqIO
from . import merge_multiple_fasta_to_single_record

_logger = logging.getLogger(__name__)


def command(args):
    """ given a mfa file, call merge_multiple_fasta_to_single_record,
        the
    """
    input_mfa = args.input_mfa
    delim = args.numberOfNs*'N'
    header = args.header

    with open(header+".fa", "w") as oh:
        genome_seq = merge_multiple_fasta_to_single_record(input_mfa, delim)
        oh.write(">"+header+"\n")
        oh.write(str(genome_seq)+"\n")


def add_arguments(parser):
    """ add arguments for parser """
    parser.add_argument('input_mfa', help="\
                         input mfa file, contains multiple fasta records")
    parser.add_argument('header', help='header for \
                         merged contigs/scaffolds/plasmids')
    parser.add_argument('-n', '--numberOfNs', type=int, default=100,
                         help="number of Ns used to seperate each fasta \
                              record, default=100")
    return parser
