"""
change all the filenames and record names of fasta files,that inside of given folder.
"""


import argparse
import logging
import os
import subprocess
from Bio import SeqIO


_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """ add arguments
    """
    parser.add_argument("inputFastaFolder",
                        help="input folder that contains all fasta files")
    parser.add_argument("lookupTable",
                        help='tab-deliminated lookup table to convert \
                            old names to new names, in "old tab new" format,\
                            if more than two columns were detected in the \
                            lookup table, only the first two columns will be used.')
    parser.add_argument('-o', '--outdir',
                        help="output directory, default='renamedFasta'",
                        default="renamedFasta")
    parser.add_argument('-l', '--length', type=int, help='maximum \
                            length of allowed short names in the lookup table. \
                            Default=9.', default=9)
    return parser


def command(args):
    """given an input folder with fasta files inside, rename each fasta file and
       it's header name, accroding to the lookup table
    """
    input_folder = args.inputFastaFolder
    lookupTable = args.lookupTable
    output_folder = args.outdir
    length = args.length

    _logger.info("Creating output folder: "+ output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])

    # read in lookup table
    lookup_dict = {}
    newNames = set()
    with open(lookupTable, "r") as ih:
        for line in ih:
            if line.startswith('#'):
                continue
            line = line.strip().split("\t")
            k, v = line[0], line[1]
            try:
                assert v not in newNames, "repeat new name was not allowed: %s"%v
                assert len(v) <= length, "the length of %s was more than %d!"%(v, length)
            except Exception as e:
                _logger.error(e)
            newNames.add(v)
            lookup_dict[k] = v

    # rename fasta records
    for fa in os.listdir(input_folder):
        if fa in lookup_dict:
            newName = lookup_dict[fa]
            _logger.info("Renaming %s to %s"%(fa, newName))
            with open(os.path.join(output_folder, newName+".fa"), 'w') as oh:
                fa_rec = SeqIO.read(os.path.join(input_folder, fa), 'fasta')
                oh.write(">"+newName+"\n")
                oh.write(str(fa_rec.seq)+"\n")
        else:
            _logger.warning("No new name was found for %s, no change"%fa)
            subprocess.check_call(['mkdir', '-p', os.path.join(output_folder, "nochange")])
            subprocess.check_call(['cp', os.path.join(input_folder, fa),
                                   os.path.join(output_folder, "nochange", fa)
                                 ])
