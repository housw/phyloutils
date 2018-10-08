"""
concatenate marker proteins found in each genome,
given an umbrella folder and proteins in genome subfolders."
"""

import argparse
import logging
import os
import subprocess

_logger = logging.getLogger(__name__)



def add_arguments(parser):

    parser.add_argument("inputUmbrellaFolder",
                            help="input folder that contains all genome subfolders,\
                               where extracted marker proteins were inside")
    parser.add_argument("inputMarkerList",
                            help='input marker list file that contains all marker names')
    parser.add_argument('-o', '--outdir',
                            help="output directory, default='MergedMarkerProteinSet'",
                            default="MergedMarkerProteinSet")


def command(args):
    """ given a list of marker names, concate all found marker proteins in each
        genome.
    """
    input_folder = args.inputUmbrellaFolder
    marker_list = args.inputMarkerList
    output_folder = args.outdir

    _logger.info("Creating output folder: " + output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])

    # read in marker names
    markers = []
    with open(marker_list, "r") as ih:
        for line in ih:
            if line.startswith('#'):
                continue
            line = line.strip()
            markers.append(line)

    # concate marker proteins
    for marker in markers:
        with open(os.path.join(output_folder, marker+".mfa"), "w") as oh:
            for genome in os.listdir(input_folder):
                for f in os.listdir(os.path.join(input_folder, genome)):
                    if marker in f and (f.endswith(".faa") or f.endswith(".fa") or f.endswith(".fasta")):
                        fasta = SeqIO.read(os.path.join(input_folder, genome, f), 'fasta')
                        oh.write(">"+genome+"\n")
                        oh.write(str(fasta.seq)+"\n")
