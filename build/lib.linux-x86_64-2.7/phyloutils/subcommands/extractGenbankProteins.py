"""
extract *_protein.faa.gz file from each genome subfolder of downloaded GenbankGenomes into one folder.
"""

import argparse
import logging
import os
import subprocess
import gzip

_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """ add arguments
    """
    parser.add_argument('input_GenbankGenomes', help="\
                        input directory path to GenbankGenomes")
    parser.add_argument('-o', '--outdir',
                        help="output directory, default='extractedGenbankProteins'",
                        default="extractedGenbankProteins")
    return parser


def command(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        extract the protein records from each genome subfolder, copy them to the
        specified folder.
    """
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    _logger.info("Creating output folder: "+ output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        _logger.info("Extracting genome proteins: "+genome_name)

        output_file = os.path.join(output_folder, genome_name+".faa")
        for f in os.listdir(os.path.join(input_folder, folder)):
            if f.endswith("_protein.faa.gz"):
                gz_faa = os.path.join(input_folder, folder, f)
                try:
                    with gzip.open(gz_faa, "r") as ih, open(output_file, "w") as oh:
                        for line in ih:
                            oh.write(line)
                except Exception as e:
                    err = "Error raised: %s"%e
                    _logger.error(err)
            break
        else:
            _logger.warning("No protein file was found for %s"%genome_name)
