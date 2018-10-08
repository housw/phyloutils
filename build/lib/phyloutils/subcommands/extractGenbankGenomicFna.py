"""
extract *_genomic.fna.gz file from each genome subfolder of downloaded GenbankGenomes into one folder.
"""

import logging
import os
import subprocess
import gzip


_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """Add arguments to the parser
    """
    parser.add_argument('input_GenbankGenomes', help="\
                                 input directory path to GenbankGenomes")
    parser.add_argument('-o', '--outdir',
                                   help="output directory, default='extractedGenbankGenomicFna'",
                                   default="extractedGenbankGenomicFna")
    return parser


def command(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        extract the *_genomic.fna records from each genome subfolder, copy them to the
        specified folder.
    """
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        _logger.info("Extracting genomic fna: "+genome_name)

        output_file = os.path.join(output_folder, genome_name+".fna")
        for f in os.listdir(os.path.join(input_folder, folder)):
            if f.endswith("_genomic.fna.gz") and ("_from_" not in f):
                gz_fna = os.path.join(input_folder, folder, f)
                #print(gz_fna)
                try:
                    with gzip.open(gz_fna, "r") as ih, open(output_file, "w") as oh:
                        for line in ih:
                            oh.write(line)
                except Exception as e:
                    err = "Error raised: %s"%e
                    _logger.error(err)
                break
        else:
            _logger.info("No protein file was found for %s"%genome_name)
            continue
