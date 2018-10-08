"""
generate tab-deliminated shortName table from NCBI Genome Assembly Report,
manually corrections might be needed."
"""

import argparse
import logging
import os

_logger = logging.getLogger(__name__)



def add_arguments(parser):
    """ add arguments to parser
    """
    parser.add_argument('input_file',
                           help='expored assembly report from NCBI Genome')
    parser.add_argument('-o', '--outdir',
                           help="output directory, default='./'",
                           default="./")
    parser.add_argument('-s', '--suffix',
                           help="suffix will be appended to the accession\
                           number, default='mfa'", default='mfa')


def command(args):
    """ generate a tab-deliminated lookUpTable of accession numbers and
        specie names
    """
    input_file = args.input_file
    outdir = args.outdir
    suffix = args.suffix
    prefix = args.prefix

    def _clean_name(sciName):
        """ replace characters to underscores or empty
        """
        sciName = sciName.replace(" ", "_").replace("'", "").replace("-", "_")
        sciName = sciName.replace("/", "_").replace("=", "_")
        sciName = sciName.replace("[", "").replace("]", "")
        sciName = sciName.replace("(", "_").replace(")", "_")
        sciName = sciName.replace("___", "_" ).replace("__", "_")
        return sciName

    # save all the names
    names = set()

    if prefix:
        output_file = prefix+"_shortNames.tab"
    else:
        basename = os.path.basename(input_file)
        filestem = os.path.splitext(basename)[0]
        output_file = filestem+"_shortNames.tab"

    # read in genome assembly report file
    with open(input_file, "r") as ih, open(output_file, "w") as oh:
        oh.write("#filename\tshortname\tgenus\tstrain\n")
        header = ih.readline().strip().lstrip("#").split("\t")
        try:
            genus_pos = header.index("Organism/Name")
            strain_pos = header.index("Strain")
            accession_pos = header.index("Assembly")
            Genbank_FTP_pos = header.index("GenBank FTP")
        except Exception as e:
            _logger.error(e)
        for line in ih:
            line = line.strip().split("\t")
            genus_full = line[genus_pos].strip().split(" ")[0].lstrip("[").rstrip("]")
            genus_short = genus_full[:5]
            accession = line[accession_pos].strip()
            assemblyID = line[Genbank_FTP_pos].strip().split("_")[-1] if "ASM" in line[-1].strip().split("_")[-1] else accession
            strain = line[strain_pos].strip().replace(" ", "")

            if strain:
                sciName = genus_short+"_"+strain+"_"+assemblyID
            else:
                sciName = genus_short+"_"+assemblyID
            sciName = _clean_name(sciName)

            accession += "."+suffix

            oh.write(accession+"\t"+sciName+"\t"+genus_full+"\t"+strain+"\n")
