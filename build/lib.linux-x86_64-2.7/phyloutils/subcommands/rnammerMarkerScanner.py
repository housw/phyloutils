"""
scan rRNA markers using rnammer, for each fasta file in given folder.
"""

import os
import sys
import subprocess
import logging

_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """ add arguments
    """

    parser.add_argument("inputFastaFolder",
                            help = "input folder that contains all genome \
                                fasta files")
    parser.add_argument('-S', '--superkingdom',
                            help="superkingdom, choose arc/bac/euk, default='bac'",
                                choices=['arc', 'bac', 'euk'], default='bac')
    parser.add_argument('-m', '--molecule', help='type of rRNAs, can be any \
                                combination of the three, \
                                default is "tsu,lsu,ssu",', default='tsu,lsu,ssu',
                            choices=['<tsu,lsu,ssu>', '<tsu,lsu>', '<tsu,ssu>',
                                '<lsu,ssu>', '<tsu>', '<lsu>', '<ssu>'])
    parser.add_argument('-o', '--outdir', help='output directory, \
                                default=ScannedrRNAs',
                            default="ScannedrRNAs")
    parser.add_argument('-f', '--force',
                            help="force to override previous result?")
    parser.add_argument('-r', '--rnammer', help='path to rnammer excutable, \
                                default=~/Software/Module_Annotation/rnammer/rnammer-1.2/rnammer',
                            default='/home/hou/Software/Module_Annotation/rnammer/rnammer-1.2/rnammer')

    return parser


def command(args):
    """ a wrapper function to call rnammer for rRNA scanning, for each genome
        sequence in the input folder, it will scan the genome and found all rRNAs
    """
    rnammer = args.rnammer
    fasta_folder = args.inputFastaFolder
    kingdom = args.superkingdom
    type = args.molecule
    force = args.force
    output_folder = args.outdir
    prefix = args.prefix

    _logger.info("Creating output folder: " + output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])


    for fa in os.listdir(fasta_folder):
        if not (fa.endswith('.fa') or fa.endswith('.fasta') or fa.endswith('.fas')):
            _logger.info(" %s is not a fasta file, ignored "%fa)
            continue
        if prefix:
            filestem = prefix
        else:
            filestem = os.path.splitext(fa)[0]
        curr_outfolder = os.path.join(output_folder, filestem)
        subprocess.check_call(['mkdir', '-p', curr_outfolder])
        input_fa = os.path.join(fasta_folder, fa)
        out_rRNA = os.path.join(curr_outfolder, filestem+".fa")
        out_gff = os.path.join(curr_outfolder, filestem+".gff")
        out_hmm = os.path.join(curr_outfolder, filestem+".hmmsearch")

        # get orf
        _logger.info(" Getting rRNAs for %s ..."%fa)
        if os.path.exists(out_rRNA):
            if not force:
                _logger.warning(" rRNA file exists, nothing to do, \
                                           use --force to override")
                continue
            else:
                _logger.warning(" rRNA file for %s was override"%fa)
        subprocess.check_call([rnammer,
                                '-S', kingdom,
                                '-m', type,
                                '-f', out_rRNA,
                                '-gff', out_gff,
                                '-h', out_hmm,
                                input_fa
                              ])
