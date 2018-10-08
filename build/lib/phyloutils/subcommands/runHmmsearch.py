"""
wrapper for running hmmsearch
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError
from Bio import AlignIO

_logger = logging.getLogger(__name__)


def add_arguments(parser):

    parser.add_argument('input_hmm_file', help='input hmm file')
    parser.add_argument('input_file_or_folder', help='input file or folder')
    parser.add_argument('-o', '--outdir', help='output directory, default=\
                            ./hmmsearchOutput', default='hmmsearchOutput')
    parser.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                            with this suffix will be processed. default=.faa',
                        default='.faa')
    parser.add_argument('-e', '--evalue', type=float,
                        help='maximum sequence E-value threshold to \
                            report a hit in the output, default=1e-7',
                        default=1e-7)
    parser.add_argument('-m', '--dom_evalue', type=float,
                        help='maximum domain E-value threshold to \
                            report a hit in the output, default=1e-7',
                        default=1e-7)
    parser.add_argument('-t', '--not_write_tblout', action='store_true',
                        help='do not write table of per-sequence hits to file')
    parser.add_argument('-d', '--not_write_domtblout', action='store_true',
                        help='do not write table of per-domain hits to file')
    parser.add_argument('-f', '--not_write_pfamtblout', action='store_true',
                        help=' do not write table of hits and domains to file \
                           in Pfam format')
    parser.add_argument('--acc', action='store_true', help='prefer accessions\
                           over names in output')
    parser.add_argument('--noali', action='store_true', help='do not output \
                           alignments, so output is smaller')
    parser.add_argument('--cpu', type=int, help='number of CPU to use, \
                           default=10', default=10)
    parser.add_argument('-r', '--hmmsearch', help='path to hmmsearch excutable, \
                           default=~/Software/Module_Annotation/hmmer/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch',
                        default='/home/hou/Software/Module_Annotation/hmmer/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch')


def command(args):
    """wrapper function of hmmsearch
    """
    hmmsearch = args.hmmsearch
    input_file_or_folder = args.input_file_or_folder
    hmm_file = args.input_hmm_file
    outdir = args.outdir
    suffix = args.suffix
    evalue = args.evalue
    dom_evalue = args.dom_evalue
    write_tblout = not args.not_write_tblout
    write_domtblout = not args.not_write_domtblout
    write_pfamtblout = not args.not_write_pfamtblout
    acc = args.acc
    noali = args.noali
    cpu = args.cpu

    # create outdir
    subprocess.check_call(['mkdir', '-p', outdir])

    if os.path.isdir(input_file_or_folder):
        for f in os.listdir(input_file_or_folder):
            if f.endswith(suffix):
                if args.prefix:
                    prefix = os.path.join(args.outdir, args.prefix)
                else:
                    prefix = os.path.join(args.outdir, os.path.splitext(f)[0])

                # run hmmsearch
                input_fasta = os.path.join(input_file_or_folder, f)
                construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, prefix,
                                         evalue, dom_evalue, write_tblout,
                                         write_domtblout, write_pfamtblout,
                                         acc, noali, cpu)

    else:
        if args.prefix:
            prefix = os.path.join(args.outdir, args.prefix)
        else:
            basename = os.path.basename(input_file_or_folder)
            prefix = os.path.join(args.outdir, os.path.splitext(basename)[0])

        # run hmmsearch
        input_fasta = input_file_or_folder
        construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, prefix,
                                 evalue, dom_evalue, write_tblout,
                                 write_domtblout, write_pfamtblout,
                                 acc, noali, cpu)


def construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, output_prefix,
                            evalue=1e-7, dom_evalue=1e-7,
                            write_tblout=None, write_domtblout=None,
                            write_pfamtblout=None,
                            acc=None, noali=None, cpu=10, *args, **kwargs):
    """ construct hmmsearch commands
    """
    # hmmsearch -o /dev/null -Z 5000 -E 1e-7 --domE 1e-7 --domtblout output_dom_file \
    # --tblout output_seq_file --pfamtblout output_pfam_file hmm_file input_fasta
    template = '%s -o /dev/null -Z 5000 -E %s --domE %s --cpu %s'%(
        hmmsearch, str(evalue), str(dom_evalue), str(cpu))

    command = template.split(" ")
    if write_tblout:
        command.extend(['--tblout', output_prefix+".tblout"])
    if write_domtblout:
        command.extend(['--domtblout', output_prefix+".domtblout"])
    if write_pfamtblout:
        command.extend(['--pfamtblout', output_prefix+".pfamtblout"])
    if acc:
        command.extend(['--acc'])
    if noali:
        command.extend(['--noali'])

    command.extend([hmm_file, input_fasta])

    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "[runHmmsearch]: Error raised when running hmmsearch:%s"%e
        _logger.error(err)
        raise SubCommandError(err)
