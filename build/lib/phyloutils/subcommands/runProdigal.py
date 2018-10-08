"""
wrapper for running prodigal
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

    parser.add_argument('input_file_or_folder', help='input file or folder')
    parser.add_argument('-o', '--outdir', help='output directory, default=\
                                    ./prodigalOutput', default='prodigalOutput')
    parser.add_argument('-a', '--not_write_proteins', action='store_true',
                                    help='do not write predicted proteins?')
    parser.add_argument('-d', '--not_write_DNAs', action='store_true',
                                    help='do not write nucleotide sequences of predicted proteins?')
    parser.add_argument('-f', '--format', choices=['gbk', 'gff'], help='output format, \
                                   for gene annotations, can be gff or gbk, default is gff',
                                   default='gff')
    parser.add_argument('-g', type=int, help='specify a translation table to use, default=11',
                                    default=11)
    parser.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                                    with this suffix will be processed. default=.fa',
                                    default='.fa')
    parser.add_argument('-m', action='store_true', help='Treat runs of N as masked sequence; \
                                    do not build genes across them.')
    parser.add_argument('-c', action='store_true', help='Closed ends. \
                                    Do not allow genes to run off edges.')
    parser.add_argument('-n', action='store_true', help='Bypass Shine-Dalgarno\
                                    trainer and force a full motif scan.')
    parser.add_argument('-r', '--prodigal', help='path to prodigal excutable, \
                                   default=~/Software/Module_Annotation/prodigal/prodigal',
                                   default='/home/hou/Software/Module_Annotation/prodigal/prodigal')

    return parser


def command(args):
    prodigal = args.prodigal
    input_file_or_folder = args.input_file_or_folder
    suffix = args.suffix

    output_fmt = args.format
    m = args.m
    c = args.c
    n = args.n

    write_proteins = not args.not_write_proteins
    write_DNAs = not args.not_write_DNAs

    # create outdir
    subprocess.check_call(['mkdir', '-p', args.outdir])

    if os.path.isdir(input_file_or_folder):
        for f in os.listdir(input_file_or_folder):
            if f.endswith(suffix):
                if args.prefix:
                    prefix = os.path.join(args.outdir, args.prefix)
                else:
                    prefix = os.path.join(args.outdir, os.path.splitext(f)[0])

                # run prodigal
                input_fasta = os.path.join(input_file_or_folder, f)
                construct_prodigal_cmds(prodigal, input_fasta, prefix, output_fmt,
                                        write_proteins, write_DNAs, m, c, n)

    else:
        if args.prefix:
            prefix = os.path.join(args.outdir, args.prefix)
        else:
            basename = os.path.basename(input_file_or_folder)
            prefix = os.path.join(args.outdir, os.path.splitext(basename)[0])

        # run prodigal
        input_fasta = input_file_or_folder
        construct_prodigal_cmds(prodigal, input_fasta, prefix, output_fmt,
                                write_proteins, write_DNAs, m, c, n)


def construct_prodigal_cmds(prodigal, input_fasta, output_prefix,
                            output_fmt='gff', write_proteins=None, write_DNAs=None,
                            m=None, c=None, n=None, *args, **kwargs):
    """ construct prodigal commands
    """
    template = '%s -i %s -o %s'%(prodigal, input_fasta, output_prefix+"."+output_fmt)

    command = template.split(" ")
    if write_proteins:
        command.extend(['-a', output_prefix+".faa"])
    if write_DNAs:
        command.extend(['-d', output_prefix+".fsa"])
    if m:
        command.extend(['-m'])
    if c:
        command.extend(['-c'])
    if n:
        command.extend(['-c'])

    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "[runProdigal]: Error raised when running prodigal:%s"%e
        raise SubCommandError(err)
