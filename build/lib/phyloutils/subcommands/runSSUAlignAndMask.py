"""
run ssu-align and mask for input rRNA sequences
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError

_logger = logging.getLogger(__name__)


def add_arguments(parser):

    parser.add_argument('inputFastaFile', help='input fasta\
                                file contains multiple rRNA sequences')
    parser.add_argument('-f', '--force', action='store_true',
                            help='override outputs')
    parser.add_argument('--no_align', action='store_true',
                            help='only search target \
                               sequence file for hits, skip alignment step')
    parser.add_argument('--no_search', action='store_true',
                            help='only align  target \
                                sequence file, skip initial search step')
    parser.add_argument('--output_RNA', action='store_true',
                            help='output alignment as RNA, default is DNA')
    parser.add_argument('--output_stk', action='store_true',
                            help='output stockholm alignments, \
                                default is fasta alignments')
    parser.add_argument('--kingdom', default='bacteria',
                            choices=['archaea', 'bacteria', 'eukarya'],
                            help='kingdom to search and align to, default is bacteria')
    parser.add_argument('--ssu_align', help='path to ssu-align excutable, \
                                default=~/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-align',
                            default='/home/hou/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-align')
    parser.add_argument('--ssu_mask', help='path to ssu-mask excutable, \
                                default=~/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-mask',
                            default='/home/hou/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-mask')
    return(parser)




def command(args):
    """ run ssu-align and ssu-mask for input sequences
    """
    ssu_align = args.ssu_align
    ssu_mask = args.ssu_mask
    input_fasta = args.inputFastaFile
    force = args.force
    prefix = args.prefix
    kingdom = args.kingdom
    no_align = args.no_align
    no_search = args.no_search
    output_RNA = args.output_RNA
    output_stk = args.output_stk

    if prefix:
        output_folder = prefix+"_ssu-align"
    else:
        basename = os.path.basename(input_fasta)
        filestem = os.path.splitext(basename)[0]
        output_folder = filestem + "_ssu-align"

    # do ssu alignment
    cmd = [ssu_align, '-n', kingdom]
    if force:
        cmd.append('-f')
    if not output_RNA:
        cmd.append('--dna')
    if no_align:
        cmd.append('--no-align')
    if no_search:
        cmd.append('--no-search')
    cmd.append(input_fasta)
    cmd.append(output_folder)
    try:
        p = subprocess.Popen(cmd)
        p.wait()
    except Exception as e:
        _logger.error("Error raised at running runSSUAlignAndMask:%s"%e)

    # do mask
    cmd=[ssu_mask]
    if not output_RNA:
        cmd.append('--dna')
    if not output_stk:
        cmd.append('--afa')
    cmd.append(output_folder)
    try:
        p = subprocess.Popen(cmd)
        p.wait()
    except Exception as e:
        err = "Error raised at running runSSUAlignAndMask:%s"%e
        _logger.error(err)
        raise SubCommandError(err)
