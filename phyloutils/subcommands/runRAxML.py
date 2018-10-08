"""
wrapper for running RAxML, with or without partitions
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError

_logger = logging.getLogger(__name__)


def add_arguments(parser):

    parser.add_argument('input_alignment', help='input alignment in phylip format')
    parser.add_argument('-r', '--raxml', help='path to raxml excutable, \
                                   default=~/Software/Module_Phylogeny/RAxML/standard-RAxML/raxmlHPC',
                               default='/home/hou/Software/Module_Phylogeny/RAxML/standard-RAxML/raxmlHPC')
    parser.add_argument("--isDNA", action="store_true",
                               help="are DNA sequence alignment? default is Protein")
    parser.add_argument('-m', '--model', help='model will be directly passed to RAxML, \
                                   default is GTRGAMMA for DNA, PROTGAMMAWAG for Protein.')
    parser.add_argument('-#', '--bootstraps', type=int, help='number of bootstraps, default=100',
                               default=100)
    parser.add_argument('-o', '--outgroup', type=str,
                               help='sequence names of outgroups, comma-delimited')
    parser.add_argument('--slow', action='store_true',
                               help = 'do slow bootstrapping instead? default is rapid bootstrapping')
    parser.add_argument('-S', '--Structure', help='path to secondary structure file')
    parser.add_argument('-q', '--Partition', help='path to partition file')
    parser.add_argument('--suffix', help='suffix string')
    parser.add_argument('--threads', type=int, default=30,
                               help='number of threads, default=30')
    return parser



def runRAxML_rapid_bootstrap(raxml, alignment, model, bootstraps, outgroup,
                             suffix, structure, partition, threads):
    command = [raxml, '-s', alignment, '-m', model, '-f', 'a', '-p', '12345',
               '-x', '12345', '-#', str(bootstraps), '-T', str(threads), '-n', suffix]
    if outgroup:
        command.extend(['-o', outgroup])
    if structure:
        command.extend(['-S', structure])
    if partition:
        command.extend(['-q', partition])
    try:
        _logger.info(" ".join(command))
        p = subprocess.Popen(command)
        p.wait()
    except Exception as e:
        err = "[runRAxML_rapid_bootstrap]: Error raised at runRAxML_rapid_bootstrap:%s"%e
        _logger.error(err)
        raise SubCommandError(err)

def runRAxML_slow_bootstrap(raxml, alignment, model, bootstraps, outgroup,
                            suffix, structure, partition, threads):
    template = '%s -s %s -p 12345 -x 12345 -m %s -n %s -T %s'%(
    raxml, alignment, model, suffix, str(threads))

    # step1: 20 ML search(-#), from a parsimony random seed(-p)
    command = template.split(" ")
    command.extend(['-#', '20'])
    if outgroup:
        command.extend(['-o', outgroup])
    if structure:
        command.extend(['-S', structure])
    if partition:
        command.extend(['-q', partition])
    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "Error raised at step1 of runRAxML_slow_bootstrap:%s"%e
        _logger.error(err)
        raise SubCommandError(err)

    # step2: 100 bootstrap search
    #/home/hou/Software/RAxML/standard-RAxML/raxmlHPC-PTHREADS -m PROTGAMMAILGF -p 12345 -b 12345 -# 100 -s tcoffee_mafft_fileset_trimAl.phy -T 21 -n PROTGAMMAILGF_100BS

    # step3: projection, draw bipartitions on the best ML tree using bootstrap replicate trees
    # for unrooted trees the correct representation is actually the one with support values assigned to branches and not nodes
    #/home/hou/Software/RAxML/standard-RAxML/raxmlHPC-PTHREADS -m PROTGAMMAILGF -p 12345 -f b -t RAxML_bestTree.PROTGAMMAILGF_20ML -z RAxML_bootstrap.PROTGAMMAILGF_100BS -n PROTGAMMAILGF_projection


def command(args):
    alignment = args.input_alignment
    raxml = args.raxml
    slow = args.slow
    bootstraps = args.bootstraps
    isDNA = args.isDNA
    if args.model:
        model = args.model
    else:
        if isDNA:
            model = "GTRGAMMA"
        else:
            model = "PROTGAMMAWAG"
    outgroup = args.outgroup
    structure = args.Structure
    partition = args.Partition
    suffix = args.suffix if args.suffix else str(int(time.time()))
    threads = args.threads

    if not slow:
        runRAxML_rapid_bootstrap(raxml, alignment, model, bootstraps, outgroup,
                                 suffix, structure, partition, threads)
    else:
        err = "[runRAxML]: Not implemented yet!!!"
        _logger.error(err)
        raise SubCommandError(err)
        runRAxML_slow_bootstrap(raxml, alignment, model, bootstraps, outgroup,
                                 suffix, structure, partition, threads)
