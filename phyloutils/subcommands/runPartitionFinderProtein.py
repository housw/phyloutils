"""
estimate models for partitions using PartitionFinderProtein
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError

_logger = logging.getLogger(__name__)


def add_arguments(parser):

    PartitionFinderProtein = '/home/hou/Software/Module_Phylogeny/PartitionFinder/PartitionFinderV1.1.1/PartitionFinderProtein.py'

    parser.add_argument('inputPhylipSuperAlign',
                            help="input supermatrix alignment in phylip-sequential format")
    parser.add_argument('inputCoordFile',
                            help='input coordinate file contains the positions of proteins')
    parser.add_argument('newWorkingFolder',
                            help='new working directory, to put alignment and config files')
    parser.add_argument('--PartitionFinderProtein', help="path to\
                                PartitionFinderProtein excutable, default=%s"%PartitionFinderProtein,
                            default=PartitionFinderProtein)
    return parser


def command(args):
    """ run PartitionFinderProtein on concatenated super alignment, find the best
        protein models and partitions.
    """
    alignment = args.inputPhylipSuperAlign
    coordinate = args.inputCoordFile
    workingdir = args.newWorkingFolder
    program = args.PartitionFinderProtein

    partition_finder_template = """## ALIGNMENT FILE ##

alignment = %s;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | <list> ##
##              for PartitionFinderProtein: all_protein | <list> ##
models = all_protein;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]
%s

## SCHEMES, search: all | user | greedy ##
[schemes]
search = greedy;

#user schemes go here if search=user. See manual for how to define.#
    """

    # copy alignment to workingdir
    subprocess.check_call(['mkdir', '-p', workingdir])
    subprocess.check_call(['cp', alignment, workingdir])
    alignment_basename = os.path.basename(alignment)

    # fill up partition_finder_template, copy to workingdir
    coord_str = ""
    with open(coordinate, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            coord_str += line.strip() + ";\n"
    config_file = os.path.join(workingdir, "partition_finder.cfg")
    with open(config_file, "w") as oh:
        partition_finder_template = partition_finder_template%(alignment_basename, coord_str)
        oh.write(partition_finder_template)

    # run PartitionFinderProtein
    #python ~/Software/Module_Phylogeny/PartitionFinder/PartitionFinderV1.1.1/PartitionFinderProtein.py ./aminoacid/ --raxml --force-restart -p 30
    try:
        p = subprocess.Popen([ 'python', program,
                            workingdir,
                            '--raxml',
                            '--force-restart',
                            '-p', '30'
                         ])
        p.wait()
    except Exception as e:
        err = "[runPartitionFinderProtein]:Error raised at running PartitionFinderProtein:%s"%e
        _logger.error(err)
        raise SubCommandError(err)
