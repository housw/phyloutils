"""
replace shortname back to scientific names in input tree file,
it will only replace the shortname strings found the tree file, should not change
anything else like topology/branchlength/bootstraps of the tree.
"""

import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError

_logger = logging.getLogger(__name__)


def add_arguments(parser):

    parser.add_argument('inputTreeFile',
                            help="input tree file generated using RAxML, in newick/nexus format")
    parser.add_argument('-f', "--format", choices=['newick', 'nexus'],
                            help="input tree format, default=newick", default='newick')
    parser.add_argument('lookUpTable',
                            help='tab-deliminated lookup table to convert \
                               short names to scientific names, this is the same file \
                               used in renameFastaFile, the second column should \
                               be shortname and third column should be scientificname, \
                               no comma/space/(/)/colon were allowed in scientificname.')
    return parser


def command(args):
    """ replace the shortname back to scientificname in the tree file, to make
        it easy for publication
    """
    inputTreeFile = args.inputTreeFile
    lookUpTable = args.lookUpTable
    format = args.format
    prefix = args.prefix

    if prefix:
        output_file = prefix+"_sci.tre"
    else:
        basename = os.path.basename(inputTreeFile)
        filestem = os.path.splitext(basename)[0]
        output_file = filestem+"_sci.tre"

    # read in lookUpTable
    short2sci = {} # {shortname: scientificname}
    with open(lookUpTable, "r") as ih:
        for line in ih:
            if line.startswith('#'):
                continue
            line = line.strip().split("\t")
            if len(line) <=2:
                _logger.warning(" No scientificname were found at: %s"%("\t".join(line)))
                continue
            shortname = line[1].strip()
            scientificname = line[2].strip()
            short2sci[shortname] = scientificname

    used = {} # {shortname: True} record this has been replaced or not
    with open(inputTreeFile, "r") as ih, open(output_file, "w") as oh:
        if format == 'newick':
            for line in ih:
                for short, sci in short2sci.iteritems():
                    # a taxon name should followed by other attributes like branchlength, seperate by colon
                    if short+":" in line:
                        if short not in used:
                            used[short] = True
                            #print line
                            line = line.replace(short+":", sci+":")
                    else:
                        _logger.warning(" %s has been used, and should be unique in tree!"%short)
                        continue
                oh.write(line)
        else:
            for line in ih:
                for short, sci in short2sci.iteritems():
                    if short in line:
                        line = line.replace(short, sci)
                oh.write(line)

        _logger.info(" Please check the names and manually correct wrong names!!! ")
