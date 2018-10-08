#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


import sys
import argparse
import logging


_logger = logging.getLogger(__name__)


class NotImplementedError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class IndexNotFoundError(Exception):
    pass


def setup_logging(loglevel):
    """Setup basic loggings

    Args:
      loglevel (str): minimum loglevel for emitting messages
    """

    loglevel = {
        'critical': logging.CRITICAL,
        'error'   : logging.ERROR,
        'warning' : logging.WARNING,
        'info'    : logging.INFO,
        'debug'   : logging.DEBUG,
    }.get(loglevel, logging.DEBUG)

    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"

    logging.basicConfig(level=loglevel, file=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main_usage(parser):
    """ display usage for main parser
    """
    _logger.error(parser.format_help())


def subparser_usage(argv, parser):
    """ display usage for subparser
    """
    cmd = argv[1]
    found = 0
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for choice, subparser in action.choices.items():
                if cmd == choice:
                    _logger.error(subparser.format_help())
                    found = 1
    if not found:
        if cmd in ("--version", "-h", "--help"):
            args = parser.parse_args()
        else:
            _logger.error("\n\nERROR:%s is not a valid command!!!\n\n" % cmd)
            main_usage(parser)


def display_help(argv, parser):
    """ display help information
    """
    if len(argv) == 1:
        main_usage(parser)
        sys.exit(1)
    elif len(argv) == 2:
        subparser_usage(argv, parser)
        sys.exit(1)
    else:
        pass