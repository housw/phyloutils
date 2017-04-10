#!/usr/bin/env python


from __future__ import division, print_function, absolute_import
import sys
import argparse
import logging
from phyloutils import subcommands
from phyloutils import __version__

__author__ = "housw"
__copyright__ = "housw"
__license__ = "GPLv3"

_logger = logging.getLogger(__name__)


class NotImplementedError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


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


def main(argv):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    # setup basic logging
    logging.basicConfig(level=logging.DEBUG, file=sys.stderr,
                    format="%(name)s:%(message)s")

    # parse subcomand and arguments
    command, arguments = parse_args(argv)

    # set up logger according to input arguments
    setup_logging(arguments.verbosity)

    # job dispatch
    _logger.debug("Starting job dispatch ...")
    return command(arguments)
    _logger.info("Job finished")


def parse_args(argv):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings, the first is subcommand

    Returns:
      :obj:`argparse.Namespace`: subcommand and arguments
    """

    # main parser
    parser = argparse.ArgumentParser(
        description="phyloutils: A set of subcommands for phylogenetic "\
        "computation.", prog="phyloutils")

    parser.add_argument('--version', action='version',
        version='phyloUtils {ver}'.format(ver=__version__))

    # parent parser, to specify shared arguments, inherited by subparsers
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-p", "--prefix", required=False,
        help="output prefix")
    parent_parser.add_argument('-v', '--verbosity', dest="verbosity",
        help="logging verbosity level, choose from ['debug', 'info', "\
            "'warning', 'error', 'critical'], default='info'",
        default='info')

    # subparsers
    subparsers = parser.add_subparsers(help='available subcommands', dest='subparser_name')

    # Add actions
    commands = {}
    for name, mod in subcommands.itermodules():
        subparser = subparsers.add_parser(name, help=mod.__doc__,
                description=mod.__doc__, parents=[parent_parser])
        subparser = mod.add_arguments(subparser)
        try:
            commands[name] = mod.command
        except Exception as err:
            _logger.error(" %s need to implement a command method!"%name)
            raise NotImplementedError(err)

    # display help
    display_help(sys.argv, parser)

    arguments = parser.parse_args(argv)
    arguments.argv = argv
    command = arguments.subparser_name

    return commands[command], arguments


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
