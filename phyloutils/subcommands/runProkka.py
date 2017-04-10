"""
wrapper for running prokka
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
                            ./prokkaOutput', default='prokkaOutput')
    parser.add_argument('-a', '--parameter_file', help="tab-deliminated "\
                            "parameter file, contains fileanme, genus and strain")
    parser.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                            with this suffix will be processed. default=.fa',
                        default='.fa')
    parser.add_argument('-c', '--cpu', type=int,
                        help='Number of CPUs to use, [0=all] default=8.',
                        default=8)
    parser.add_argument('--nogenes', action='store_true',
                        help="Do not add 'gene' features for each 'CDS' feature")
    parser.add_argument('--donotforce', action='store_true',
                        help="Do not force remove existing folders")
    parser.add_argument('--addmrna', action='store_true',
                        help="Add 'mRNA' features for each 'CDS' feature (default OFF)")
    parser.add_argument('--locustag',
                        help="Locus tag prefix (default 'PROKKA' or filestem of \
                            input files when input is a folder).",
                        default='PROKKA')
    parser.add_argument('--increment', type = int,
                        help="Locus tag counter increment (default '1').",
                        default=1)
    parser.add_argument('--kingdom',
                        help="Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')",
                        default='Bacteria')
    parser.add_argument('--ignoreError', action="store_true",
                        help="ignore Errors raised when annotating multiple genomes?")
    parser.add_argument('-r', '--prokka', help='path to prokka excutable, \
                           default=~/Software/Module_Annotation/prokka/prokka-1.12/bin/prokka',
                        default='/home/hou/Software/Module_Annotation/prokka/prokka-1.12/bin/prokka')

    return parser



def command(args):
    """ wrapper function of prokka
    """
    prokka = args.prokka
    input_file_or_folder = args.input_file_or_folder
    outdir = args.outdir
    parameter_file = args.parameter_file
    suffix = args.suffix
    cpu = args.cpu
    addgenes = not args.nogenes
    addmrna = args.addmrna
    locustag = args.locustag
    increment = args.increment
    kingdom = args.kingdom
    force = not args.donotforce
    ignoreError = args.ignoreError

    # read in parameter_file, if exists
    file2parameters = {}
    if parameter_file:
        with open(parameter_file, "r") as ih:
            header = ih.readline().lstrip("#").strip().split("\t")
            filename_pos = header.index("filename")
            genus_pos = header.index("genus")
            strain_pos = header.index("strain")
            for line in ih:
                line = line.strip().split("\t")
                filename = line[filename_pos]
                genus = line[genus_pos]
                strain = line[strain_pos]
                if filename not in file2parameters:
                    file2parameters[filename] = {"filename":filename,
                                                 "genus":genus,
                                                 "strain":strain}
                else:
                    _logger.warning(" found sequence %s more than once!"%filename)
                    continue

    if os.path.isdir(input_file_or_folder):
        for f in os.listdir(input_file_or_folder):
            if f.endswith(suffix):
                if not args.prefix:
                    prefix = os.path.splitext(f)[0]
                    if f in file2parameters:
                        genus = file2parameters[f].get('genus', 'genus')
                        species = file2parameters[f].get('species', 'species')
                        strain = file2parameters[f].get('strain', 'strain')
                        locustag = strain
                    else:
                        locustag = prefix
                        genus = "genus"
                        species = "species"
                        strain = "strain"

                # run prokka
                input_fasta = os.path.join(input_file_or_folder, f)
                construct_prokka_cmds(prokka, input_fasta, outdir, prefix,
                                      locustag, genus, species, strain, increment,
                                      kingdom, cpu, addgenes, addmrna, force, ignoreError)
    else:
        if not args.prefix:
            basename = os.path.basename(input_file_or_folder)
            prefix = os.path.splitext(basename)[0]

            if f in file2parameters:
                genus = file2parameters[f].get('genus', 'genus')
                species = file2parameters[f].get('species', 'species')
                strain = file2parameters[f].get('strain', 'strain')
                locustag = strain
            else:
                locustag = prefix
                genus = "genus"
                species = "species"
                strain = "strain"

        # run prokka
        input_fasta = input_file_or_folder
        construct_prokka_cmds(prokka, input_fasta, outdir, prefix,
                              locustag, genus, species, strain, increment,
                              kingdom, cpu, addgenes, addmrna, force, ignoreError)


def construct_prokka_cmds(prokka, input_fasta, outdir, output_prefix,
                            locustag="PROKKA", genus="genus", species="species",
                            strain="strain", increment=1, kingdom="Bacteria", cpu=8,
                            addgenes=True, addmrna=False, force=True, ignoreError=False,
                             *args, **kwargs):
    """ construct prokka commands
    """
    template = "%s %s --outdir %s --prefix %s \
    --cpus %d --locustag %s --kingdom %s --increment %d"%(prokka, input_fasta, outdir, output_prefix,
    cpu, locustag, kingdom, increment)

    command = template.split(" ")
    if genus != "genus":
        command.extend(['--genus', genus, '--usegenus'])
    if species != "species":
        command.extend(['--species', species])
    if strain != "strain":
        command.extend(['--strain', strain])
    if addgenes:
        command.extend(['--addgenes'])
    if addmrna:
        command.extend(['--addmrna'])
    if force:
        command.extend(['--force'])
    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "[runProkka]: Error raised when running prokka:%s"%e
        _logger.error(err)
        if not ignoreError:
            raise SubCommandError(err)
