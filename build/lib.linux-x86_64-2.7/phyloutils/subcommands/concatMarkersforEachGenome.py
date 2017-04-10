"""
given a folder with trimmed marker protein alignment inside,
concatenate markers for each genome, gaps will be placed if no marker found for one genome.
"""


import argparse
import logging
import os
import subprocess
import sys
from . import SubCommandError


_logger = logging.getLogger(__name__)


def add_arguments(parser):
    parser.add_argument('inputMarkerAlignmentFolder',
                            help="input folder contains marker protein alignments")
    parser.add_argument('-o', '--outdir',
                            help='output directory, default is input directory')
    parser.add_argument('-s', '--suffix', default='_trimAl.msf',
                            help='suffix to identify input trimmed protein alignment file, default=_trimAl.msf')


def command(args):
    """given a folder with marker protein alignments inside, concatenate markers
       for each genome, gaps will be placed if no marker found for one genome.
    """
    input_folder = args.inputMarkerAlignmentFolder
    suffix = args.suffix
    if not args.outdir:
        output_folder = input_folder
    else:
        output_folder = args.outdir

    # summarize all markers, marker length and all genome names
    genome2markers = {} # {genome:{marker1:seq1, marker2:seq2, ...}}
    marker2length = {} # {marker1:length1, marker2:length2}
    for f in os.listdir(input_folder):
        if f.endswith(suffix):
            marker = f.split("_")[0]
            aln_records = SeqIO.parse(os.path.join(input_folder, f), 'fasta')
            for aln in aln_records:
                genome = aln.name
                seq = str(aln.seq).strip()
                marker_len = len(seq)
                if marker not in marker2length:
                    marker2length[marker] = marker_len
                else:
                    assert marker_len == marker2length[marker], "Error: Are you\
                    sure the input file is an alignment? The length is different!"
                if genome not in genome2markers:
                    genome2markers[genome] = {marker:seq}
                else:
                    genome2markers[genome].update({marker:seq})

    # write out concated markers for each genome
    all_genomes = sorted(genome2markers.keys())
    all_markers = sorted(marker2length.keys())
    with open(os.path.join(output_folder, "concatMarkersforEachGenome_supermatrix.msf"), 'w') as oh:
        for genome in all_genomes:
            header = ">"+genome+"\n"
            seq = ""
            for marker in all_markers:
                marker_seq = genome2markers[genome].get(marker, None)
                if not marker_seq:
                    marker_seq = "-"*marker2length[marker]
                seq += marker_seq.strip()
            oh.write(header)
            oh.write(seq+"\n")

    # write out marker partitions
    with open(os.path.join(output_folder, "concatMarkersforEachGenome_supermatrix.coord"), 'w') as coord:
        accum_start = 0
        for marker in all_markers:
            start = accum_start + 1
            stop = start + marker2length[marker]-1
            coord.write("%s = %d-%d\n"%(marker, start, stop))
            accum_start = stop
