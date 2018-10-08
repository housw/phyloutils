"""
merge marker rRNAs found in each genome into one file,
given an umbrella folder and rRNA sequences in genome subfolders.
"""

import argparse
import logging
import os
import subprocess


_logger = logging.getLogger(__name__)


def add_arguments(parser):

    parser.add_argument("inputUmbrellaFolder",
                            help="input folder that contains all genome subfolders,\
                               where extracted marker rRNAs were inside")
    parser.add_argument("-m", '--molecule', default='16s',
                            choices=['5s', '16s', '23s', '8s', '18s', '28s'],
                            help='molecule type you want to merge, default=16s')
    parser.add_argument('-s', '--strategy', default='all',
                            choices = ['all', 'best'], help = 'strategy to select\
                               marker rRNA when multiple markers were found, default=all')
    parser.add_argument('-l', '--length_cutoff', default=1200, type=int,
                            help='length cutoff, to exclude short sequences, default=1200')
    parser.add_argument('-n', '--numberOfNs', default=50, type=int,
                            help='maximum Ns allowed in the sequences, default=50')
    parser.add_argument('-o', '--outdir',
                            help="output directory, default='MergedrRNASet'",
                            default="MergedrRNASet")
    return parser


def command(args):
    """ given the name of marker rRNA, concate this marker rRNA for all genomes.
    """
    input_folder = args.inputUmbrellaFolder
    marker = args.molecule+"_rRNA"
    strategy = args.strategy       # all or best
    output_folder = args.outdir
    length_cutoff = args.length_cutoff
    numberOfNs = args.numberOfNs

    _logger.info("Creating output folder: " + output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])

    # merge marker rRNAs
    with open(os.path.join(output_folder, marker+"_"+strategy+".fa"), "w") as oh:
        for genome_folder in os.listdir(input_folder):
            if not os.path.isdir(os.path.join(input_folder, genome_folder)):
                continue
            for f in os.listdir(os.path.join(input_folder, genome_folder)):
                if f.endswith(".fa") or f.endswith(".fasta") or f.endswith(".fas"):
                    genome_name = genome_folder
                    marker_seq = os.path.join(input_folder, genome_folder, f)
                    if strategy == 'best':
                        best = None
                        best_score = 0
                        for rec in SeqIO.parse(marker_seq, "fasta"):
                            if marker in rec.description:
                                if len(rec.seq) < length_cutoff:
                                    continue
                                if rec.seq.count('N') + rec.seq.count('n') >= numberOfNs:
                                    continue
                                if not best:
                                    best = rec
                                    best_score = float(rec.description.split("=")[-1])
                                else:
                                    current = rec
                                    current_score = float(rec.description.split("=")[-1])
                                    if current_score > best_score:
                                        best_score = current_score
                                        best = current
                        if best:
                            if best.seq.count('N') >=50:
                                continue
                            oh.write(">"+genome_name+"\n")
                            oh.write(str(best.seq)+"\n")
                    else:
                        count = 0
                        for rec in SeqIO.parse(marker_seq, 'fasta'):
                            seq_set = set()
                            if marker in rec.description:
                                if rec.seq.count('N') + rec.seq.count('n') >= numberOfNs:
                                    continue
                                if len(rec.seq) < length_cutoff:
                                    continue
                                if str(rec.seq) in seq_set:
                                    continue
                                count += 1
                                seq_set.add(str(rec.seq))
                                oh.write(">"+genome_name+"_"+str(count)+"\n")
                                oh.write(str(rec.seq)+"\n")
