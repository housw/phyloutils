"""
scan markers using hmmsearch, for each fasta file in given folder, a python implementation of AMPHORE2.
"""

import argparse
import logging
import os
import subprocess

_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """ add arguments
    """
    parser.add_argument("inputFastaFolder",
                            help="input folder that contains all genome \
                                fasta files")
    parser.add_argument("inputMarkerFile",
                            help='input marker file that contains hmm models')
    parser.add_argument('-e', '--evalue', default=1e-7,
                            help="evalue cutoff, default=1e-7")
    parser.add_argument('-o', '--outdir',
                            help="output directory, default='ScannedMarkers'",
                            default="ScannedMarkers")
    parser.add_argument('-f', '--force',
                            help="force to override previous result?")
    return parser


def command(args):
    """ a simplified python implementation of MarkerScanner.pl of AMPHORE2, given
        a folder with fasta sequences inside, a hmm file with hmm models inside,
        search markers within input fasta sequences
    """
    fasta_folder = args.inputFastaFolder
    input_marker = args.inputMarkerFile
    force = args.force
    evalue = args.evalue
    output_folder = args.outdir

    _logger.info("Creating output folder: " + output_folder)
    subprocess.check_call(['mkdir', '-p', output_folder])


    for fa in os.listdir(fasta_folder):
        if not (fa.endswith('.fa') or fa.endswith('.fasta') or fa.endswith('.fas')):
            _logger.info(" %s is not a fasta file, ignored"%fa)
            continue
        filestem = os.path.splitext(fa)[0]
        curr_outfolder = os.path.join(output_folder, filestem)
        subprocess.check_call(['mkdir', '-p', curr_outfolder])
        input_fa = os.path.join(fasta_folder, fa)
        out_orf = os.path.join(curr_outfolder, filestem+".orf")
        out_hmm = os.path.join(curr_outfolder, filestem+".hmmsearch")

        # get orf
        _logger.info(" Getting orfs of %s ..."%fa)
        if os.path.exists(out_orf):
            if not force:
                _logger.warning(" ORF exists, nothing to do, \
                                           use --force to override")
                continue
        subprocess.check_call(['getorf',
                                '-sequence', input_fa,
                                '-outseq', out_orf,
                                '-table', '1',
                                '-minsize', '100'
                              ])
        # run hmmsearch
        _logger.info(" hmmsearch of markers for %s ..."%fa)
        subprocess.check_call(['hmmsearch',
                                '-Z', '5000',
                                '--cpu', '40',
                                '-E', str(evalue),
                                '--domE', str(evalue),
                                '--domtblout', out_hmm,
                                '-o', '/dev/null',
                                input_marker,
                                out_orf
                              ])

        # get protein fasta belonging to each marker
        _logger.info(" write marker proteins for %s ..."%fa)
        get_fasta_from_hmmsearch_result(out_hmm, out_orf, curr_outfolder)


def get_fasta_from_hmmsearch_result(hmmsearch_domtblout, query_fasta, output_folder):
    """given a query fasta, and hmmsearch_domtblout, find the best hit for each
       marker, and write it out to output_folder with marker name as prefix
    """

    marker2queries = {} # {marker:[query1, query2]}
    query_marker_evalue = {} # {query: {marker: evalue}}
    query_match = {} # {marker:{query1:query_match_len}}
    marker_match = {} # {marker:{query1:marker_match_len}}
    len_dict = {} # {query: query_length, marker:marker_length}

    with open(hmmsearch_domtblout, "r") as ih:
        for line in ih:
            #print line
            if line.startswith("#"):
                continue
            line = line.strip().split()
            query = line[0]
            query_len = int(line[2])
            marker = line[3]
            marker_len = int(line[5])
            full_evalue = float(line[6])
            #dom_evalue = float(line[11])
            marker_match_len = int(line[16]) - int(line[15])
            query_match_len = int(line[18]) - int(line[17])

            len_dict[query] = query_len
            len_dict[marker] = marker_len
            if marker not in query_match:
                query_match[marker] = {}
                query_match[marker][query] = query_match_len
            else:
                if not query in query_match[marker]:
                    query_match[marker][query] = query_match_len
                else:
                    query_match[marker][query] += query_match_len
            if marker not in marker_match:
                marker_match[marker] = {}
                marker_match[marker][query] = marker_match_len
            else:
                if not query in marker_match[marker]:
                    marker_match[marker][query] = marker_match_len
                else:
                    marker_match[marker][query] += marker_match_len
            if not query in query_marker_evalue:
                query_marker_evalue[query] = {}
                query_marker_evalue[query][marker] = full_evalue
            else:
                if marker in query_marker_evalue[query]:
                    if full_evalue < query_marker_evalue[query]:
                        query_marker_evalue[query][marker] = full_evalue
            if marker not in marker2queries:
                marker2queries[marker] = [query]
            else:
                marker2queries[marker].append(query)

    marker2bestquery = {} # {marker: best_query}
    for marker, queries in marker2queries.items():
        #print marker, queries
        for query in queries:
            query_len_pct = query_match[marker][query] / float(len_dict[query])
            marker_len_pct = marker_match[marker][query] / float(len_dict[marker])
            if (query_len_pct < 0.7) and (marker_len_pct < 0.7):
                continue
            else:
                if marker not in marker2bestquery:
                    marker2bestquery[marker] = query
                else:
                    old_query = marker2bestquery[marker]
                    if query_marker_evalue[query] < query_marker_evalue[old_query]:
                        marker2bestquery[marker] = query

    # query2fasta
    query2fasta = {}
    for query_rec in SeqIO.parse(query_fasta, 'fasta'):
        query = query_rec.name
        query2fasta[query] = query_rec.seq

    for marker, query in marker2bestquery.items():
        #print marker, query
        with open(os.path.join(output_folder, marker+".faa"), "w") as oh:
            name = query.rsplit('_', 1)[0]
            oh.write(">"+name+"\n")
            oh.write(str(query2fasta[query])+"\n")
