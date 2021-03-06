#!/usr/bin/env python

# Copyright (C) 2016  Shengwei Hou
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import argparse
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from collections import Counter
from ftplib import FTP
import ftputil
import subprocess
import gzip
import urllib2


def main_usage(parser):
    """ display usage for main parser
    """
    print >>sys.stderr, parser.format_help()


def subparser_usage(argv, parser):
    """ display usage for subparser
    """
    cmd = argv[1]
    found = 0
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for choice, subparser in action.choices.items():
                if cmd == choice:
                    print >>sys.stderr, subparser.format_help()
                    found = 1
    if not found:
        if cmd in ("-v", "--version", "-h", "--help"):
            args = parser.parse_args()
        else:
            print >>sys.stderr, "\n\nERROR:%s is not a valid command!!!\n\n" % cmd
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


class SubCommandError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def html2text(strText):
    """
        This is a function to convert html file into plain text. Please credit
        the original code:
        http://stackoverflow.com/questions/14694482/converting-html-to-text-with-python
    """
    str1 = strText
    int2 = str1.lower().find("<body")
    if int2 > 0:
        str1 = str1[int2:]
    int2 = str1.lower().find("</body>")
    if int2 > 0:
        str1 = str1[:int2]
    list1 = ['<br>',  '<tr',  '<td', '</p>', 'span>', 'li>', '</h', 'div>']
    list2 = [chr(13), chr(13), chr(9), chr(13), chr(13),
             chr(13), chr(13), chr(13)]
    bolFlag1 = True
    bolFlag2 = True
    strReturn = ""
    for int1 in range(len(str1)):
        str2 = str1[int1]
        for int2 in range(len(list1)):
            if str1[int1:int1+len(list1[int2])].lower() == list1[int2]:
                strReturn = strReturn + list2[int2]
        if str1[int1:int1+7].lower() == '<script' or str1[int1:int1+9].lower() == '<noscript':
            bolFlag1 = False
        if str1[int1:int1+6].lower() == '<style':
            bolFlag1 = False
        if str1[int1:int1+7].lower() == '</style':
            bolFlag1 = True
        if str1[int1:int1+9].lower() == '</script>' or str1[int1:int1+11].lower() == '</noscript>':
            bolFlag1 = True
        if str2 == '<':
            bolFlag2 = False
        if bolFlag1 and bolFlag2 and (ord(str2) != 10):
            strReturn = strReturn + str2
        if str2 == '>':
            bolFlag2 = True
        if bolFlag1 and bolFlag2:
            strReturn = strReturn.replace(chr(32)+chr(13), chr(13))
            strReturn = strReturn.replace(chr(9)+chr(13), chr(13))
            strReturn = strReturn.replace(chr(13)+chr(32), chr(13))
            strReturn = strReturn.replace(chr(13)+chr(9), chr(13))
            strReturn = strReturn.replace(chr(13)+chr(13), chr(13))
    strReturn = strReturn.replace(chr(13), '\n')
    return strReturn


def merge_multiple_fasta_to_single_record(input_mfa, delim="N"*100):
    """ given a mfa file, read the contents of this mfa.
        If only one fasta, return this fasta sequence;
        If more than one, connect these sequences using delim then return.
    """
    # handle gzip compressed file
    if input_mfa.endswith(".gz"):
        ih = gzip.open(input_mfa, "r")
    else:
        ih = open(input_mfa, "r")

    # merge fasta records in mfa
    single_seq = ""
    fasta_records = SeqIO.parse(ih, "fasta")
    for fasta in fasta_records:
        # if single_seq already assigned, first add delim, then add seq
        if single_seq:
            single_seq += delim
            single_seq += str(fasta.seq)
        # if single_seq is not assigned, then assign fasta seq to it
        else:
            single_seq = str(fasta.seq)
    return single_seq


def mergeMultipleFastaRecords(args):
    """ given a mfa file, call merge_multiple_fasta_to_single_record,
        the
    """
    input_mfa = args.input_mfa
    delim = args.numberOfNs*'N'
    header = args.header

    with open(header+".fa", "w") as oh:
        genome_seq = merge_multiple_fasta_to_single_record(input_mfa, delim)
        oh.write(">"+header+"\n")
        oh.write(str(genome_seq)+"\n")


def ftp_download_assemblies(assembly2ftp, outdir):
    """ assembly2ftp is a dictionary contains assemblyID and it's ftp address
    """
    # login NCBI ftp site
    host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous',
                           'housw2010@gmail.com')

    # download all the files in the directory
    for assembly, ftp_path in assembly2ftp.items():

        # create a local folder
        local_folder = os.path.join(outdir, assembly)
        subprocess.check_call(['mkdir', '-p', local_folder])

        # move to the ftp_path folder on host
        server_folder = ftp_path.split("ftp.ncbi.nlm.nih.gov")[-1]
        host.chdir(server_folder)
        file_list = host.listdir(host.curdir)
        for remote_file in file_list:
            if host.path.isfile(remote_file) and \
                (not host.path.islink(remote_file)) and \
                    (remote_file.endswith('.gz')):
                print "[getGenomesFromGenbank]: Downloading %s ..." % \
                    remote_file
                local_file = os.path.join(local_folder, remote_file)
                try:
                    # arguments: remote filename, local filename, callback
                    # ftputil.error.FTPOSError: Debugging info: ftputil 3.2,
                    # Python 2.7.12 (linux2)
                    host.download(remote_file, local_file)
                except Exception as e:
                    print "[getGenomesFromGenbank]:ERROR: CANNOT DOWNLOAD %s"%remote_file
                    print "[getGenomesFromGenbank]:ERROR: MESSAGE IS: %s"%e


def download_assemblies(assembly2ftp, outdir):
    """ assembly2url is a dictionary contains assemblyID and it's NCBI ftp file address,
    this is an alternative function of ftp_download_assemblies, in case that
    doesn't work
    """
    # download all the files in the directory
    for assembly, ftp_path in assembly2ftp.items():

        # create a local folder
        local_folder = os.path.join(outdir, assembly)
        remote_folder = os.path.split(ftp_path)[-1]
        subprocess.check_call(['mkdir', '-p', local_folder])

        try:
            ftp_folder = urllib2.urlopen(ftp_path)
        except ValueError as e:
            print "[getGenomesFromGenbank]: ERROR: can NOT download %s using address: %s"%(assembly, ftp_path)
            print "[getGenomesFromGenbank]: The error message is : %s"%e
            continue
        for f in ftp_folder:
            # here convert the html content to text, and split into lines
            for line in html2text(f).split("\n"):
                line = line.strip().split()
                if len(line) <=1:
                    continue
                remote_file = line[-1]
                if remote_file.startswith(remote_folder) and remote_file.endswith(".gz"):
                    print "[getGenomesFromGenbank]: Downloading %s ..."%remote_file
                    local_file = os.path.join(local_folder, remote_file)
                    try:
                        success = False
                        tries = 0
                        cmd="wget %s -O %s"%(os.path.join(ftp_path, remote_file), local_file)
                        while not success:
                            tries += 1
                            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            stdout, stderr = p.communicate()
                            if int(p.returncode) == 0:
                                success = True
                            if tries >=10:
                                print "ERROR: %s was NOT downloaded after 10 tries!"%remote_file
                                break
                    except Exception as e:
                        print "[getGenomesFromGenbank]:ERROR: CANNOT DOWNLOAD %s"%remote_file
                        print "[getGenomesFromGenbank]:ERROR: MESSAGE IS: %s"%e





def getGenomesFromGenbank(args):
    """download genome sequences from NCBI, given an input assembly report file
    """

    input_file = args.input_file
    type = args.type

    if type == "report":
        #  parse assembly report to get assembly ID and ftp path
        assembly2ftp = {} # {assembly: ftp}
        with open(input_file, "r") as ih:
            header = ih.readline().strip().split('\t')
            assemby_idx = header.index('Assembly')
            ftp_idx = header.index('GenBank FTP')
            for line in ih:
                line = line.strip().split("\t")
                assembly = line[assemby_idx].strip()
                ftp = line[ftp_idx].strip()
                assert assembly not in assembly2ftp, \
                "Repeating genbank assembly ID %s has been found!"%assembly
                print "[getGenomesFromGenbank]: %s ==> %s"%(assembly, ftp)
                assembly2ftp.update({assembly:ftp})

        ftp_download_assemblies(assembly2ftp, args.outdir)

    else:
        # read assembly summary from genbank, store ftp path for each assembly
        ID2path = {}
        assembly_summary_genbank = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
        response = urllib2.urlopen(assembly_summary_genbank)
        for line in response:
            # consume comments
            if line.startswith("#"):
                continue
            sline = line.strip().split("\t")
            ID = sline[0].strip()
            path = sline[-1].strip()
            if path:
                ID2path[ID] = path

        # read input genbank ID file, get IDs
        assembly2ftp = {} # {genbank_ID: ftp}
        with open(input_file, "r") as ih:
            for line in ih:
                if line.startswith("GenBank Assembly ID") or line.startswith('#'):
                    continue
                line = line.strip().split("\t")
                assembly = line[0].strip()
                ftp = ID2path[assembly]
                print "[getGenomesFromGenbank]: %s ==> %s"%(assembly, ftp)
                assembly2ftp[assembly] = ftp

        download_assemblies(assembly2ftp, args.outdir)


def generateShortNameTable(args):
    """ generate a tab-deliminated lookUpTable of accession numbers and
        specie names
    """
    input_file = args.input_file
    outdir = args.outdir
    suffix = args.suffix
    prefix = args.prefix

    def _clean_name(sciName):
        """ replace characters to underscores or empty
        """
        sciName = sciName.replace(" ", "_").replace("'", "").replace("-", "_")
        sciName = sciName.replace("/", "_").replace("=", "_")
        sciName = sciName.replace("[", "").replace("]", "")
        sciName = sciName.replace("(", "_").replace(")", "_")
        sciName = sciName.replace("___", "_" ).replace("__", "_")
        return sciName

    # save all the names
    names = set()

    if prefix:
        output_file = prefix+"_shortNames.tab"
    else:
        basename = os.path.basename(input_file)
        filestem = os.path.splitext(basename)[0]
        output_file = filestem+"_shortNames.tab"

    # read in genome assembly report file
    with open(input_file, "r") as ih, open(output_file, "w") as oh:
        oh.write("#filename\tshortname\tgenus\tstrain\n")
        header = ih.readline().strip().split("\t")
        genus_pos = header.index("Organism/Name")
        strain_pos = header.index("Strain")
        accession_pos = header.index("Assembly")
        Genbank_FTP_pos = header.index("GenBank FTP")
        for line in ih:
            line = line.strip().split("\t")
            genus_full = line[genus_pos].strip().split(" ")[0].lstrip("[").rstrip("]")
            genus_short = genus_full[:5]
            accession = line[accession_pos].strip()
            assemblyID = line[Genbank_FTP_pos].strip().split("_")[-1] if "ASM" in line[-1].strip().split("_")[-1] else accession
            strain = line[strain_pos].strip().replace(" ", "")

            if strain:
                sciName = genus_short+"_"+strain+"_"+assemblyID
            else:
                sciName = genus_short+"_"+assemblyID
            sciName = _clean_name(sciName)

            accession += "."+suffix

            oh.write(accession+"\t"+sciName+"\t"+genus_full+"\t"+strain+"\n")


def mergeGenbankGenomes(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        merge fasta records for each assembly into one fasta record.
    """
    delim = args.numberOfNs*'N'
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        print "[mergeGenbankGenomes]: Merging genome sequence: "+genome_name

        with open(os.path.join(output_folder, genome_name+".mfa"), "w") as oh:
            for f in os.listdir(os.path.join(input_folder, folder)):
                if f.endswith("_genomic.fna.gz") and \
                   ("_cds_from_genomic" not in f) and \
                   ("_rna_from_genomic" not in f):
                    gz_fna = os.path.join(input_folder, folder, f)
                    genome_seq = merge_multiple_fasta_to_single_record(gz_fna, delim)

                    # write out fasta record for current folder
                    if genome_seq:
                        oh.write(">"+genome_name+"\n")
                        oh.write(str(genome_seq)+"\n")
                    else:
                        print "No genome fna was found for %s"%genome_name


def extractGenbankGenomicFna(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        extract the *_genomic.fna records from each genome subfolder, copy them to the
        specified folder.
    """
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        print "[extractGenbankGenomicFna]: Extracting genomic fna: "+genome_name

        output_file = os.path.join(output_folder, genome_name+".fna")
        for f in os.listdir(os.path.join(input_folder, folder)):
            if f.endswith("_genomic.fna.gz") and ("_from_" not in f):
                gz_fna = os.path.join(input_folder, folder, f)
                #print(gz_fna)
                try:
                    with gzip.open(gz_fna, "r") as ih, open(output_file, "w") as oh:
                        for line in ih:
                            oh.write(line)
                except Exception as e:
                    err = "[extractGenbankGenomicFna]: Error raised: %s"%e
                    raise SubCommandError(err)
                break
        else:
            print "[extractGenbankGenomicFna]: No protein file was found for %s"%genome_name
            continue


def extractGenbankProteins(args):
    """ given the genbank genome records downloaded using getGenomesFromGenbank,
        extract the protein records from each genome subfolder, copy them to the
        specified folder.
    """
    input_folder = args.input_GenbankGenomes
    output_folder = args.outdir

    subprocess.check_call(['mkdir', '-p', output_folder])

    # input folder has a lot of subfolders, each subfolder contains its genome information
    for folder in os.listdir(input_folder):

        # take the fodler name as its genome name
        genome_name = folder
        print "[extractGenbankProteins]: Extracting genome proteins: "+genome_name

        output_file = os.path.join(output_folder, genome_name+".faa")
        for f in os.listdir(os.path.join(input_folder, folder)):
            if f.endswith("_protein.faa.gz"):
                gz_faa = os.path.join(input_folder, folder, f)
                try:
                    with gzip.open(gz_faa, "r") as ih, open(output_file, "w") as oh:
                        for line in ih:
                            oh.write(line)
                except Exception as e:
                    err = "[extractGenbankProteins]: Error raised: %s"%e
                    raise SubCommandError(err)
            else:
                print "[extractGenbankProteins]: No protein file was found for %s"%genome_name
                continue


def renameFastaFile(args):
    """given an input folder with fasta files inside, rename each fasta file and
       it's header name, accroding to the lookup table
    """
    input_folder = args.inputFastaFolder
    lookupTable = args.lookupTable
    output_folder = args.outdir
    length = args.length
    subprocess.check_call(['mkdir', '-p', output_folder])


    # read in lookup table
    lookup_dict = {}
    newNames = set()
    with open(lookupTable, "r") as ih:
        for line in ih:
            if line.startswith('#'):
                continue
            line = line.strip().split("\t")
            k, v = line[0], line[1]
            assert v not in newNames, "repeat new name was not allowed: %s"%v
            assert len(v) <= length, "the length of %s was more than %d!"%(v, length)
            newNames.add(v)
            lookup_dict[k] = v

    # rename fasta records
    for fa in os.listdir(input_folder):
        if fa in lookup_dict:
            newName = lookup_dict[fa]
            print "[renameFastaFile]: renaming %s to %s"%(fa, newName)
            with open(os.path.join(output_folder, newName+".fa"), 'w') as oh:
                fa_rec = SeqIO.read(os.path.join(input_folder, fa), 'fasta')
                #fa_rec.name = newName
                #SeqIO.write(fa_rec, oh, 'fasta')
                oh.write(">"+newName+"\n")
                oh.write(str(fa_rec.seq)+"\n")
        else:
            print "[renameFastaFile]: Warining, no new name was found for %s, no change"%fa
            subprocess.check_call(['mkdir', '-p', os.path.join(output_folder, "nochange")])
            subprocess.check_call(['cp', os.path.join(input_folder, fa),
                                   os.path.join(output_folder, "nochange", fa)
                                 ])


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


def AMPHORE2HmmMarkerScanner(args):
    """ a simplified python implementation of MarkerScanner.pl of AMPHORE2, given
        a folder with fasta sequences inside, a hmm file with hmm models inside,
        search markers within input fasta sequences
    """
    fasta_folder = args.inputFastaFolder
    input_marker = args.inputMarkerFile
    force = args.force
    evalue = args.evalue
    output_folder = args.outdir
    subprocess.check_call(['mkdir', '-p', output_folder])


    for fa in os.listdir(fasta_folder):
        if not (fa.endswith('.fa') or fa.endswith('.fasta') or fa.endswith('.fas')):
            print "[AMPHORE2HmmMarkerScanner]: %s is not a fasta file, pass"%fa
            continue
        filestem = os.path.splitext(fa)[0]
        curr_outfolder = os.path.join(output_folder, filestem)
        subprocess.check_call(['mkdir', '-p', curr_outfolder])
        input_fa = os.path.join(fasta_folder, fa)
        out_orf = os.path.join(curr_outfolder, filestem+".orf")
        out_hmm = os.path.join(curr_outfolder, filestem+".hmmsearch")

        # get orf
        print "[AMPHORE2HmmMarkerScanner]: getting orfs of %s ..."%fa
        if os.path.exists(out_orf):
            if not force:
                print "[AMPHORE2HmmMarkerScanner]: orf exists, nothing to do, \
                                           use --force to override"
                continue
        subprocess.check_call(['getorf',
                                '-sequence', input_fa,
                                '-outseq', out_orf,
                                '-table', '1',
                                '-minsize', '100'
                              ])
        # run hmmsearch
        print "[AMPHORE2HmmMarkerScanner]: hmmsearch of markers for %s ..."%fa
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
        print "[AMPHORE2HmmMarkerScanner]: write marker proteins for %s ..."%fa
        get_fasta_from_hmmsearch_result(out_hmm, out_orf, curr_outfolder)



def rnammerMarkerScanner(args):
    """ a wrapper function to call rnammer for rRNA scanning, for each genome
        sequence in the input folder, it will scan the genome and found all rRNAs
    """
    rnammer = args.rnammer
    fasta_folder = args.inputFastaFolder
    kingdom = args.superkingdom
    type = args.molecule
    force = args.force
    output_folder = args.outdir
    prefix = args.prefix
    subprocess.check_call(['mkdir', '-p', output_folder])


    for fa in os.listdir(fasta_folder):
        if not (fa.endswith('.fa') or fa.endswith('.fasta') or fa.endswith('.fas')):
            print "[rnammerMarkerScanner]: %s is not a fasta file, pass"%fa
            continue
        if prefix:
            filestem = prefix
        else:
            filestem = os.path.splitext(fa)[0]
        curr_outfolder = os.path.join(output_folder, filestem)
        subprocess.check_call(['mkdir', '-p', curr_outfolder])
        input_fa = os.path.join(fasta_folder, fa)
        out_rRNA = os.path.join(curr_outfolder, filestem+".fa")
        out_gff = os.path.join(curr_outfolder, filestem+".gff")
        out_hmm = os.path.join(curr_outfolder, filestem+".hmmsearch")

        # get orf
        print "[rnammerMarkerScanner]: getting rRNAs for %s ..."%fa
        if os.path.exists(out_rRNA):
            if not force:
                print "[rnammerMarkerScanner]: rRNA file exists, nothing to do, \
                                           use --force to override"
                continue
            else:
                print "[rnammerMarkerScanner]: rRNA file for %s was override"%fa
        subprocess.check_call([rnammer,
                                '-S', kingdom,
                                '-m', type,
                                '-f', out_rRNA,
                                '-gff', out_gff,
                                '-h', out_hmm,
                                input_fa
                              ])



def mergeEachMarkerProteinSet(args):
    """ given a list of marker names, concate all found marker proteins in each
        genome.
    """
    input_folder = args.inputUmbrellaFolder
    marker_list = args.inputMarkerList
    output_folder = args.outdir
    subprocess.check_call(['mkdir', '-p', output_folder])

    # read in marker names
    markers = []
    with open(marker_list, "r") as ih:
        for line in ih:
            if line.startswith('#'):
                continue
            line = line.strip()
            markers.append(line)

    # concate marker proteins
    for marker in markers:
        with open(os.path.join(output_folder, marker+".mfa"), "w") as oh:
            for genome in os.listdir(input_folder):
                for f in os.listdir(os.path.join(input_folder, genome)):
                    if marker in f and (f.endswith(".faa") or f.endswith(".fa") or f.endswith(".fasta")):
                        fasta = SeqIO.read(os.path.join(input_folder, genome, f), 'fasta')
                        oh.write(">"+genome+"\n")
                        oh.write(str(fasta.seq)+"\n")


def mergerRNASet(args):
    """ given the name of marker rRNA, concate this marker rRNA for all genomes.
    """
    input_folder = args.inputUmbrellaFolder
    marker = args.molecule+"_rRNA"
    strategy = args.strategy       # all or best
    output_folder = args.outdir
    length_cutoff = args.length_cutoff
    numberOfNs = args.numberOfNs
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


def runSSUAlignAndMask(args):
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
        print "Error raised at running runSSUAlignAndMask:%s"%e
        sys.exit(1)

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
        err = "[runSSUAlignAndMask]: Error raised at running runSSUAlignAndMask:%s"%e
        raise SubCommandError(err)


def runMSAonMarkers(args):
    """ given an input folder contains marker protein sequences, run multiple sequence alignment
        on each marker
    """
    input_folder = args.inputMarkerProteinsFolder
    suffix = args.suffix
    if args.outdir:
        output_folder = args.outdir
    else:
        output_folder = input_folder

    for f in os.listdir(input_folder):
        if f.endswith(suffix):
            filestem = os.path.splitext(f)[0]
            # run msa
            try:
                align = subprocess.Popen([ 'muscle',
                                        '-in', os.path.join(input_folder, f),
                                        '-out', os.path.join(output_folder, filestem+".msf")
                                        ])
                align.wait()

                trim = subprocess.Popen(['trimal',
                                         '-in', os.path.join(output_folder, filestem+".msf"),
                                         '-out', os.path.join(output_folder, filestem+"_trimAl.msf"),
                                         '-gt', '0.9',
                                         '-cons', '60',
                                         '-w', '3'
                                        ])
                trim.wait()
            except Exception as e:
                err = "[runMSAonMarkers]: Error raised at runMSAonMarkers: %s"%e
                raise SubCommandError(err)

def concatMarkersforEachGenome(args):
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


def convertFasta2Phylip(input_fasta, output_prefix, format='phylip-sequential'):
    """ convert input_fasta to phylip-relaxed format using Biopython
    """
    AlignIO.convert(input_fasta, 'fasta', output_prefix+".phy", format)


def convertMSFtoPhylip(args):
    """convert multiple sequence alignments in multiple fasta format to phylip-sequential format
    """
    input_folder = args.inputSuperAlignmentFolder
    suffix = args.suffix
    phylip_relaxed = args.phylip_relaxed

    if args.outdir:
        output_folder = args.outdir
    else:
        output_folder = input_folder

    outfmt = "phylip-sequential"
    if phylip_relaxed:
        outfmt = 'phylip-relaxed'

    # convert to phylip-relaxed format
    for f in os.listdir(input_folder):
        if f.endswith(suffix):
            input_fasta = os.path.join(input_folder, f)
            basename = os.path.splitext(f)[0]
            output_prefix = os.path.join(output_folder, basename)
            print "[convertMSFtoPhylip]: converting %s to %s format ..."%(f, outfmt)

            # it's compulsory to use phylip-sequential format for PartitionFinderProtein
            convertFasta2Phylip(input_fasta, output_prefix, format=outfmt)


def runPartitionFinderProtein(args):
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
        raise SubCommandError(err)


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
        print  " ".join(command)
        p = subprocess.Popen(command)
        p.wait()
    except Exception as e:
        err = "[runRAxML_rapid_bootstrap]: Error raised at runRAxML_rapid_bootstrap:%s"%e
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
        raise SubCommandError(err)

    # step2: 100 bootstrap search
    #/home/hou/Software/RAxML/standard-RAxML/raxmlHPC-PTHREADS -m PROTGAMMAILGF -p 12345 -b 12345 -# 100 -s tcoffee_mafft_fileset_trimAl.phy -T 21 -n PROTGAMMAILGF_100BS

    # step3: projection, draw bipartitions on the best ML tree using bootstrap replicate trees
    # for unrooted trees the correct representation is actually the one with support values assigned to branches and not nodes
    #/home/hou/Software/RAxML/standard-RAxML/raxmlHPC-PTHREADS -m PROTGAMMAILGF -p 12345 -f b -t RAxML_bestTree.PROTGAMMAILGF_20ML -z RAxML_bootstrap.PROTGAMMAILGF_100BS -n PROTGAMMAILGF_projection


def runRAxML(args):
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
        raise SubCommandError(err)
        runRAxML_slow_bootstrap(raxml, alignment, model, bootstraps, outgroup,
                                 suffix, structure, partition, threads)


def replaceNameInTree(args):
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
                print "[replaceNameInTree]: WARNING: no scientificname were found here: %s"%("\t".join(line))
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
                        print "[replaceNameInTree]: WARNING: %s has been used, and should be unique in tree!"%short
                        continue
                oh.write(line)
        else:
            for line in ih:
                for short, sci in short2sci.iteritems():
                    if short in line:
                        line = line.replace(short, sci)
                oh.write(line)

        print "[replaceNameInTree]: Please check the names and manually correct wrong names!!! "


def construct_prodigal_cmds(prodigal, input_fasta, output_prefix,
                            output_fmt='gff', write_proteins=None, write_DNAs=None,
                            m=None, c=None, n=None, *args, **kwargs):
    """ construct prodigal commands
    """
    template = '%s -i %s -o %s'%(prodigal, input_fasta, output_prefix+"."+output_fmt)

    command = template.split(" ")
    if write_proteins:
        command.extend(['-a', output_prefix+".faa"])
    if write_DNAs:
        command.extend(['-d', output_prefix+".fsa"])
    if m:
        command.extend(['-m'])
    if c:
        command.extend(['-c'])
    if n:
        command.extend(['-c'])

    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "[runProdigal]: Error raised when running prodigal:%s"%e
        raise SubCommandError(err)


def runProdigal(args):
    """ wrapper function of prodigal
    """
    prodigal = args.prodigal
    input_file_or_folder = args.input_file_or_folder
    suffix = args.suffix

    output_fmt = args.format
    m = args.m
    c = args.c
    n = args.n

    write_proteins = not args.not_write_proteins
    write_DNAs = not args.not_write_DNAs

    # create outdir
    subprocess.check_call(['mkdir', '-p', args.outdir])

    if os.path.isdir(input_file_or_folder):
        for f in os.listdir(input_file_or_folder):
            if f.endswith(suffix):
                if args.prefix:
                    prefix = os.path.join(args.outdir, args.prefix)
                else:
                    prefix = os.path.join(args.outdir, os.path.splitext(f)[0])

                # run prodigal
                input_fasta = os.path.join(input_file_or_folder, f)
                construct_prodigal_cmds(prodigal, input_fasta, prefix, output_fmt,
                                        write_proteins, write_DNAs, m, c, n)

    else:
        if args.prefix:
            prefix = os.path.join(args.outdir, args.prefix)
        else:
            basename = os.path.basename(input_file_or_folder)
            prefix = os.path.join(args.outdir, os.path.splitext(basename)[0])

        # run prodigal
        input_fasta = input_file_or_folder
        construct_prodigal_cmds(prodigal, input_fasta, prefix, output_fmt,
                                write_proteins, write_DNAs, m, c, n)


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
        if ignoreError:
            print err
        else:
            raise SubCommandError(err)


def runProkka(args):
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
                    print("[runProkka]: found sequence %s more than once!"%filename)
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

def construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, output_prefix,
                            evalue=1e-7, dom_evalue=1e-7,
                            write_tblout=None, write_domtblout=None,
                            write_pfamtblout=None,
                            acc=None, noali=None, cpu=10, *args, **kwargs):
    """ construct hmmsearch commands
    """
    # hmmsearch -o /dev/null -Z 5000 -E 1e-7 --domE 1e-7 --domtblout output_dom_file \
    # --tblout output_seq_file --pfamtblout output_pfam_file hmm_file input_fasta
    template = '%s -o /dev/null -Z 5000 -E %s --domE %s --cpu %s'%(
        hmmsearch, str(evalue), str(dom_evalue), str(cpu))

    command = template.split(" ")
    if write_tblout:
        command.extend(['--tblout', output_prefix+".tblout"])
    if write_domtblout:
        command.extend(['--domtblout', output_prefix+".domtblout"])
    if write_pfamtblout:
        command.extend(['--pfamtblout', output_prefix+".pfamtblout"])
    if acc:
        command.extend(['--acc'])
    if noali:
        command.extend(['--noali'])

    command.extend([hmm_file, input_fasta])

    try:
        p1 = subprocess.Popen(command)
        p1.wait()
    except Exception as e:
        err = "[runHmmsearch]: Error raised when running hmmsearch:%s"%e
        raise SubCommandError(err)


def runHmmsearch(args):
    """wrapper function of hmmsearch
    """
    hmmsearch = args.hmmsearch
    input_file_or_folder = args.input_file_or_folder
    hmm_file = args.input_hmm_file
    outdir = args.outdir
    suffix = args.suffix
    evalue = args.evalue
    dom_evalue = args.dom_evalue
    write_tblout = not args.not_write_tblout
    write_domtblout = not args.not_write_domtblout
    write_pfamtblout = not args.not_write_pfamtblout
    acc = args.acc
    noali = args.noali
    cpu = args.cpu

    # create outdir
    subprocess.check_call(['mkdir', '-p', outdir])

    if os.path.isdir(input_file_or_folder):
        for f in os.listdir(input_file_or_folder):
            if f.endswith(suffix):
                if args.prefix:
                    prefix = os.path.join(args.outdir, args.prefix)
                else:
                    prefix = os.path.join(args.outdir, os.path.splitext(f)[0])

                # run hmmsearch
                input_fasta = os.path.join(input_file_or_folder, f)
                construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, prefix,
                                         evalue, dom_evalue, write_tblout,
                                         write_domtblout, write_pfamtblout,
                                         acc, noali, cpu)

    else:
        if args.prefix:
            prefix = os.path.join(args.outdir, args.prefix)
        else:
            basename = os.path.basename(input_file_or_folder)
            prefix = os.path.join(args.outdir, os.path.splitext(basename)[0])

        # run hmmsearch
        input_fasta = input_file_or_folder
        construct_hmmsearch_cmds(hmmsearch, hmm_file, input_fasta, prefix,
                                 evalue, dom_evalue, write_tblout,
                                 write_domtblout, write_pfamtblout,
                                 acc, noali, cpu)




def main():

    # main parser
    parser = argparse.ArgumentParser(description="A set of subcommands for phylogenetic computation")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.0.1")


    # parent parser, to specify shared arguments, inherited by subparsers
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-p", "--prefix", required=False, help="output prefix")



    # subparsers
    subparsers = parser.add_subparsers(help='available subcommands')

    # --------------------- #
    # getGenomesFromGenbank #
    # --------------------- #
    getGenomesFromGenbank_desc="ftp download genbank genome sequences from \
                                   NCBI Genome Assembly Report."
    parser_getGenomesFromGenbank = subparsers.add_parser('getGenomesFromGenbank',
                                   parents=[parent_parser],
                                   help=getGenomesFromGenbank_desc,
                                   description=getGenomesFromGenbank_desc)
    parser_getGenomesFromGenbank.add_argument('input_file',
                                   help='either expored assembly report from NCBI Genome,\
                                    or exported ID table from NCBI Assemlby')
    parser_getGenomesFromGenbank.add_argument("-t", '--type', required=True,
                                   choices=['report', 'ID'], default='ID',
                                   help="tab-deliminated NCBI Genome assembly report file, \
                                   or NCBI Assembly ID table")
    parser_getGenomesFromGenbank.add_argument('-o', '--outdir',
                                   help="output directory, default='GenbankGenomes'",
                                   default="GenbankGenomes")
    parser_getGenomesFromGenbank.set_defaults(func=getGenomesFromGenbank)


    # ---------------------- #
    # generateShortNameTable #
    # ---------------------- #
    generateShortNameTable_desc="generate tab-deliminated shortName table from \
                                   NCBI Genome Assembly Report, manually corrections \
                                   might be needed."
    parser_generateShortNameTable = subparsers.add_parser('generateShortNameTable',
                                   parents=[parent_parser],
                                   help=generateShortNameTable_desc,
                                   description=generateShortNameTable_desc)
    parser_generateShortNameTable.add_argument('input_file',
                                   help='expored assembly report from NCBI Genome')
    parser_generateShortNameTable.add_argument('-o', '--outdir',
                                   help="output directory, default='./'",
                                   default="./")
    parser_generateShortNameTable.add_argument('-s', '--suffix',
                                   help="suffix will be appended to the accession\
                                   number, default='.mfa'", default='.mfa')
    parser_generateShortNameTable.set_defaults(func=generateShortNameTable)


    # ------------------- #
    # mergeGenbankGenomes #
    # ------------------- #
    mergeGenbankGenomes_desc="merge contigs/scaffolds/plasmids for each folder \
                                 in GenbankGenomes into one fasta."
    parser_mergeGenbankGenomes = subparsers.add_parser('mergeGenbankGenomes',
                                 parents=[parent_parser],
                                 help=mergeGenbankGenomes_desc,
                                 description=mergeGenbankGenomes_desc)
    parser_mergeGenbankGenomes.add_argument('input_GenbankGenomes', help="\
                                 input directory path to GenbankGenomes")
    parser_mergeGenbankGenomes.add_argument('-n', '--numberOfNs', type=int, default=100,
                                 help="number of Ns used to seperate each fasta \
                                 record, default=100")
    parser_mergeGenbankGenomes.add_argument('-o', '--outdir',
                                   help="output directory, default='MergedGenbankGenomes'",
                                   default="MergedGenbankGenomes")
    parser_mergeGenbankGenomes.set_defaults(func=mergeGenbankGenomes)


    # ------------------------ #
    # extractGenbankGenomicFna #
    # ------------------------ #
    extractGenbankGenomicFna_desc="extract *_genomic.fna.gz file from each genome subfolder \
                                 of downloaded GenbankGenomes into one folder."
    parser_extractGenbankGenomicFna = subparsers.add_parser('extractGenbankGenomicFna',
                                 parents=[parent_parser],
                                 help=extractGenbankGenomicFna_desc,
                                 description=extractGenbankGenomicFna_desc)
    parser_extractGenbankGenomicFna.add_argument('input_GenbankGenomes', help="\
                                 input directory path to GenbankGenomes")
    parser_extractGenbankGenomicFna.add_argument('-o', '--outdir',
                                   help="output directory, default='extractedGenbankGenomicFna'",
                                   default="extractedGenbankGenomicFna")
    parser_extractGenbankGenomicFna.set_defaults(func=extractGenbankGenomicFna)


    # ---------------------- #
    # extractGenbankProteins #
    # ---------------------- #
    extractGenbankProteins_desc="extract *_protein.faa.gz file from each genome subfolder \
                                 of downloaded GenbankGenomes into one folder."
    parser_extractGenbankProteins = subparsers.add_parser('extractGenbankProteins',
                                 parents=[parent_parser],
                                 help=extractGenbankProteins_desc,
                                 description=extractGenbankProteins_desc)
    parser_extractGenbankProteins.add_argument('input_GenbankGenomes', help="\
                                 input directory path to GenbankGenomes")
    parser_extractGenbankProteins.add_argument('-o', '--outdir',
                                   help="output directory, default='extractedGenbankProteins'",
                                   default="extractedGenbankProteins")
    parser_extractGenbankProteins.set_defaults(func=extractGenbankProteins)

    # ------------------------- #
    # mergeMultipleFastaRecords #
    # ------------------------- #
    mergeMultipleFastaRecords_desc="merge contigs/scaffolds/plasmids for one mfa"

    parser_mergeMultipleFastaRecords = subparsers.add_parser('mergeMultipleFastaRecords',
                                 parents=[parent_parser],
                                 help=mergeMultipleFastaRecords_desc,
                                 description=mergeMultipleFastaRecords_desc)
    parser_mergeMultipleFastaRecords.add_argument('input_mfa', help="\
                                 input mfa file, contains multiple fasta records")
    parser_mergeMultipleFastaRecords.add_argument('header', help='header for \
                                 merged contigs/scaffolds/plasmids')
    parser_mergeMultipleFastaRecords.add_argument('-n', '--numberOfNs', type=int, default=100,
                                 help="number of Ns used to seperate each fasta \
                                 record, default=100")
    parser_mergeMultipleFastaRecords.set_defaults(func=mergeMultipleFastaRecords)

    # --------------- #
    # renameFastaFile #
    # --------------- #
    renameFastaFile_desc="change all the filenames and record names of fasta files,\
                                   that inside of given folder."
    parser_renameFastaFile = subparsers.add_parser('renameFastaFile',
                                   parents=[parent_parser],
                                   help=renameFastaFile_desc,
                                   description=renameFastaFile_desc)
    parser_renameFastaFile.add_argument("inputFastaFolder",
                                   help="input folder that contains all fasta files")
    parser_renameFastaFile.add_argument("lookupTable",
                                   help='tab-deliminated lookup table to convert \
                                   old names to new names, in "old tab new" format,\
                                   if more than two columns were detected in the \
                                   lookup table, only the first two columns will be used.')
    parser_renameFastaFile.add_argument('-o', '--outdir',
                                   help="output directory, default='renamedFasta'",
                                   default="renamedFasta")
    parser_renameFastaFile.add_argument('-l', '--length', type=int, help='maximum \
                                   length of allowed short names in the lookup table. \
                                   Default=9.', default=9)
    parser_renameFastaFile.set_defaults(func=renameFastaFile)


    # ------------------------ #
    # AMPHORE2HmmMarkerScanner #
    # ------------------------ #
    AMPHORE2HmmMarkerScanner_desc="scan markers using hmmsearch, for each fasta file in \
                                   given folder, a python implementation of AMPHORE2."
    parser_AMPHORE2HmmMarkerScanner = subparsers.add_parser('AMPHORE2HmmMarkerScanner',
                                   parents=[parent_parser],
                                   help=AMPHORE2HmmMarkerScanner_desc,
                                   description=AMPHORE2HmmMarkerScanner_desc)
    parser_AMPHORE2HmmMarkerScanner.add_argument("inputFastaFolder",
                                   help="input folder that contains all genome \
                                   fasta files")
    parser_AMPHORE2HmmMarkerScanner.add_argument("inputMarkerFile",
                                   help='input marker file that contains hmm models')
    parser_AMPHORE2HmmMarkerScanner.add_argument('-e', '--evalue', default=1e-7,
                                   help="evalue cutoff, default=1e-7")
    parser_AMPHORE2HmmMarkerScanner.add_argument('-o', '--outdir',
                                   help="output directory, default='ScannedMarkers'",
                                   default="ScannedMarkers")
    parser_AMPHORE2HmmMarkerScanner.add_argument('-f', '--force',
                                   help="force to override previous result?")
    parser_AMPHORE2HmmMarkerScanner.set_defaults(func=AMPHORE2HmmMarkerScanner)


    # -------------------- #
    # rnammerMarkerScanner #
    # -------------------- #
    rnammerMarkerScanner_desc="scan rRNA markers using rnammer, for each fasta \
                                   file in given folder."
    parser_rnammerMarkerScanner = subparsers.add_parser('rnammerMarkerScanner',
                                   parents = [parent_parser],
                                   help = rnammerMarkerScanner_desc,
                                   description = rnammerMarkerScanner_desc)
    parser_rnammerMarkerScanner.add_argument("inputFastaFolder",
                                   help = "input folder that contains all genome \
                                   fasta files")
    parser_rnammerMarkerScanner.add_argument('-S', '--superkingdom', help="superkingdom, choose arc/bac/euk, \
                                   default='bac'",
                                   choices=['arc', 'bac', 'euk'], default='bac')
    parser_rnammerMarkerScanner.add_argument('-m', '--molecule', help='type of rRNAs, can be any \
                                   combination of the three, \
                                   default is "tsu,lsu,ssu",', default='tsu,lsu,ssu',
                                   choices=['<tsu,lsu,ssu>', '<tsu,lsu>', '<tsu,ssu>', '<lsu,ssu>',
                                   '<tsu>', '<lsu>', '<ssu>'])
    parser_rnammerMarkerScanner.add_argument('-o', '--outdir', help='output directory, \
                                   default=ScannedrRNAs', default="ScannedrRNAs")
    parser_rnammerMarkerScanner.add_argument('-f', '--force',
                                   help="force to override previous result?")
    parser_rnammerMarkerScanner.add_argument('-r', '--rnammer', help='path to rnammer excutable, \
                                   default=~/Software/Module_Annotation/rnammer/rnammer-1.2/rnammer',
                                   default='/home/hou/Software/Module_Annotation/rnammer/rnammer-1.2/rnammer')
    parser_rnammerMarkerScanner.set_defaults(func=rnammerMarkerScanner)


    # ------------------------- #
    # mergeEachMarkerProteinSet #
    # ------------------------- #
    mergeEachMarkerProteinSet_desc="concatenate marker proteins found in each genome, \
                                   given an umbrella folder and proteins in genome subfolders."
    parser_mergeEachMarkerProteinSet = subparsers.add_parser('mergeEachMarkerProteinSet',
                                   parents=[parent_parser],
                                   help=mergeEachMarkerProteinSet_desc,
                                   description=mergeEachMarkerProteinSet_desc)
    parser_mergeEachMarkerProteinSet.add_argument("inputUmbrellaFolder",
                                   help="input folder that contains all genome subfolders,\
                                   where extracted marker proteins were inside")
    parser_mergeEachMarkerProteinSet.add_argument("inputMarkerList",
                                   help='input marker list file that contains all marker names')
    parser_mergeEachMarkerProteinSet.add_argument('-o', '--outdir',
                                   help="output directory, default='MergedMarkerProteinSet'",
                                   default="MergedMarkerProteinSet")
    parser_mergeEachMarkerProteinSet.set_defaults(func=mergeEachMarkerProteinSet)


    # ------------ #
    # mergerRNASet #
    # ------------ #
    mergerRNASet_desc="merge marker rRNAs found in each genome into one file,\
                                   given an umbrella folder and rRNA sequences in genome subfolders."
    parser_mergerRNASet = subparsers.add_parser('mergerRNASet',
                                   parents=[parent_parser],
                                   help=mergerRNASet_desc,
                                   description=mergerRNASet_desc)
    parser_mergerRNASet.add_argument("inputUmbrellaFolder",
                                   help="input folder that contains all genome subfolders,\
                                   where extracted marker rRNAs were inside")
    parser_mergerRNASet.add_argument("-m", '--molecule', default='16s',
                                   choices=['5s', '16s', '23s', '8s', '18s', '28s'],
                                   help='molecule type you want to merge, default=16s')
    parser_mergerRNASet.add_argument('-s', '--strategy', default='all',
                                   choices = ['all', 'best'], help = 'strategy to select\
                                   marker rRNA when multiple markers were found, default=all')
    parser_mergerRNASet.add_argument('-l', '--length_cutoff', default=1200, type=int,
                                   help='length cutoff, to exclude short sequences, default=1200')
    parser_mergerRNASet.add_argument('-n', '--numberOfNs', default=50, type=int,
                                   help='maximum Ns allowed in the sequences, default=50')
    parser_mergerRNASet.add_argument('-o', '--outdir',
                                   help="output directory, default='MergedrRNASet'",
                                   default="MergedrRNASet")
    parser_mergerRNASet.set_defaults(func=mergerRNASet)


    # ------------------ #
    # runSSUAlignAndMask #
    # ------------------ #
    runSSUAlignAndMask_desc = "run ssu-align and mask for input rRNA sequences"
    parser_runSSUAlignAndMask = subparsers.add_parser('runSSUAlignAndMask',
                                   parents=[parent_parser],
                                   help=runSSUAlignAndMask_desc,
                                   description=runSSUAlignAndMask_desc)
    parser_runSSUAlignAndMask.add_argument('inputFastaFile', help='input fasta\
                                   file contains multiple rRNA sequences')
    parser_runSSUAlignAndMask.add_argument('-f', '--force', action='store_true',
                                   help='override outputs')
    parser_runSSUAlignAndMask.add_argument('--no_align', action='store_true',
                                   help='only search target \
                                   sequence file for hits, skip alignment step')
    parser_runSSUAlignAndMask.add_argument('--no_search', action='store_true',
                                   help='only align  target \
                                   sequence file, skip initial search step')
    parser_runSSUAlignAndMask.add_argument('--output_RNA', action='store_true',
                                   help='output alignment as RNA, default is DNA')
    parser_runSSUAlignAndMask.add_argument('--output_stk', action='store_true',
                                   help='output stockholm alignments, \
                                   default is fasta alignments')
    parser_runSSUAlignAndMask.add_argument('--kingdom', default='bacteria',
                                   choices=['archaea', 'bacteria', 'eukarya'],
                                   help='kingdom to search and align to, default is bacteria')
    parser_runSSUAlignAndMask.add_argument('--ssu_align', help='path to ssu-align excutable, \
                                   default=~/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-align',
                                   default='/home/hou/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-align')
    parser_runSSUAlignAndMask.add_argument('--ssu_mask', help='path to ssu-mask excutable, \
                                   default=~/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-mask',
                                   default='/home/hou/Software/Module_Phylogeny/ssu-align-0.1.1/src/ssu-mask')
    parser_runSSUAlignAndMask.set_defaults(func=runSSUAlignAndMask)



    # --------------- #
    # runMSAonMarkers #
    # --------------- #
    runMSAonMarkers_desc="given a folder with marker proteins inside, \
                          run multiple resquence alignment on each marker."
    parser_runMSAonMarkers = subparsers.add_parser('runMSAonMarkers',
                                   parents = [parent_parser],
                                   help=runMSAonMarkers_desc,
                                   description=runMSAonMarkers_desc)
    parser_runMSAonMarkers.add_argument('inputMarkerProteinsFolder',
                                   help="input folder contains identified marker proteins")
    parser_runMSAonMarkers.add_argument('-o', '--outdir',
                                   help='output directory, default is input directory')
    parser_runMSAonMarkers.add_argument('-s', '--suffix', default='.mfa',
                                   help='suffix to identify input protein fasta file, default=.mfa')
    parser_runMSAonMarkers.set_defaults(func=runMSAonMarkers)


    # -------------------------- #
    # concatMarkersforEachGenome #
    # -------------------------- #
    concatMarkersforEachGenome_desc="given a folder with trimmed marker protein alignment inside, \
                          concatenate markers for each genome, gaps will be placed if no marker found for one genome."
    parser_concatMarkersforEachGenome = subparsers.add_parser('concatMarkersforEachGenome',
                                   parents = [parent_parser],
                                   help=concatMarkersforEachGenome_desc,
                                   description=concatMarkersforEachGenome_desc)
    parser_concatMarkersforEachGenome.add_argument('inputMarkerAlignmentFolder',
                                   help="input folder contains marker protein alignments")
    parser_concatMarkersforEachGenome.add_argument('-o', '--outdir',
                                   help='output directory, default is input directory')
    parser_concatMarkersforEachGenome.add_argument('-s', '--suffix', default='_trimAl.msf',
                                   help='suffix to identify input trimmed protein alignment file, default=_trimAl.msf')
    parser_concatMarkersforEachGenome.set_defaults(func=concatMarkersforEachGenome)


    # ------------------ #
    # convertMSFtoPhylip #
    # ------------------ #
    convertMSFtoPhylip_desc="convert multiple sequence alignments in multiple fasta format \
                             to phylip-relaxed format."
    parser_convertMSFtoPhylip = subparsers.add_parser('convertMSFtoPhylip',
                                   parents = [parent_parser],
                                   help=convertMSFtoPhylip_desc,
                                   description=convertMSFtoPhylip_desc)
    parser_convertMSFtoPhylip.add_argument('inputSuperAlignmentFolder',
                                   help="input folder contains supermatrix alignments")
    parser_convertMSFtoPhylip.add_argument('-o', '--outdir',
                                   help='output directory, default is input directory')
    parser_convertMSFtoPhylip.add_argument('--phylip_relaxed', action='store_true',
                                   help='output phylip-relaxed\
                                   format, default is phylip-sequential')
    parser_convertMSFtoPhylip.add_argument('-s', '--suffix', default='_supermatrix.msf',
                                   help='suffix to identify input super alignment file, default=_supermatrix.msf')
    parser_convertMSFtoPhylip.set_defaults(func=convertMSFtoPhylip)


    # ------------------------- #
    # runPartitionFinderProtein #
    # ------------------------- #
    runPartitionFinderProtein_desc="estimate models for partitions using PartitionFinderProtein."
    PartitionFinderProtein = '/home/hou/Software/Module_Phylogeny/PartitionFinder/PartitionFinderV1.1.1/PartitionFinderProtein.py'
    parser_runPartitionFinderProtein = subparsers.add_parser('runPartitionFinderProtein',
                                   parents = [parent_parser],
                                   help=runPartitionFinderProtein_desc,
                                   description=runPartitionFinderProtein_desc)
    parser_runPartitionFinderProtein.add_argument('inputPhylipSuperAlign',
                                   help="input supermatrix alignment in phylip-sequential format")
    parser_runPartitionFinderProtein.add_argument('inputCoordFile',
                                   help='input coordinate file contains the positions of proteins')
    parser_runPartitionFinderProtein.add_argument('newWorkingFolder',
                                   help='new working directory, to put alignment and config files')
    parser_runPartitionFinderProtein.add_argument('--PartitionFinderProtein', help="path to\
                                   PartitionFinderProtein excutable, default=%s"%PartitionFinderProtein,
                                   default=PartitionFinderProtein)
    parser_runPartitionFinderProtein.set_defaults(func=runPartitionFinderProtein)


    # -------- #
    # runRAxML #
    # -------- #
    parser_runRAxML = subparsers.add_parser('runRAxML', parents=[parent_parser],
                                   help='wrapper for running RAxML, with or without partitions',
                                   description='wrapper for running RAxML, with or without partitions')
    parser_runRAxML.add_argument('input_alignment', help='input alignment in phylip format')
    parser_runRAxML.add_argument('-r', '--raxml', help='path to raxml excutable, \
                                   default=~/Software/Module_Phylogeny/RAxML/standard-RAxML/raxmlHPC',
                                   default='/home/hou/Software/Module_Phylogeny/RAxML/standard-RAxML/raxmlHPC')
    parser_runRAxML.add_argument("--isDNA", action="store_true",
                               help="are DNA sequence alignment? default is Protein")
    parser_runRAxML.add_argument('-m', '--model', help='model will be directly passed to RAxML, \
                                   default is GTRGAMMA for DNA, PROTGAMMAWAG for Protein.')
    parser_runRAxML.add_argument('-#', '--bootstraps', type=int, help='number of bootstraps, default=100',
                                   default=100)
    parser_runRAxML.add_argument('-o', '--outgroup', type=str,
                                   help='sequence names of outgroups, comma-delimited')
    parser_runRAxML.add_argument('--slow', action='store_true',
                                   help = 'do slow bootstrapping instead? default is rapid bootstrapping')
    parser_runRAxML.add_argument('-S', '--Structure', help='path to secondary structure file')
    parser_runRAxML.add_argument('-q', '--Partition', help='path to partition file')
    parser_runRAxML.add_argument('--suffix', help='suffix string')
    parser_runRAxML.add_argument('--threads', type=int, default=30,
                                   help='number of threads, default=30')
    parser_runRAxML.set_defaults(func=runRAxML)

    # ----------------- #
    # replaceNameInTree #
    # ----------------- #
    replaceNameInTree_desc="replace shortname back to scientific names in input tree file,\
    it will only replace the shortname strings found the tree file, should not change \
    anything else like topology/branchlength/bootstraps of the tree."
    parser_replaceNameInTree = subparsers.add_parser('replaceNameInTree',
                                   parents = [parent_parser],
                                   help=replaceNameInTree_desc,
                                   description=replaceNameInTree_desc)
    parser_replaceNameInTree.add_argument('inputTreeFile',
                                   help="input tree file generated using RAxML, in newick/nexus format")
    parser_replaceNameInTree.add_argument('-f', "--format", choices=['newick', 'nexus'],
                                   help="input tree format, default=newick", default='newick')
    parser_replaceNameInTree.add_argument('lookUpTable',
                                   help='tab-deliminated lookup table to convert \
                                   short names to scientific names, this is the same file \
                                   used in renameFastaFile, the second column should \
                                   be shortname and third column should be scientificname, \
                                   no comma/space/(/)/colon were allowed in scientificname.')

    parser_replaceNameInTree.set_defaults(func=replaceNameInTree)


    # ----------- #
    # runProdigal #
    # ----------- #
    parser_runProdigal = subparsers.add_parser('runProdigal', parents=[parent_parser],
                                   help='wrapper for running prodigal',
                                   description='wrapper for running runProdigal')
    parser_runProdigal.add_argument('input_file_or_folder', help='input file or folder')
    parser_runProdigal.add_argument('-o', '--outdir', help='output directory, default=\
                                    ./prodigalOutput', default='prodigalOutput')
    parser_runProdigal.add_argument('-a', '--not_write_proteins', action='store_true',
                                    help='do not write predicted proteins?')
    parser_runProdigal.add_argument('-d', '--not_write_DNAs', action='store_true',
                                    help='do not write nucleotide sequences of predicted proteins?')
    parser_runProdigal.add_argument('-f', '--format', choices=['gbk', 'gff'], help='output format, \
                                   for gene annotations, can be gff or gbk, default is gff',
                                   default='gff')
    parser_runProdigal.add_argument('-g', type=int, help='specify a translation table to use, default=11',
                                    default=11)
    parser_runProdigal.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                                    with this suffix will be processed. default=.fa',
                                    default='.fa')
    parser_runProdigal.add_argument('-m', action='store_true', help='Treat runs of N as masked sequence; \
                                    do not build genes across them.')
    parser_runProdigal.add_argument('-c', action='store_true', help='Closed ends. \
                                    Do not allow genes to run off edges.')
    parser_runProdigal.add_argument('-n', action='store_true', help='Bypass Shine-Dalgarno\
                                    trainer and force a full motif scan.')
    parser_runProdigal.add_argument('-r', '--prodigal', help='path to prodigal excutable, \
                                   default=~/Software/Module_Annotation/prodigal/prodigal',
                                   default='/home/hou/Software/Module_Annotation/prodigal/prodigal')
    parser_runProdigal.set_defaults(func=runProdigal)

    # --------- #
    # runProkka #
    # --------- #
    parser_runProkka = subparsers.add_parser('runProkka', parents=[parent_parser],
                                   help='wrapper for running prokka',
                                   description='wrapper for running prokka')
    parser_runProkka.add_argument('input_file_or_folder', help='input file or folder')
    parser_runProkka.add_argument('-o', '--outdir', help='output directory, default=\
                                    ./prokkaOutput', default='prokkaOutput')
    parser_runProkka.add_argument('-a', '--parameter_file', help="tab-deliminated "\
                                  "parameter file, contains fileanme, genus and strain")
    parser_runProkka.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                                    with this suffix will be processed. default=.fa',
                                    default='.fa')
    parser_runProkka.add_argument('-c', '--cpu', type=int,
                                    help='Number of CPUs to use, [0=all] default=8.',
                                    default=8)
    parser_runProkka.add_argument('--nogenes', action='store_true',
                                    help="Do not add 'gene' features for each 'CDS' feature")
    parser_runProkka.add_argument('--donotforce', action='store_true',
                                    help="Do not force remove existing folders")
    parser_runProkka.add_argument('--addmrna', action='store_true',
                                    help="Add 'mRNA' features for each 'CDS' feature (default OFF)")
    parser_runProkka.add_argument('--locustag',
                                    help="Locus tag prefix (default 'PROKKA' or filestem of \
                                    input files when input is a folder).",
                                    default='PROKKA')
    parser_runProkka.add_argument('--increment', type = int,
                                    help="Locus tag counter increment (default '1').",
                                    default=1)
    parser_runProkka.add_argument('--kingdom',
                                    help="Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')",
                                    default='Bacteria')
    parser_runProkka.add_argument('--ignoreError', action="store_true",
                                    help="ignore Errors raised when annotating multiple genomes?")
    parser_runProkka.add_argument('-r', '--prokka', help='path to prokka excutable, \
                                   default=~/Software/Module_Annotation/prokka/prokka-1.12/bin/prokka',
                                   default='/home/hou/Software/Module_Annotation/prokka/prokka-1.12/bin/prokka')
    parser_runProkka.set_defaults(func=runProkka)

    # ------------ #
    # runHmmsearch #
    # ------------ #
    parser_runHmmsearch = subparsers.add_parser('runHmmsearch', parents=[parent_parser],
                                   help='wrapper for running hmmsearch',
                                   description='wrapper for running runHmmsearch')
    parser_runHmmsearch.add_argument('input_hmm_file', help='input hmm file')
    parser_runHmmsearch.add_argument('input_file_or_folder', help='input file or folder')
    parser_runHmmsearch.add_argument('-o', '--outdir', help='output directory, default=\
                                    ./hmmsearchOutput', default='hmmsearchOutput')
    parser_runHmmsearch.add_argument('-s', '--suffix', help='if the input is a folder, only files\
                                    with this suffix will be processed. default=.faa',
                                    default='.faa')
    parser_runHmmsearch.add_argument('-e', '--evalue', type=float,
                                     help='maximum sequence E-value threshold to \
                                     report a hit in the output, default=1e-7',
                                     default=1e-7)
    parser_runHmmsearch.add_argument('-m', '--dom_evalue', type=float,
                                     help='maximum domain E-value threshold to \
                                     report a hit in the output, default=1e-7',
                                     default=1e-7)
    parser_runHmmsearch.add_argument('-t', '--not_write_tblout', action='store_true',
                                     help='do not write table of per-sequence hits to file')
    parser_runHmmsearch.add_argument('-d', '--not_write_domtblout', action='store_true',
                                     help='do not write table of per-domain hits to file')
    parser_runHmmsearch.add_argument('-f', '--not_write_pfamtblout', action='store_true',
                                     help=' do not write table of hits and domains to file \
                                     in Pfam format')
    parser_runHmmsearch.add_argument('--acc', action='store_true', help='prefer accessions\
                                     over names in output')
    parser_runHmmsearch.add_argument('--noali', action='store_true', help='do not output \
                                     alignments, so output is smaller')
    parser_runHmmsearch.add_argument('--cpu', type=int, help='number of CPU to use, \
                                      default=10', default=10)
    parser_runHmmsearch.add_argument('-r', '--hmmsearch', help='path to hmmsearch excutable, \
                                   default=~/Software/Module_Annotation/hmmer/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch',
                                   default='/home/hou/Software/Module_Annotation/hmmer/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch')
    parser_runHmmsearch.set_defaults(func=runHmmsearch)


    # ----------------------- #
    # parse arguments and run #
    # ----------------------- #
    # display help
    display_help(sys.argv, parser)

    # parse args
    args = parser.parse_args()

    # run commands
    args.func(args)

if __name__ == "__main__":
    sys.exit(main())
