"""
download genbank genome sequences from NCBI Genome Assembly Report via ftp.
"""


import argparse
import logging
import urllib2
import ftputil
import os
import subprocess
from . import SubCommandError


_logger = logging.getLogger(__name__)


def add_arguments(parser):
    """
    Add arguments to the parser
    """
    parser.add_argument('input_file',
                                   help='either expored assembly report from NCBI Genome,\
                                    or exported ID table from NCBI Assemlby')
    parser.add_argument("-t", '--type', required=True,
                                   choices=['report', 'ID'], default='ID',
                                   help="tab-deliminated NCBI Genome assembly report file, \
                                   or NCBI Assembly ID table")
    parser.add_argument('-o', '--outdir',
                                   help="output directory, default='GenbankGenomes'",
                                   default="GenbankGenomes")
    return parser


def command(args):
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
                _logger.info("%s ==> %s"%(assembly, ftp))
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
                _logger.info("%s ==> %s"%(assembly, ftp))
                assembly2ftp[assembly] = ftp

        download_assemblies(assembly2ftp, args.outdir)


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
                _logger.info("Downloading %s ..." % \
                    remote_file)
                local_file = os.path.join(local_folder, remote_file)
                try:
                    # arguments: remote filename, local filename, callback
                    # ftputil.error.FTPOSError: Debugging info: ftputil 3.2,
                    # Python 2.7.12 (linux2)
                    host.download(remote_file, local_file)
                except Exception as e:
                    _logger.error("ERROR: CANNOT DOWNLOAD %s"%remote_file)
                    _logger.error("ERROR: MESSAGE IS: %s"%e)
                    raise SubCommandError(e)


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
            _logger.error("ERROR: can NOT download %s using address: %s"%(assembly, ftp_path))
            _logger.error("The error message is : %s"%e)
            continue
        for f in ftp_folder:
            # here convert the html content to text, and split into lines
            for line in html2text(f).split("\n"):
                line = line.strip().split()
                if len(line) <=1:
                    continue
                remote_file = line[-1]
                if remote_file.startswith(remote_folder) and remote_file.endswith(".gz"):
                    _logger.info("Downloading %s ..."%remote_file)
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
                                _logger.error("ERROR: %s was NOT downloaded after 10 tries!"%remote_file)
                                break
                    except Exception as e:
                        _logger.error("ERROR: CANNOT DOWNLOAD %s"%remote_file)
                        _logger.error("ERROR: MESSAGE IS: %s"%e)
                        raise SubCommandError(e)


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
