import gzip
from Bio import SeqIO

class SubCommandError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


commands = ['getGenomesFromGenbank', 'mergeGenbankGenomes',
            'mergeMultipleFastaRecords', 'generateShortNameTable',
            'extractGenbankGenomicFna', 'extractGenbankProteins',
            'renameFastaFile', 'AMPHORE2HmmMarkerScanner',
            'mergeEachMarkerProteinSet', 'concatMarkersforEachGenome',
            'convertMSFtoPhylip', 'runMSAonMarkers',
            'rnammerMarkerScanner', 'runSSUAlignAndMask', 'mergerRNASet',
            'runPartitionFinderProtein', 'runRAxML', 'replaceNameInTree',
            'runProdigal', 'runProkka', 'runHmmsearch'
            ]


def itermodules(root=__name__):
    for command in commands:
        yield (command, __import__('%s.%s' % (root, command), fromlist=['command']))


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
