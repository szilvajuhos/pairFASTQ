import sys
import re
import click

# For general desription see the method assigned as  __main__ below.

def getReads(fh,chunkSize,reads):
    """Read a chunk of reads and put it into the dictionary
    """
    items = 0
    haveMoreReads = True
    while chunkSize > items and haveMoreReads:
        isFirstInPair = True                    # always assuming first in pair
        id_line = fh.readline().rstrip()        # the ID line we are going to split
        # for some reason checking the file handler is not enough to find EOF
        # so we are asking for the length of the string read
        if len(id_line) == 0:
            haveMoreReads = False
            break
        seq = fh.readline().rstrip()            # the sequence
        fh.readline()                           # dummy + in FASTQ
        qual = fh.readline().rstrip()           # quality values
        split_line = re.split('[/# ]',id_line)
        readID = split_line[0]
        idx = ''
        if id_line.find("#"):                   # something like @HWUSI-EAS100R:6:73:941:1973#0/1
            isFirstInPair = idx.endswith("1")   # index is something like 0/1 
            idx = "#"+split_line[1]
        elif id_line.find(" "):                 # like @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
            isFirstInPair = idx.startswith("1") # index is something like 1:Y:18:ATCACG
            idx = " "+split_line[1]
        elif id_line.find("/"):                 # like @C09DFACXX111207:1:1101:3739:194647/1
            isFirstInPair = idx.startswith("1") # index is something like 1
            idx = "/"+split_line[1]

        if readID not in reads.keys():          # it is a new read
            reads[readID] = []

        if isFirstInPair:
            reads[readID].insert(0,(idx,seq,qual) ) # we are always inserting into the beginning
        else:
            reads[readID].append( (idx,seq,qual) ) # always append to the end
        items+=1                

    print "read " + str(items) + " reads"
    return reads, haveMoreReads

def doPairing(r1,r2,o1,o2):
    """Topmost method for read pairing
    The data structure is like:
    
    reads{ID} = [(R1_idx_string, sequence, quality), (R2_idx_string, sequence, quality)]
    
    We are writing out any new records where both pairs are present, still, 
    it is not guaranteed that we are not eating up all the memory.
    """
    reads = {}
    chunkSize = 210     # number or reads to get for one go
    while True:
        (reads,haveMoreR1Reads) = getReads(r1,chunkSize,reads)    # read a set of reads from R1 
        (reads,haveMoreR2Reads) = getReads(r2,chunkSize,reads)    # ditto from R1 
        # write the reads and delete the written ones
        for readID in reads.keys():
            writeReadPair(o1, o2, readID, reads[readID])
            #writeReadPair(sys.stdout, sys.stdout, readID, reads[readID])
            del reads[readID]
        if not haveMoreR1Reads and not haveMoreR2Reads:
            break


def writeSingleRead(fh,readID,read):
    """Writes out a single read only
    """
    (idx,seq,qual) = read
    fh.write(readID + idx + "\n")
    fh.write(seq + "\n")
    fh.write("+\n")
    fh.write(qual + "\n")


def writeReadPair(fh_r1,fh_r2,readID,read):
    """Writes a pair of reads
    R1 first, than R2
    """

    if len(read) == 2:
        for (fh,i) in zip([fh_r1,fh_r2],[0,1]):
            writeSingleRead(fh,readID, read[i])

@click.command()
@click.option('--R1',       '-1', type=str, help='The R1 file')
@click.option('--outR1',    '-r', type=str, help='R1 reads written into this file')
@click.option('--R2',       '-2', type=str, help='The R2 file')
@click.option('--outR2',    '-l', type=str, help='R2 reads are written to this')
@click.option('--unpaired', '-u', type=str, help='Unpaired reads are written to this single file', default='unpaired.fastq')
# optionally dummy pairs can be written
def initPairing(r1,r2, outr1,outr2, unpaired):
    """Pairing Illumina reads
    The problem to solve is that for Illumina sequencing the results are usually two files, each containing one
    fragment of a sequence pair. Sometimes one of the pairs is missing, so there order of pairs is broken.
    In most of the cases the nth R1 read in the first file is the pair of the nth read in the R2 file: sometimes it is 
    not the case, and aligners and downstream failing, i.e. for BWA mem you are getting something like: 
    [mem_sam_pe] paired reads have different names: "ST-E00201:119:HTTMMCCXX:5:1112:18182:26378", "ST-E00201:119:HTTMMCCXX:5:1112:15402:26378"
    
    To resolve this, we are reading both files, storing the read IDs in a hash, and writing out reads in pairs
    TODO:
    - write a naive pairing software first
    - find out all the Illumina string formats (space or # delimiter)
    - do it parallel/producer/consumer scheme if slow

    For SRA please, please use fastq-dump --origfmt --split-3 SRR0123456 (see https://en.wikipedia.org/wiki/FASTQ_format#NCBI_Sequence_Read_Archive)
    """

    with open(r1,"r") as r1_fh, \
            open(r2,"r") as r2_fh, \
            open(outr1,"w") as out_r1_fh, \
            open(outr2,"w") as out_r2_fh, \
            open(unpaired,"w") as unpaired_fh:
        doPairing(r1_fh, r2_fh, out_r1_fh, out_r2_fh)

if __name__ == "__main__":
    initPairing()
