import sys
import click

# For general desription see the method assigned as  __main__ below.

def doPairing(r1,r2,o1,o2):
    """Topmost method for read pairing
    The data structure is like:
    
    reads{ID} = [(R1_idx_string, sequence, quality), (R2_idx_string, sequence, quality)]
    
    We are writing out any new records where both pairs are present, still, 
    it is not guaranteeing that we are not eating up all the memory.
    """
    reads = {}
    for bugger in ["egy","ketto","harom"]:
        R1 = (" 1:N:18:1", "ACTACGACTGC", "eeefffggasd")
        R2 = (" 2:N:18:1", "ACTTTGACCCC", "eAAAefEgasd")
        reads["ST-E00201:119:HTTMMCCXX:5:1112:15402:"+str(hash(bugger))] = [R1,R2]
    for readID in reads.keys():
        #writeReadPair(o1, o2, readID, reads[readID])
        writeReadPair(sys.stdout, sys.stdout, readID, reads[readID])


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
    """
    r1_fh = open(r1,"r")
    r2_fh = open(r2,"r")
    out_r1_fh = open(outr1,"w")
    out_r2_fh = open(outr2,"w")
    unpaired_fh = open(unpaired,"w")

    doPairing(r1_fh, r2_fh, out_r1_fh, out_r2_fh)

    r1_fh.close()
    r2_fh.close()
    out_r1_fh.close()
    out_r2_fh.close()
    unpaired_fh.close()

if __name__ == "__main__":
    initPairing()
