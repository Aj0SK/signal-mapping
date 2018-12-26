import sys
import os
import h5py
import numpy as np
import mappy as mp

if len(sys.argv) != 3:
    print("Usage is python refmatch.py <sequence>")

# get reference string from referenceFile given as first parameter
referenceFile = open(sys.argv[1], 'r')
referenceString = referenceFile.read().rstrip('\n')
referenceFile.close()

# read sequence file and obtain details about it

sequenceFile = h5py.File(sys.argv[2], 'r')
basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
basecallString = basecallOut[1]
basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]

# create fasta file from reference string a load it into

fastaReferenceFile = open("test/reference.fa", 'w')
fastaReferenceFile.write(">ref\n")
fastaReferenceFile.write(referenceString)
fastaReferenceFile.close()

referenceIndex = mp.Aligner("test/reference.fa")  
if not referenceIndex: raise Exception("ERROR: failed to load/build reference index")

# create fasta file from sequence string and load it

fastaSequenceFile = open("test/sequence.fa", 'w')
fastaSequenceFile.write(">seq1\n")
fastaSequenceFile.write(basecallString)
fastaSequenceFile.close()

#sequenceIndex = mp.Aligner("test/sequence.fa")  
#if not sequenceIndex: raise Exception("ERROR: failed to load/build reference index")

for name, seq, qual in mp.fastx_read("test/sequence.fa"): # read a fasta/q sequence
        for hit in referenceIndex.map(seq): # traverse alignments
            print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            #s = referenceFile.seq("MT_human", hit.r_st, hit.r_en)    
            #print(s)


# print(a[i], end = '')

os.remove("test/reference.fa")
os.remove("test/sequence.fa")
