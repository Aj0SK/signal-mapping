import sys
import os
import glob

import h5py
import numpy as np
import mappy as mp

def myF(table, hit_beg, hit_end):
    x = []
    j = -1
    
    for index in range(0, len(table)):
        j = j + table[index][5]
        if (j>= hit_beg and j<=hit_end) or (j+5> hit_beg and j+5<=hit_end):
            x.append(table[index][1])
            print(table[index][4])
        
    return x

if len(sys.argv) != 3:
    print("Usage is python refmatch.py <reference.fa>")

# read sequence files and obtain details about it

sequenceFiles = []

fast5Files = glob.glob(sys.argv[2] + '/**/*.fast5', recursive=True)
for file in fast5Files:
    sequenceFiles.append(h5py.File(file, 'r'))

# create fasta file from sequence string and load it

fastaSequenceFile = open("test/sequence.fa", 'w')

for i in range(0, len(sequenceFiles)):
    fastaSequenceFile.write(">" + fast5Files[i] + "\n")
    basecallOut = str(sequenceFiles[i]['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
    basecallString = basecallOut[1]
    fastaSequenceFile.write(basecallString + "\n")

fastaSequenceFile.close()

sequenceIndex = mp.Aligner("test/sequence.fa")  
if not sequenceIndex: raise Exception("ERROR: failed to load/build reference index")

'''
basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]

# referenceIndex = mp.Aligner("test/reference.fa")
# if not referenceIndex: raise Exception("ERROR: failed to load/build reference index")
'''

for name, seq, qual in mp.fastx_read(sys.argv[1]): # read a fasta/q sequence
        for hit in sequenceIndex.map(seq): # traverse alignments
            print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            #pole = myF(basecallEventTable, hit.r_st, hit.r_en)
            #print("Pole je" + str(pole))
            #s = referenceFile.seq("MT_human", hit.r_st, hit.r_en)    
            #print(s)

os.remove("test/sequence.fa")
