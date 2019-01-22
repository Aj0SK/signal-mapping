import sys
import os
import glob

import h5py
import numpy as np
import mappy as mp


# find signal corresponding to string [hit_beg, hit_end] in sequenceFile with name sequenceFileName 

def myF(sequenceFileName, hit_beg, hit_end):
    # open file
    sequenceFile = h5py.File(sequenceFileName, 'r')
    
    # get raw data and basecallTable
    
    basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]
    readName = str(list(sequenceFile['Raw/Reads'])[0])
    rawData = sequenceFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    
    # output list with raw-start and raw-end
    out = []
    # index to basecall string
    j = -1
    
    # as we are passing through table, did we pass begg or end?
    begg = 0
    endd = 0
    
    for index in range(0, len(basecallEventTable)):
        # we are moving through basecalled string using table
        j = j + basecallEventTable[index][5]
        
        # if our string of length 5 passed beginning, we mark beginning
        if begg == 0 and j+5 > hit_beg:
            begg = 1
            out.append(basecallEventTable[index][1].item())
        
        # if our string of length 5 passed ending, we mark ending
        if begg == 1 and endd == 0 and j>=hit_end:
            endd = 1;
            out.append(basecallEventTable[index][1].item())
            
    sequenceFile.close()
    return out

################################################################################


if len(sys.argv) != 3:
    print("Usage is python refmatch.py <reference.fa> <search_directory>")

# store all temporary data here
if not os.path.exists("temporary"):
    os.makedirs("temporary")

if not os.path.isfile(sys.argv[1]):
    print("Referencny subor neexistuje.")
    exit(1)

# recursively find all fast5 files in directory

fast5Files = glob.glob(sys.argv[2] + '/**/*.fast5', recursive=True)

print("Najdenych " + str(len(fast5Files)) + " suborov")

# create fasta file from sequence strings

fastaSequenceFile = open("temporary/sequence.fa", 'w')

for file in fast5Files:
    
    sequenceFile = h5py.File(file, 'r')
    basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
    basecallString = basecallOut[1]
    
    fastaSequenceFile.write(">" + file + "\n")
    fastaSequenceFile.write(basecallString + "\n")
    sequenceFile.close()

fastaSequenceFile.close()

# load created fasta file

sequenceIndex = mp.Aligner("temporary/sequence.fa")  
if not sequenceIndex: raise Exception("ERROR: failed to load/build reference index")

# create out-file and fill it with hits

outFile = open("out.txt", "w")

for name, seq, qual in mp.fastx_read(sys.argv[1]): # read a fasta sequence
        for hit in sequenceIndex.map(seq): # traverse alignments
            # hit.q_st, hit.q_en  <- hit in reference
            # hit_r_st, hir.r_en  <- hit in sequence
            
            queryIndex = myF(hit.ctg, hit.r_st, hit.r_en)
            
            if len(queryIndex) == 2:
                outFile.write(hit.ctg + " " + str(hit.q_st) + " " + str(hit.q_en) + " " + str(queryIndex[0]) + " " + str(queryIndex[1]) + "\n\n")

os.remove("temporary/sequence.fa")
os.rmdir("temporary")
