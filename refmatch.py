import sys
import os
import glob
import tempfile
import argparse

import numpy as np
import h5py
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

minimal_match_size = 50

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--referenceFile", help="location of reference .fa file", type=str)
parser.add_argument("-s", "--sequenceFolder", help="folder with sequence files", type=str)
parser.add_argument("-mm", "--minimalMatch", help="minimal size of match that will can be on output", type=int, default = 50)

args = parser.parse_args()

if not os.path.isfile(args.referenceFile):
    print("Referencny subor neexistuje.")
    exit(1)

# recursively find all fast5 files in directory

fast5Files = glob.glob(args.sequenceFolder + '/**/*.fast5', recursive=True)

print("Najdenych " + str(len(fast5Files)) + " suborov")

# create fasta file from sequence strings

fastaSequenceFile = tempfile.NamedTemporaryFile(mode = 'w', suffix = '.fa', delete=False)

for file in fast5Files:
    
    sequenceFile = h5py.File(file, 'r')
    basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
    basecallString = basecallOut[1]
    fastaSequenceFile.write(">" + file + "\n")
    fastaSequenceFile.write(basecallString + "\n")
    sequenceFile.close()

fastaSequenceFile.close()

# load created fasta file

sequenceIndex = mp.Aligner(args.referenceFile)
if not sequenceIndex: raise Exception("ERROR: failed to load/build reference index")

# create out-file and fill it with hits

outFile = open("out.txt", "w")

# header

outFile.write("# <sequenceFileName> <reference_begin> <reference_end> <sequence_begin> <sequence_end>\n\n")

for name, seq, qual in mp.fastx_read(fastaSequenceFile.name): # read a fasta sequence
        for hit in sequenceIndex.map(seq): # traverse alignments
            
            if (hit.r_en - hit.r_st < args.minimalMatch):
                continue
            
            # hit.q_st, hit.q_en  <- hit in reference
            # hit_r_st, hir.r_en  <- hit in sequence
            
            queryIndex = myF(name, hit.q_st, hit.q_en)
            
            if len(queryIndex) == 2:
                outFile.write(name + '\t')
                outFile.write(str(hit.r_st) + '\t' + str(hit.r_en) + '\t')
                outFile.write(str(queryIndex[0]) + "\t" + str(queryIndex[1]) + "\n\n")
