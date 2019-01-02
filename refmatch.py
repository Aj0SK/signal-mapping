import sys
import os
import glob

import h5py
import numpy as np
import mappy as mp

def myF(sequenceFile, hit_beg, hit_end):
    basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]
    readName = str(list(sequenceFile['Raw/Reads'])[0])
    rawData = sequenceFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    
    out = []
    j = -1
    
    for index in range(0, len(basecallEventTable)):
        j = j + basecallEventTable[index][5]
        if (j>= hit_beg and j<=hit_end) or (j+5> hit_beg and j+5<=hit_end):
            out.append(rawData[basecallEventTable[index][1].item()].item())
        
    return out




if len(sys.argv) != 3:
    print("Usage is python refmatch.py <reference.fa>")

# read sequence files and obtain details about it

sequenceFiles = []
sequenceFilesIndexer = {}

fast5Files = glob.glob(sys.argv[2] + '/**/*.fast5', recursive=True)
for file in fast5Files:
    sequenceFiles.append(h5py.File(file, 'r'))
    sequenceFilesIndexer[file] = len(sequenceFilesIndexer)

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

outFile = open("out.txt", "w")

for name, seq, qual in mp.fastx_read(sys.argv[1]): # read a fasta/q sequence
        for hit in sequenceIndex.map(seq): # traverse alignments
            #print("{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en))
            index = sequenceFilesIndexer[hit.ctg]
            pole = myF(sequenceFiles[index], hit.r_st, hit.r_en)
            outFile.write(hit.ctg)
            if len(pole) != 0:
                for elem in pole:
                    outFile.write(" " + str(elem))
                outFile.write("\n")
                    

os.remove("test/sequence.fa")
