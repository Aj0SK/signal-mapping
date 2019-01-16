import sys
import os
import glob

import h5py
import numpy as np
import mappy as mp

def myF(sequenceFileName, hit_beg, hit_end):
    sequenceFile = h5py.File(sequenceFileName, 'r')
    
    basecallEventTable = sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Events'][()]
    readName = str(list(sequenceFile['Raw/Reads'])[0])
    rawData = sequenceFile['Raw/Reads/' + readName + "/" + "Signal"][()]
    
    out = []
    j = -1
    
    zaciatok = 0
    koniec = 0
    
    for index in range(0, len(basecallEventTable)):
        j = j + basecallEventTable[index][5]
        
        if zaciatok == 0 and j+5 > hit_beg:
            zaciatok = 1
            out.append(basecallEventTable[index][1].item())
        
        if zaciatok == 1 and koniec == 0 and j>=hit_end:
            koniec = 1;
            out.append(basecallEventTable[index][1].item())
            
    sequenceFile.close()
    return out

################################################################################


if len(sys.argv) != 3:
    print("Usage is python refmatch.py <reference.fa> <search_directory>")

if not os.path.exists("temporary"):
    os.makedirs("temporary")

if not os.path.isfile(sys.argv[1]):
    print("Referencny subor neexistuje.")
    exit(1)

sequenceFilesNames = []
sequenceFilesIndexer = {}

# create fasta file from sequence string and load it

fastaSequenceFile = open("temporary/sequence.fa", 'w')

fast5Files = glob.glob(sys.argv[2] + '/**/*.fast5', recursive=True)

print("Najdenych " + str(len(fast5Files)) + " suborov")

for file in fast5Files:
    
    sequenceFilesNames.append(file)
    sequenceFilesIndexer[file] = len(sequenceFilesIndexer)
    
    sequenceFile = h5py.File(file, 'r')
    basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
    basecallString = basecallOut[1]
    
    fastaSequenceFile.write(">" + file + "\n")
    fastaSequenceFile.write(basecallString + "\n")
    sequenceFile.close()

fastaSequenceFile.close()

sequenceIndex = mp.Aligner("temporary/sequence.fa")  
if not sequenceIndex: raise Exception("ERROR: failed to load/build reference index")

outFile = open("out.txt", "w")

for name, seq, qual in mp.fastx_read(sys.argv[1]): # read a fasta/q sequence
        for hit in sequenceIndex.map(seq): # traverse alignments
            #print("{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en))
            index = sequenceFilesIndexer[hit.ctg]
            pole = myF(sequenceFilesNames[index], hit.r_st, hit.r_en)
            outFile.write(hit.ctg)
            if len(pole) != 0:
                for elem in pole:
                    outFile.write(" " + str(elem))
                outFile.write("\n")

os.remove("temporary/sequence.fa")
os.rmdir("temporary")
