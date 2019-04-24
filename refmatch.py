import sys
import os
import glob
import tempfile
import argparse
import numpy as np
import h5py
import mappy as mp
import re

out_header = "# <sequenceFileName> <reference_begin> <reference_end> <sequence_begin> <sequence_end>"

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--referenceFile", help="location of reference .fa file", type=str)
parser.add_argument("-s", "--sequenceFolder", help="folder with sequence files", type=str)
parser.add_argument("-mm", "--minimalMatch", help="minimal size of match that will can be on output", type=int, default = 50)
parser.add_argument("-o", "--outputFile", help="specifies output file", type=str, default = "out.txt")
parser.add_argument('-raw', action='store_true', help="output raw signal")
parser.add_argument('-fake', action='store_true', help="create fake read from hit")


class Table_Iterator:

    def __init__(self, basecallEventTable):
        self.table = basecallEventTable
        self.tableindex = 0
        self.localindex = 0
    def __iter__(self):
        return self

    def __next__(self):

        while self.localindex == 5:
            if self.tableindex + 1 != len(self.table):
                self.tableindex += 1
                self.localindex = 5-int(self.table[self.tableindex][5])
            else:
                raise StopIteration

        self.localindex += 1
        return self.table[self.tableindex][4][self.localindex-1]

def myF(args, sequenceFileName, hit_beg, hit_end, csSequence = None):
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
    begg, endd = False, False
    
    for index in range(0, len(basecallEventTable)):
        # we are moving through basecalled string using table
        j = j + basecallEventTable[index][5]
        
        # if our string of length 5 passed beginning, we mark beginning
        if begg == False and j+5 > hit_beg:
            begg = True
            out.append(basecallEventTable[index][1].item())
        
        # if our string of length 5 passed ending, we mark ending
        if begg == True and j>=hit_end:
            out.append(basecallEventTable[index][1].item())
            endd = True
            break

    if begg == True and endd == False:
        out.append(basecallEventTable[index][1].item())

    for i in range(out[0], out[1]):
        out.append(rawData[i].item())

    '''
    sequenceIterator = Table_Iterator(basecallEventTable)

    csParsed = re.findall(r':[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+', csSequence)

    #print(str(csParsed))

    index, current, reference_index = 0, 0, 0
    for parsed in csParsed:
        if parsed[0] == ':':

        elif csParsed[index][0] == '-':
            reference_index += len(parsed)-1
    '''
    sequenceFile.close()
    return out

################################################################################

args = parser.parse_args()

assert os.path.isfile(args.referenceFile), "Reference file not exists."

# recursively find all fast5 files in directory

fast5Files = glob.glob(args.sequenceFolder + '/**/*.fast5', recursive=True)

print("Found %d .fast5 files\n" % (len(fast5Files)))

# create fasta file from sequence strings

fastaSequenceFile = tempfile.NamedTemporaryFile(mode = 'w', suffix = '.fa', delete = False)

for file in fast5Files:
    
    sequenceFile = h5py.File(file, 'r')
    basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')
    basecallString = basecallOut[1]
    fastaSequenceFile.write(">" + file + "\n")
    fastaSequenceFile.write(basecallString + "\n")
    sequenceFile.close()

fastaSequenceFile.close()

# index reference File

sequenceIndex = mp.Aligner(args.referenceFile)
assert sequenceIndex, "failed to load/build reference index"

# create out-file and fill it with hits

outFile = open(args.outputFile, "w")

# header

outFile.write(out_header + "\n\n")

for name, seq, qual in mp.fastx_read(fastaSequenceFile.name): # read a fasta sequence
        for hit in sequenceIndex.map(seq, cs = True): # traverse alignments
            
            if (hit.r_en - hit.r_st < args.minimalMatch):
                continue
            
            # hit.r_st, hit.r_en  <- hit in reference
            # hit.q_st, hir.q_en  <- hit in sequence


            queryIndex = myF(args, name, hit.q_st, hit.q_en, hit.cs)
            
            if len(queryIndex) >= 2:
                outFile.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (name, hit. ctg, hit.r_st, hit.r_en, queryIndex[0], queryIndex[1], hit.strand) )
            else:
                print("Match not found in sequence file, could be caused by corruption of data.")
                exit(1)

            if args.raw:
                for signal in queryIndex[2:]:
                    outFile.write(str(signal) + ' ')
                outFile.write("\n")
            
            
            outFile.write("\n")


outFile.close()
