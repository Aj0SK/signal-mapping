import sys
import h5py
import numpy as np
import mappy as mp

if len(sys.argv) != 2:
    print("Usage is python refmatch.py <sequence>")

sequenceFile = h5py.File(sys.argv[1], 'r')

basecallOut = str(sequenceFile['/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'][()]).split('\\n')

basecallString = basecallOut[1]

print(basecallString)


referenceFile = mp.Aligner("test/MT-human.fa")  # load or build index
if not referenceFile: raise Exception("ERROR: failed to load/build index")

for name, seq, qual in mp.fastx_read("test/MT-orang.fa"): # read a fasta/q sequence
        for hit in referenceFile.map(seq): # traverse alignments
            print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            #s = referenceFile.seq("MT_human", hit.r_st, hit.r_en)    
            #print(s)

#print(a[i], end = '')
