#! /usr/bin/env python 

from Bio import AlignIO
import numpy
import sys
import getopt
import statistics as sts

arguments = sys.argv[1:]

## code mode default off
code = 0


try:  
    opts, args = getopt.getopt(arguments, "-s:-m:-o:-t:-h-c",["sequence","mask","output","threshold","help","code"])
except getopt.GetoptError:  
    # print help information and exit: 
    print("Your options have some error")
    print('''Modified by Shaolin.Xu (2020-07-01) from script from Schwentner, Martin (DOI:https://doi.org/10.5061/)
        -s (--sequence) for input sequence, require string
        -m (--mask) for mask file, require string
        -t (--threshold) for threshold, require string
        -o (--output) for output file, require string
        -c (--code) use code mode
        -h (--help) for help info''')
    exit()

for opt, arg in opts:
    if opt in ['-h','--help']:
        print('''Modified by Shaolin.Xu (2020-07-01) from script from Schwentner, Martin (DOI:https://doi.org/10.5061/)
        -s (--sequence) for input sequence, require string
        -m (--mask) for mask file, require string
        -t (--threshold) for threshold, require string
        -o (--output) for output file, require string
        -c (--code) use code mode
        -h (--help) for help info''')
        exit()
    elif opt in ['-s','--sequence']:
        seq = arg
    elif opt in ['-m','--mask']:
        mask = arg
    elif opt in ['-o','--output']:
        out = arg
    elif opt in ['-t','--threshold']:
        threshold = float(arg)
    elif opt in ['-c','--code']:
        code = 1

#threshold = 5 #float(raw_input('What HMM value do you which to exclude? '))

alignment = AlignIO.read(seq, 'fasta')
align_array = numpy.array([list(rec) for rec in alignment], numpy.character, order="F")


## get the number of nucleotides in the alignment from mask file
LN = 0
with open(mask, 'r') as f:
    for line in f:
        LN += 1


### The number of codon
codonN = int(LN/3)

zorro_handle = open(mask, 'r')
count = 0
cut = align_array[:,:1]

if code==1:
        for j in range(codonN):
                codon1 = float(zorro_handle.readline().rstrip("\n"))
                codon2 = float(zorro_handle.readline().rstrip("\n"))
                codon3 = float(zorro_handle.readline().rstrip("\n"))
                codon123 = (codon1+codon2+codon3)/3
                if codon123 >= threshold:
                        cut = numpy.column_stack((cut, align_array[:,(j*3):(j*3+3)]))
        cut = cut[:,1:]
        zorro_handle.close()
else:
        for j in zorro_handle:
                if float(j.rstrip('\n')) >= threshold:
                        cut = numpy.column_stack((cut, align_array[:,count]))
                count += 1
        cut = cut[:,1:]
        zorro_handle.close()        

outfile = open(out, 'w')
for q in range(len(alignment)):
        seq = ''
        for p in cut[q,:]:
                seq += bytes.decode(p)
        outfile.write('>' + alignment[q].id + '\n' + seq + '\n')
outfile.close()
