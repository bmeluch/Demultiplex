#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt
import argparse
import gzip

# Demultiplex Assignment the First
# Part 1 Question 2
# Generate a per base distribution of quality scores for read1, read2, index1, and index2. Average the quality scores at each position for all reads and generate a per nucleotide mean distribution.

# accept filename as argument
parser = argparse.ArgumentParser(description="Get FASTQ filename")
parser.add_argument("-f", help="FASTQ file to run this script on", required=True)
parser.add_argument("-o", help="Output file path and name", required=True)
args = parser.parse_args()
fname: str = args.f
outname: str = args.o

# test file
# fname = "/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-first/test.fastq"

# init line counter
lcount: int = 0
# init sequence length
seqlen: int = 0

# get sequence length
with gzip.open(fname, 'rt') as fq:
    for line in fq:
        lcount += 1
        if lcount == 2:
            seqlen = len(line.strip('\n'))
            break

# create empty lists equal to sequence length
qscores: list = [0]*seqlen
mean: list = [0.0]*seqlen

# reset line counter
lcount = 0
# open the file
with gzip.open(fname,'rt') as fq:
    for line in fq:
        lcount+=1
        # for all quality lines
        if lcount%4==0:
            # iterate over characters in the quality line
            for position in range(len(line.strip())):
                # convert the character to phred score
                # add each new score to the sublist in the corresponding position
                qscores[position] += bioinfo.convert_phred(line[position])
        if lcount%4000000 == 0:
            print("Working on read: ", lcount/4)

# calculate means
for position in range(len(qscores)):
    mean[position] = qscores[position]/(lcount/4)

# plot histogram
plt.title("Average quality score across position in read")
plt.xlabel("Base position")
plt.ylabel("Average quality score")
plt.bar(range(len(mean)), mean)
plt.savefig(outname)
plt.show()