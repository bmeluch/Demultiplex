#!/usr/bin/env python

import bioinfo
import matplotlib.pyplot as plt
import argparse
import gzip
import math
import numpy as np
import itertools

# assumption: the records are in the same order in each file

############################## FUNCTIONS ##############################

def write_demux_fastq(read1_record: list, read2_record: list, out_files: tuple):
    """Takes record information and tuple of 2 file objects
    and writes record to the correct demultiplexed FASTQ file.""" 
    for i in range(0,4):
        out_files[0].write(read1_record[i]+'\n')
        out_files[1].write(read2_record[i]+'\n')
    pass

############################## GET INPUT FILES ##############################

# # argparse filenames and settings#
# parser = argparse.ArgumentParser(description="Get FASTQ filename")
# parser.add_argument("-r1", help="R1 FASTQ file (read1)", required=True)
# parser.add_argument("-r2", help="R2 FASTQ file (index1)", required=True)
# parser.add_argument("-r3", help="R3 FASTQ file (index2)", required=True)
# parser.add_argument("-r4", help="R4 FASTQ file (read2)", required=True)
# parser.add_argument("-o", help="Output folder path and name, INCLUDE /", required=True)
# args = parser.parse_args()
# read1_fname: str = args.r1
# index1_fname: str = args.r2
# index2_fname: str = args.r3
# read2_fname: str = args.r4
# out_fold: str = args.o

# local test files
read1_fname = "/home/bmeluch/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/unit-test_read1_R1.fastq.gz"
index1_fname = "/home/bmeluch/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/unit-test_index1_R2.fastq.gz"
index2_fname = "/home/bmeluch/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/unit-test_index2_R3.fastq.gz"
read2_fname = "/home/bmeluch/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/unit-test_read2_R4.fastq.gz"
index_fname = "/home/bmeluch/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/unit-test_known-index.txt"
out_fold = "/home/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output/"


############################## IMPORT INDICES ##############################
# Read file with list of known indices, save 24 index sequences into a set
known_indices = set()
with open(index_fname, 'r') as indices:
    for line in indices:
        i = line.strip().split()
        # if the 5th column is a DNA sequence, add it to set of known indices
        if bioinfo.validate_base_seq(i[4]):
            known_indices.add(i[4])
print("Known indices (should be 24): ", len(known_indices))
print(known_indices)

# use itertools product to produce tuples of all possible index pairs
# pair_counts dictionary: keys = tuples of index pairs; values: count of observed pairs (starts 0)
pair_counts = dict.fromkeys(itertools.product(known_indices, repeat=2), 0)
print("Potential index combinations (should be 576): ", len(pair_counts))


############################## OPEN OUTPUT FILES ##############################
# output filename format: Matched_INDEX-INDEX_Read1.fastq.gz, Hopped_Read1.fastq.gz, Unknown_Read1.fastq.gz

# out_names dictionary: keys = condition; values = tuple (Read1, Read2) file objects
# Conditions: <index> (if matched), "Unknown", "Hopped")
out_names: dict = {}
out_names["Unknown"] = (gzip.open(out_fold+"Unknown_Read1.fastq.gz", 'wt'), gzip.open(out_fold+"Unknown_Read2.fastq.gz", 'wt'))
out_names["Hopped"] = (gzip.open(out_fold+"Hopped_Read1.fastq.gz", 'wt'), gzip.open(out_fold+"Hopped_Read2.fastq.gz", 'wt'))

# open files labeled with matched index pairs
for i in known_indices:
    out_names[i] = (gzip.open(out_fold+"Matched_"+i+"-"+i+"_Read1.fastq.gz", 'wt'), \
        gzip.open(out_fold+"Matched_"+i+"-"+i+"_Read2.fastq.gz", 'wt'))
print("Output file names (should be 26): ", len(out_names))


############################## COUNTERS AND HOLDERS ##############################
# read category counts
matchcount: int = 0
unkcount: int = 0
qualcount: int = 0
hoppedcount: int = 0
# total record counter
rcount: int = 0
# quality score cutoff
qual_cutoff: int = 30

# holder variables for records: list of four strings
# 0: header, 1: seq, 2: +, 3: qual
read1_record: list = [""]*4
read2_record: list = [""]*4
index1_record: list = [""]*4
index2_record: list = [""]*4


############################## MAIN LOGIC ##############################
# Open four input files
with gzip.open(read1_fname, 'rt') as r1, gzip.open(index1_fname, 'rt') as i1, gzip.open(index2_fname, 'rt') as i2, gzip.open(read2_fname, 'rt') as r2:
    # Read lines in files simultaneously
    while True:

        # store the next four lines of each file as a list
        read1_record = list(itertools.islice(r1, 4))
        read2_record = list(itertools.islice(r2, 4))
        index1_record = list(itertools.islice(i1, 4))
        index2_record = list(itertools.islice(i2, 4))

        # If you read in nothing, stop the while loop
        if not read1_record:
            break

        # strip newlines from each record line
        read1_record = [rec.strip() for rec in read1_record]
        read2_record = [rec.strip() for rec in read2_record]
        index1_record = [rec.strip() for rec in index1_record]
        index2_record = [rec.strip() for rec in index2_record]

        # reverse complement index2
        index2_record[1] = bioinfo.reverse_complement(index2_record[1])
        # reverse index2 quality line so the qual characters stay in line with the bases
        index2_record[3] = index2_record[3][::-1]

        # optional - error correcting - if single n in index, correct with same character from other index
        # bring over quality score as well
        erri1 = 0
        erri2 = 0
        for b in range(len(index2_record[1])):
            if index2_record[1][b] == "N":
                erri2 += 1
                if erri2 > 1:
                    break
                index2_record[1] = index2_record[1][:b]+index1_record[1][b]+index2_record[1][(b+1):]
                index2_record[3] = index2_record[3][:b]+index1_record[3][b]+index2_record[3][(b+1):]
            elif index1_record[1][b] == "N":
                erri1 += 1
                if erri1 > 1:
                    break
                index1_record[1] = index1_record[1][:b]+index2_record[1][b]+index1_record[1][(b+1):]
                index1_record[3] = index1_record[3][:b]+index2_record[3][b]+index1_record[3][(b+1):]

        # append indices to headers
        read1_record[0] = read1_record[0]+"_"+index1_record[1]+"-"+index2_record[1]
        read2_record[0] = read2_record[0]+"_"+index1_record[1]+"-"+index2_record[1]

        # 1. Are both indices known? (compare to set of indices, inherently checks for N)
        if index1_record[1] not in known_indices and index2_record[1] not in known_indices:
            # if both are known, continue to next question
            # if either index is unknown, the read is unknown
            write_demux_fastq(read1_record, read2_record, out_names["Unknown"])
            unkcount += 1

        # 2. Are they high enough quality?
        elif not bioinfo.qual_pass(index1_record[3], qual_cutoff) or not bioinfo.qual_pass(index2_record[3], qual_cutoff):
            # if both pass quality, continue to next question
            # if either index has a low quality base, the read quality fails
            write_demux_fastq(read1_record, read2_record, out_names["Unknown"])
            qualcount += 1

        # 3. Do they match?
        elif index1_record[1] == index2_record[1]:
            # if indices match, the read is matched
            write_demux_fastq(read1_record, read2_record, out_names[index1_record[1]])
            matchcount += 1
            pair_counts[(index1_record[1], index2_record[1])] += 1
        else:
            # if indices do not match, the read is hopped
            write_demux_fastq(read1_record, read2_record, out_names["Hopped"])
            hoppedcount += 1
            pair_counts[(index1_record[1], index2_record[1])] += 1

        # count that a record has been filed
        rcount += 1

# close all output files
for f in out_names:
    out_names[f][0].close()
    out_names[f][1].close()


############################## PRETTY OUTPUT ##############################

print("matchcount: ", matchcount)
print("unkcount: ", unkcount)
print("qualcount: ", qualcount)
print("hoppedcount: ", hoppedcount)
print("readcount: ", rcount)

with open(out_fold+"Demux_Output.txt", 'w') as dot, open(out_fold+"Index_Counts.tsv", 'w') as ict:
    dot.write("Total reads with matched indices:\t"+str(matchcount)+'\n')
    dot.write("Total reads with hopped indices:\t"+str(hoppedcount)+'\n')
    dot.write("Total reads with unknown indices:\t"+str(unkcount+qualcount)+'\n')
    dot.write("\tTotal reads with unrecognized indices:\t"+str(unkcount)+'\n')
    dot.write("\tTotal reads with quality failed indices:\t"+str(qualcount)+'\n')
    dot.write("Total reads written to files:\t"+str(rcount)+'\n')

    for i in pair_counts:
        ict.write(i[0]+'\t'+i[1]+'\t'+str(pair_counts[i])+'\n')