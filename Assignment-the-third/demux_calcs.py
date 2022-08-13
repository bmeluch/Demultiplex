#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math

############################## IMPORT RESULTS ##############################

# read in talapas index file to get sample info
index_file: str = "/projects/bgmp/shared/2017_sequencing/indexes.txt"
# index info array columns: sample, group, treatment, index, sequence
index_info = np.empty((24,5), dtype=object)

with open(index_file, 'r') as index_file:
    for lcount, line in enumerate(index_file):
        # skip the first line
        if lcount > 0:
            index_info[lcount-1, :] = line.strip().split()

# read in counts tsv
paircounts_file: str = "/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output/Index_Counts.tsv"
# keys = tuple of index pair, value = count
paircounts: dict = {}

with open(paircounts_file, 'r') as paircounts_file:
    for line in paircounts_file:
        i1, i2, c = line.strip().split()
        paircounts[(i1, i2)] = c
print("Pair counts length (should be 576): ", len(paircounts))

# SORT PAIR COUNTS BY ORDER IN INDEX FILE

# HOW DO I DO THIS


# read in file length info from bash script
linecounts_file: str = "/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output/FASTQ_Line_Counts.txt"
linecounts: dict = {}
with open(linecounts_file, 'r') as linecounts_file:
    for enum, line in enumerate(linecounts_file):
        lc, fname = line.strip().split()
        # all the output file counts have ./ at the beginning, skip input counts then remove first two characters
        if enum > 2:
            fname = fname[2:]
        rc = int(lc)/4
        linecounts[fname] = (lc, rc)
print("Line counts length (number of files, should be 54): ", len(linecounts))

# YOU NOW HAVE
# index_info: 24x5 array, columns: sample, group, treatment, index, sequence
# paircounts: 576 element dict. key = (index1, index2) value = record count
# linecounts: 54 element dictionary. key = filename, value = (line count, record count)

############################## CHECK COUNTS ##############################

# counter variables from Demux_Output.txt 
matchcount: int = 268193924
hoppedcount: int = 416466
unkcount: int = 26964596
qualcount: int = 67671749
rcount: int = 363246735

print("Percentage of reads with matched indices:", round(matchcount/rcount*100, 2), "%")
print("Percentage of reads with hopped indices:", round(hoppedcount/rcount*100, 2), "%")
print("Percentage of reads with unknown indices:", round((unkcount+qualcount)/rcount*100, 2), "%")
print('\t', "Percentage of reads with unrecognized indices:", round(unkcount/rcount*100, 2), "%")
print('\t', "Percentage of reads with quality failed indices:", round(qualcount/rcount*100, 2), "%")


# maths for use in summary markdown
total_input_records: int = linecounts["1294_S1_L008_R1_001.fastq.gz"][1]+linecounts["1294_S1_L008_R4_001.fastq.gz"][1]
print("Total records in input FASTQ files:", int(total_input_records))

total_output_records: int = 0
for fname in linecounts:
    if "Read" in fname:
        total_output_records += linecounts[fname][1]
print("Total records in output FASTQ files:", int(total_output_records))

############################## MATCHED COUNTS ##############################

# matched counts plot
matched_pairs: dict = {}
for key in paircounts:
    # if indices of pair match, save to matched dict
    if key[0] == key[1]:
        matched_pairs[key[0]] = int(paircounts[key])

fig1, ax1 = plt.subplots()
ax1.bar(matched_pairs.keys(), matched_pairs.values())
ax1.set_xlabel("Index")
ax1.set_ylabel("Number of reads")
ax1.set_title("Number of reads per sample")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output/bar_matched.png")
plt.show()
plt.close(fig1)

# matched reads counts
print("Index", "Read count", "Percentage of matched reads")
for i in matched_pairs:
    print(i, matched_pairs[i], round(matched_pairs[i]/matchcount*100, 2), sep=" | ")

############################## INDEX HOPPING ##############################

# heatmap array to populate from paircounts
heatmap = np.zeros((24,24), dtype=int)

# fill heatmap array with counts of each index pair
for r in range(0,24):
    for c in range(0,24):
        # leave out matched pairs so they don't break the color scale
        if r != c:
            heatmap[r,c] = paircounts[(index_info[r,4], index_info[c,4])]

fig2, ax2 = plt.subplots()
ax2.imshow(heatmap)
ax2.set_xlabel("Index 2")
ax2.set_xticks(range(len(index_info[:,4])), labels=index_info[:,4])
ax2.set_ylabel("Index 1")
ax2.set_yticks(range(len(index_info[:,4])), labels=index_info[:,4])
ax2.set_title("Heatmap of reads by index combination (excluding matched reads)")
plt.xticks(rotation=90, ha='right')
plt.tight_layout()
plt.savefig("/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output/hopped_heatmap.png")
plt.show()
plt.close(fig2)
















# Matched total count (matchcount)
#         should equal total records in matched r1 files
#         should equal total records in matched r2 files
#         should equal total of counts for matched pairs in Index_Counts.tsv
#         include breakdown of matched record counts for each of 24 indices

#     Unknown total count
#         should equal unkcount + qualcount
#         should equal total records in r1 unknown file
#         should equal total records in r2 unknown file

#     Index hopped total count
#         should equal hoppedcount
#         should equal number of records in index hopped r1 file
#         should equal number of records in index hopped r2 file
#         should equal total of counts for hopped pairs in Index_Counts.tsv

#     Total number of records in all output files
#         should equal total number of records in input files
#         should equal matchcount + unkcount + qualcount + hoppedcount

############################## PLOTS ##############################

# pie chart of 3 (4) read categories





# bar chart of total reads in each matched index
# y axis on one side in total reads, other side in percent of total matched reads

# heatmap of index pairs