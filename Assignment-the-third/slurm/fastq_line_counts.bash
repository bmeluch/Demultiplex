#!/bin/bash
cd /projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-third/output
for filename in ./*.fastq.gz; do
    echo "$(zcat $filename | wc -l) $filename" >> FASTQ_Line_Counts.txt
done