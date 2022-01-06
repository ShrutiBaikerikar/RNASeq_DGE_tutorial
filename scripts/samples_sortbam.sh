#!/bin/bash

SAMPLES="SRR12852623 SRR12852624 SRR12852625 SRR12852626 SRR12852627 SRR12852628"

for SAMPLE in $SAMPLES; do
    samtools sort -@ 4 -o ./assembly/sorted_bam_reads/${SAMPLE}_sorted.bam ./alignment/output/${SAMPLE}_aligned.sam
	echo "$SAMPLE sorted"

done