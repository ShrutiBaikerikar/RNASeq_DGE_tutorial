#!/bin/bash

SAMPLES="SRR12852623 SRR12852624 SRR12852625 SRR12852626 SRR12852627 SRR12852628"

for SAMPLE in $SAMPLES; do
    stringtie -e -p 4 -B -G ./assembly/merged_transcript/stringtie_merged.gtf -o ./quantification/${SAMPLE}/${SAMPLE}_requant.gtf ./assembly/sorted_bam_reads/${SAMPLE}_sorted.bam
	echo "$SAMPLE reestimated"

done