#!/bin/bash

SAMPLES="SRR12852623 SRR12852624 SRR12852625 SRR12852626 SRR12852627 SRR12852628"

for SAMPLE in $SAMPLES; do
    stringtie -p 4 -G ./assembly/annotation/Homo_sapiens_chr.gtf -o ./assembly/output/${SAMPLE}.gtf -l ${SAMPLE} ./assembly/sorted_bam_reads/${SAMPLE}_sorted.bam
	echo "$SAMPLE assembled"

done