#!/bin/bash

SAMPLES="SRR12852623 SRR12852624 SRR12852625 SRR12852626 SRR12852627 SRR12852628"

for SAMPLE in $SAMPLES; do
    echo "$SAMPLE alignment started"
	hisat2 -p 4 -f -x ./alignment/index/genome -q -U ./alignment/input_reads/${SAMPLE}_trimmed.fastq -S ./alignment/output/${SAMPLE}_aligned.sam
	echo "$SAMPLE alignment successful"

done