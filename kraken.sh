#!/bin/bash

sample=SRR1944846

kk_database=kk2_GMR


#kraken2-build --build --threads 12 --db .
kraken2 --use-names  --threads 16  --db ${kk_database}  --confidence 0  --output data/kraken/${sample}.output.txt  --report data/kraken/${sample}.report.txt  --report-zero-counts  --gzip-compressed  --paired data/${sample}_1.fastq.gz data/${sample}_2.fastq.gz


#bracken-build -d ${kk_database} -t 8 -k 35 -l 150
bracken -d ${kk_database} -i data/kraken/${sample}.report.txt -o data//bracken/${sample}.txt  -r 150 -l R7



