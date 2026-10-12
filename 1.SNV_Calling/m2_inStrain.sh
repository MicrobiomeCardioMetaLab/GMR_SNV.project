#!/bin/bash

threads=32
sample=SRR1944846


## step_1. genome profile
bowtie2 -p ${threads} -x ${db}/GMR_DB -1 data/${sample}_1.fastq.gz -2 data//${sample}_2.fastq.gz | samtools sort -@ ${threads} -O bam -o - > data/${sample}.bam

## step_2 get SNV file
#mkdir -p ${dir}/sg_results/${sample}

inStrain profile data//${sample}.bam ${db}/GMR_DB.fasta -o data/ -p ${threads} -g ${db}/GMR_DB_prodigal_genes.fna -s ${db}/GMR_DB.stb --database_mode --use_full_fasta_header --force_compress --skip_plot_generation


