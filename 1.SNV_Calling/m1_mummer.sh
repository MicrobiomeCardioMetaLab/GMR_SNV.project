#!/bin/bash

ref_genome=data/GT_GMRSGB03426_ref.fasta
input_genome=data/GT_GMRSGB03426_000002.fasta

####step1. nucmer generates nucleotide alignments between two FASTA files

nucmer -p ${input_genome} ${ref_genome}.fasta ${input_genome}.fa -t 2

####step2 delta-filte filter
delta-filter -q -r ${input_genome}.delta > ${input_genome}_qr.filter


####step3 show snp
show-snps -CIlr -x 1 -T ${input_genome}_qr.filter > ${input_genome}_qr.snps

