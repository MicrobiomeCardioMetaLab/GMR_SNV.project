#!/bin/bash


ref=data/GT_GMRSGB03426_pan_genome.fa.gz
query=data/GT_GMRSGB03426_000002.fasta

snippy --reference ${ref} --ctgs ${query} --outdir ./  --cpus 2 --mincov 5 --force

