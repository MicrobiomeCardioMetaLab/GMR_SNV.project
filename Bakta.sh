#!/bin/bash
sample=GT_GMRSGB03426_000002
threads=16
db_dir=bakta/db


mkdir -p data/${sample}

bakta --db ${db_dir} \
 --output data/${sample} \
 --threads ${threads} \
 --prefix ${sample} \
 --tmp-dir data/${sample} \
 --locus-tag ${sample} \
 --force \
 --keep-contig-headers \
 --meta \
 data/${sample}.fasta

pigz * -p ${threads}


