#!/bin/bash

sample=GT_GMRSGB03426_000002


mkdir -p data/eggNOG/${sample}

emapper.py  -i data/${sample}.faa.gz  -o ${sample}  --temp_dir data/eggNOG/${sample}  --output_dir data/eggNOG/${sample}  --cpu 16




