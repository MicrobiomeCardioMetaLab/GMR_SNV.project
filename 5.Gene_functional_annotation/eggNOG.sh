#!/usr/bin/env bash
#
# Functional annotation with eggNOG-mapper (v2.1.12, eggNOG DB 5.0.2).
#
# Usage:
#   ./eggNOG.sh <proteins.faa[.gz]> [sample_name]
#
# Optional environment variables:
#   EGGNOG_DATA_DIR  directory holding the eggNOG database (default: /data/eggNOG)
#   THREADS          number of threads (default: 16)
#   OUTDIR           output directory  (default: data/eggNOG)
#
set -euo pipefail

PROT="$1"
SAMPLE="${2:-$(basename "${PROT}")}"
SAMPLE="${SAMPLE%.gz}"; SAMPLE="${SAMPLE%.faa}"

EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-/data/eggNOG}"
THREADS="${THREADS:-16}"
OUTDIR="${OUTDIR:-data/eggNOG}"

mkdir -p "${OUTDIR}/${SAMPLE}"

emapper.py -i "${PROT}" \
    -o "${SAMPLE}" \
    --data_dir "${EGGNOG_DATA_DIR}" \
    --temp_dir "${OUTDIR}/${SAMPLE}" \
    --output_dir "${OUTDIR}/${SAMPLE}" \
    --cpu "${THREADS}"