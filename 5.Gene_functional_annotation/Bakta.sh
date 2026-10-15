#!/usr/bin/env bash
#
# Genome annotation with Bakta.  Needs the Bakta database (~/bakta/db, v6.0).
#
# Usage:
#   ./Bakta.sh <genome.fasta> [sample_name]
#
# Optional environment variables:
#   BAKTA_DB   path to the Bakta database      (default: bakta/db)
#   THREADS    number of threads               (default: 16)
#   OUTDIR     output directory                (default: data)
#
set -euo pipefail

GENOME="$1"
SAMPLE="${2:-$(basename "${GENOME}")}"
SAMPLE="${SAMPLE%.fna}"; SAMPLE="${SAMPLE%.fa}"; SAMPLE="${SAMPLE%.fasta}"
SAMPLE="${SAMPLE%.gz}";  SAMPLE="${SAMPLE%.fna}"; SAMPLE="${SAMPLE%.fa}"; SAMPLE="${SAMPLE%.fasta}"

BAKTA_DB="${BAKTA_DB:-bakta/db}"
THREADS="${THREADS:-16}"
OUTDIR="${OUTDIR:-data}"

mkdir -p "${OUTDIR}/${SAMPLE}"

bakta --db "${BAKTA_DB}" \
    --output "${OUTDIR}/${SAMPLE}" \
    --threads "${THREADS}" \
    --prefix "${SAMPLE}" \
    --tmp-dir "${OUTDIR}/${SAMPLE}" \
    --locus-tag "${SAMPLE}" \
    --force \
    --keep-contig-headers \
    --meta \
    "${GENOME}"

# compress the Bakta tables that pNpS.py / snv_frequency.py expect
find "${OUTDIR}/${SAMPLE}" -type f \( -name '*.tsv' -o -name '*.ffn' -o -name '*.faa' -o -name '*.fna' \) \
    -exec gzip -f {} +