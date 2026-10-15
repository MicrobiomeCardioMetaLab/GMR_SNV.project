#!/usr/bin/env bash
#
# Taxonomic profiling with Kraken2 + Bracken.
#
# Usage:
#   ./kraken.sh <sample_name> <reads_1.fastq.gz> <reads_2.fastq.gz>
#
# Optional environment variables:
#   KRAKEN_DB  Kraken2 database directory (default: kk2_GMR)
#   READ_LEN   read length for Bracken    (default: 150)
#   THREADS    number of threads          (default: 16)
#   OUTDIR     output directory           (default: data)
#
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <sample> <reads_1.fastq.gz> <reads_2.fastq.gz>" >&2
    exit 1
fi

SAMPLE="$1"; R1="$2"; R2="$3"

KRAKEN_DB="${KRAKEN_DB:-kk2_GMR}"
READ_LEN="${READ_LEN:-150}"
THREADS="${THREADS:-16}"
OUTDIR="${OUTDIR:-data}"

mkdir -p "${OUTDIR}/kraken" "${OUTDIR}/bracken"

# build the database once (uncomment):
# kraken2-build --build --threads 12 --db "${KRAKEN_DB}"
# bracken-build -d "${KRAKEN_DB}" -t 8 -k 35 -l 150

kraken2 --use-names --threads "${THREADS}" --db "${KRAKEN_DB}" --confidence 0 \
    --output "${OUTDIR}/kraken/${SAMPLE}.output.txt" \
    --report "${OUTDIR}/kraken/${SAMPLE}.report.txt" \
    --report-zero-counts --gzip-compressed --paired "${R1}" "${R2}"

bracken -d "${KRAKEN_DB}" -i "${OUTDIR}/kraken/${SAMPLE}.report.txt" \
    -o "${OUTDIR}/bracken/${SAMPLE}.txt" -r "${READ_LEN}" -l R7