#!/usr/bin/env bash
#
# SNV calling with MUMmer4 (nucmer -> delta-filter -q -r -> show-snps).
#
# Every MAG of one species is aligned against that species' reference genome,
# one run per MAG.  The resulting <sample>_qr.snps tables are the input of
# 3.Genetic_distance_matrix/1_vcf_2_matrix.py.
#
# Usage:
#   ./m1_mummer.sh <reference.fasta> <sample1.fa> [sample2.fa ...]
#

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <reference.fasta> <sample.fa> [sample2.fa ...]" >&2
    exit 1
fi

REF="$1"
shift

OUTDIR="${OUTDIR:-MUMmer}"
mkdir -p "${OUTDIR}"

for QUERY in "$@"; do
    SAMPLE="$(basename "${QUERY}")"
    SAMPLE="${SAMPLE%.fa}"
    SAMPLE="${SAMPLE%.fasta}"
    SAMPLE="${SAMPLE%.fna}"

    echo "[INFO] ${SAMPLE}: nucmer vs $(basename "${REF}")"
    nucmer -p "${OUTDIR}/${SAMPLE}"  "${REF}" "${QUERY}"

    echo "[INFO] ${SAMPLE}: delta-filter -q -r"
    delta-filter -q -r "${OUTDIR}/${SAMPLE}.delta" > "${OUTDIR}/${SAMPLE}_qr.filter"

    echo "[INFO] ${SAMPLE}: show-snps"
    show-snps -CIlr -x 1 -T "${OUTDIR}/${SAMPLE}_qr.filter" > "${OUTDIR}/${SAMPLE}_qr.snps"

    gzip -f "${OUTDIR}/${SAMPLE}.delta" "${OUTDIR}/${SAMPLE}_qr.filter" "${OUTDIR}/${SAMPLE}_qr.snps"
done

echo "[DONE] MUMmer output in ${OUTDIR}/"