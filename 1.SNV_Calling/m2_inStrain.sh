#!/usr/bin/env bash
#
# Read-based SNV calling with inStrain (bowtie2 -> samtools sort -> inStrain profile).
#
# Usage:
#   ./m2_inStrain.sh <sample_name> <read_1.fastq.gz> <read_2.fastq.gz>
#
#   GMR_DB        bowtie2 index prefix, e.g. /data/GMR_DB
#   GMR_DB_FASTA  matching FASTA,           e.g. /data/GMR_DB.fasta     (inStrain reference)
#   GMR_DB_GENES  prodigal gene FASTA,      e.g. /data/GMR_DB_prodigal_genes.fna
#   GMR_DB_STB    scaffold-to-bin file,     e.g. /data/GMR_DB.stb
#
# Optional:
#   THREADS (default 32), OUTDIR (default data)
#
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <sample> <reads_1.fastq.gz> <reads_2.fastq.gz>" >&2
    exit 1
fi

SAMPLE="$1"
R1="$2"
R2="$3"

THREADS="${THREADS:-32}"
OUTDIR="${OUTDIR:-data}"

: "${GMR_DB:?set GMR_DB to the bowtie2 index prefix}"
: "${GMR_DB_FASTA:?set GMR_DB_FASTA to the reference FASTA}"
: "${GMR_DB_GENES:?set GMR_DB_GENES to the prodigal gene FASTA}"
: "${GMR_DB_STB:?set GMR_DB_STB to the scaffold-to-bin file}"

mkdir -p "${OUTDIR}"

## step 1. read mapping
bowtie2 -p "${THREADS}" -x "${GMR_DB}" -1 "${R1}" -2 "${R2}" \
    | samtools sort -@ "${THREADS}" -O bam -o "${OUTDIR}/${SAMPLE}.bam" -

## step 2. inStrain profile
inStrain profile "${OUTDIR}/${SAMPLE}.bam" "${GMR_DB_FASTA}" \
    -o "${OUTDIR}/" \
    -p "${THREADS}" \
    -g "${GMR_DB_GENES}" \
    -s "${GMR_DB_STB}" \
    --database_mode \
    --use_full_fasta_header \
    --force_compress \
    --skip_plot_generation