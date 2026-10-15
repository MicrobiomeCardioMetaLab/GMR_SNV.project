#!/usr/bin/env bash
#
# Variant calling against a pan-genome reference with Snippy.
#
# Usage:
#   ./m3.2_snippy.sh <reference.fasta[.gz]> <contigs.fasta> [outdir]
#
# Optional: THREADS (default 2), MINCOV (default 5)
#
set -euo pipefail

REF="$1"
QUERY="$2"
THREADS="${THREADS:-2}"
MINCOV="${MINCOV:-5}"

mkdir -p ../results/snippy/${QUERY}

../soft/snippy/bin/snippy --reference "${REF}" --ctgs "${QUERY}" --outdir /results/snippy/${QUERY} \
    --cpus "${THREADS}" --mincov "${MINCOV}" --force