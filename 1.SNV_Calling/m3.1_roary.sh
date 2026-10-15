#!/usr/bin/env bash
#
# Pan-genome construction with Roary (only needed to rebuild the pan-genome
# reference; the submitted Example already contains pan_genome_reference.fa).
#
# Usage:
#   ./m3.1_roary.sh <gff_dir> [outdir]
#
# Optional: THREADS (default 64)
#
set -euo pipefail

GFF_DIR="${1:-gff}"
OUTDIR="${2:-roary}"
THREADS="${THREADS:-64}"

roary -e --mafft -i 90 -cd 90 -s "${GFF_DIR}"/*.gff -f "${OUTDIR}" -p "${THREADS}"