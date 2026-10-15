#!/usr/bin/env bash
#
# pairwise ML distances with IQ-TREE (model K2P).
#
# Usage:
#   ./3_iqtree.sh [alignment.fasta] [prefix]
#

species=$1

cd ../results/snv_distance/${species}
iqtree3 -s snv_mtx2base_filter.fasta -m K2P --prefix snv_mtx2base_filter -T 1