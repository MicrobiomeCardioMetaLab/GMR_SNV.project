# Single-nucleotide variations in prevalent gut microbes are lineage-specific and associated with functional divergence relevant to host physiology



This repository contains the SNV calling, processed intermediate data, randomization analyses and association scripts associated with the study:

Single-nucleotide variations of prevalent gut microbes are lineage-specific functional signatures and associations with host traits


---

## Overview

The analysis is organized around four major components:

1. **SNV Calling**
   - Identification of SNV in the MAG-based GMR catalogue using MUMmer.
   - Identification of SNV in the reads-based CGMR cohort using inStrain.


2. **Species evolutionary constraints**
   - Calculate pN / pS ratio, snv_frequency, Shannon entropy across species.
   - Functional interpretation using KEGG Orthology profiles.
   - KO-SNV Presence/absence profiling across gut microbial species.
   - Jaccard-distance analyses and phylogenetic/taxonomic structure.


3. **SNVsite-Sample matrix**
   - Generate basetype-sample matrix.
   - Calculate within-species genetic distance.
   - Generate phylogenetic tree using IQTREE.


4. **Association with phenotypes**
   - Pheno adnois snvDis.
   - snvDis association with phenoDis using Mantel.
   - snvDis association with phenoDis using Spearman.
   - Fitted phylogenetic generalized linear mixed models (PGLMMs) using anpan.


---

## Repository structure

```text
25-CF-IBD-MGX/
├── README.md
├── Example/                    # An example of all analysis
├── data/                       # Species SNV-profile summary data
├── 1.SNV_Calling/              # SNV Calling pipeline
├── 2.evolutionary constraints/ # Calculate evolutionary constraints index
├── 3.Genetic_distance_matrix/  # SNVsite-Sample matrix
├── 4.Association_with_pheno/   # Association with phenotypes

```


---

## Public cohorts

| Cohort | BioProject |
|---|---|
| Zachary_2022 | PRJNA834801 |


---


## Software environment

Analyses were performed under:

- Ubuntu 24.04.4
- Python 3.12.3
- R 4.5.3


### Command-line software

| Software | Version / notes |
|---|---|
| Anpan |  |
| DIAMOND | 2.0.11 |
| eggNOG-mapper | 2.1.12 |
| eggNOG database | 5.0.2 |
| MUMmer | 4.0.0rc1 |
| inStrain | 1.9.0 |
| Bakta | 1.12.0 |
| Bakta database | 6.0 |
| Prodigal | 2.6.3 |
| Bowtie2 | 2.2.5 |
| Samtools | 1.15.1 |
| Kraken2 | 2.0.8 |
| Bracken | 3.0.1 |
| IQtree | 3.0.1 |
| iTOL | 7.3 |
| seqkit | 0.10.0 |
| SNP-site | 2.5.1 |
| VCFtools | 0.1.17 |

### Major R packages

| Package | Version |
|---|---|
| vegan | 2.7 |
| ape | 5.8 |
| clusterProfiler | 4.17 |
| ggplot2 | 3.5.2 |


### Major Python packages

| Package | Version |
|---|---|
| NumPy | 2.1.3 |
| pandas | 2.2.3 |

---

## External resources

This project relies on external reference resources that are not fully redistributed here:

- GMR genome catalogue and metadata;
- GMR Kraken2 database;
- eggNOG v5.0;
- Bakta v6.0;
- public sequencing datasets listed above.

---

