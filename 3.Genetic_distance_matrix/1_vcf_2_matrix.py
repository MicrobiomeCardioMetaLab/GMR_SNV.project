import pandas as pd
import os
import time
import numpy as np
import sys
import argparse
from functools import reduce

parser = argparse.ArgumentParser(description="calculate distance matrix from mummer result")


parser.add_argument('-i',"--input", required = True,
		help="the input species")


args = parser.parse_args()
species = args.input


output_file = 'snv_matrix_raw.xls'
allSnvSite = 'allSnvSite.txt'
base2_matrix = 'snv_matrix2base.xls'
base2_filtermatrix = 'snv_matrix2base_filter.xls'
base234_matrix = 'snv_matrix234base.xls'
aln_file = 'allAln.fasta'
base2aln_file = 'base2Aln.fasta'
base2_filteraln_file = 'base2_filterAln.fasta'
base234aln_file = 'base234Aln.fasta'

all_sites = set()
sample_files = os.listdir(os.path.join(mummer_dir,'mummer'))
for file in sample_files:
    vcf_file_path = os.path.join(mummer_dir, 'mummer', file, file + '_qr.snps')
    if not os.path.isfile(vcf_file_path):
        continue
    snp_file = pd.read_table(vcf_file_path)
    snp_file['snv_site'] = (
        snp_file['[FRM]'].astype(str) + '__' +
        snp_file['[P1]'].astype(str) + '__' +
        snp_file['[SUB]'].astype(str)
    )
    all_sites.update(snp_file['snv_site'])


all_sites = sorted(all_sites)
with open(allSnvSite, 'w') as f:
    for site in all_sites:
        f.write(f"{site}\n")



chunk_size = 10000
chunks = [all_sites[i:i+chunk_size] for i in range(0, len(all_sites), chunk_size)]



if os.path.exists(output_file):
    os.remove(output_file)
if os.path.exists(base2_matrix):
    os.remove(base2_matrix)
if os.path.exists(base234_matrix):
    os.remove(base234_matrix)
open(output_file, 'w').close()
open(base2_matrix, 'w').close()
open(base234_matrix, 'w').close()


for ci, sites_chunk in enumerate(chunks, start=1):

    dfs = []
    for file in sample_files:
        vcf_file_path = os.path.join(mummer_dir, species, 'mummer', file, file + '.csv')
        if not os.path.isfile(vcf_file_path):
            continue
            snp_file = pd.read_table(vcf_file_path)
            snp_file['snv_site'] = (
                snp_file['[FRM]'].astype(str) + '__' +
                snp_file['[P1]'].astype(str) + '__' +
                snp_file['[SUB]'].astype(str)
            )
        snp_file = snp_file[snp_file['snv_site'].isin(sites_chunk)]
        df = snp_file[['snv_site', 'ALT']]
        df.rename(columns={'ALT': file}, inplace=True)
        dfs.append(df)

    if not dfs:
        continue

    chunk_df = reduce(lambda left, right: pd.merge(left, right, on='snv_site', how='outer'), dfs)
    chunk_df.fillna('-', inplace=True)
    chunk_df.to_csv(output_file, mode='a', index=False, sep='\t', header=(ci == 1))


    chunk_df.set_index('snv_site', inplace=True)
    unique_counts = chunk_df.apply(lambda x: x.nunique(), axis=1)
    #print(unique_counts)
    snv2 = chunk_df[unique_counts == 3]
    snv234 = chunk_df[unique_counts.isin([3,4, 5])]


    snv2.to_csv(base2_matrix, mode='a', sep='\t', header=(ci == 1))
    snv234.to_csv(base234_matrix, mode='a', sep='\t', header=(ci == 1))



dat = pd.read_csv(output_file, sep='\t')
dat.set_index('snv_site', inplace=True)

with open(aln_file, 'w+') as output:
    for index, col in dat.iteritems():
        if len(set(col)) > 1:
            output.write(f'>{index}\n{"".join(col)}\n')

snv2 = pd.read_csv(base2_matrix, sep='\t')
snv2.set_index('snv_site', inplace=True)
with open(base2aln_file, 'w+') as output:
    for index, col in snv2.iteritems():
        if len(set(col)) > 1:
            output.write(f'>{index}\n{"".join(col)}\n')

snv2_filter1 = snv2.replace('-', '')
snv2_filter2 = snv2_filter1[snv2_filter1.count(axis=1) / snv2_filter1.shape[1] >= 0.05]
counts = snv2_filter2.apply(lambda x: x.value_counts(dropna=True, normalize=True), axis=1)
mask = (counts.max(axis=1) <= 0.95) & (counts.apply(lambda x: x.nlargest(2)[-1], axis=1) >= 0.05)
snv2_filter3 = snv2_filter2[mask]
snv2_filter3.to_csv(base2_filtermatrix, sep='\t')
snv2_filter3.fillna('-', inplace=True)
with open(base2_filteraln_file, 'w+') as output:
    for index, col in snv2_filter3.iteritems():
        if len(set(col)) > 1:
            output.write(f'>{index}\n{"".join(col)}\n')


snv234 = pd.read_csv(base234_matrix, sep='\t')
snv234.set_index('snv_site', inplace=True)
with open(base234aln_file, 'w+') as output:
    for index, col in snv234.iteritems():
        if len(set(col)) > 1:
            output.write(f'>{index}\n{"".join(col)}\n')








