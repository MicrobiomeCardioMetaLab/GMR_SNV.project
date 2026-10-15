import pandas as pd
import numpy as np
import argparse
import os
import gzip


parser = argparse.ArgumentParser(description="calculate distance matrix from inStrain snv matrix")


parser.add_argument('-i',"--input", required = True,
                help="the species name")


args = parser.parse_args()
species = args.input



raw_matrix = f'../results/snv_distance/{species}/snv_matrix_raw.xls'
raw_aln = f'../results/snv_distance/{species}/snv_raw.fasta'

base2_matrix = f'../results/snv_distance/{species}/snv_matrix2base_filter.xls'
base2_aln = f'../results/snv_distance/{species}/snv_mtx2base_filter.fasta'

base234_matrix = f'../results/snv_distance/{species}/snv_matrix234base.xls'
base234_aln = f'../results/snv_distance/{species}/snv_mtx234base.fasta'


base2_df = pd.read_table(base2_matrix, index_col=0, low_memory=False)
base2_df.fillna('-', inplace=True)
base2_df_filter = base2_df.loc[:, ~(base2_df == '-').all()]



raw_df = pd.read_table(raw_matrix, index_col=0, low_memory=False)
raw_df.fillna('-', inplace=True)
raw_df_filter = raw_df.loc[:, ~(raw_df == '-').all()]


base234_df = pd.read_table(base234_matrix, index_col=0, low_memory=False)
base234_df.fillna('-', inplace=True)
base234_df_filter = base234_df.loc[:, ~(base234_df == '-').all()]


with open(raw_aln, 'w') as output:
    for sample_name in raw_df_filter.columns:  
        sequence = ''.join(raw_df_filter[sample_name].values)  
        output.write(f'>{sample_name}\n{sequence}\n')



with open(base234_aln, 'w') as output:
    for sample_name in base234_df_filter.columns:
        sequence = ''.join(base234_df_filter[sample_name].values)
        output.write(f'>{sample_name}\n{sequence}\n')



with open(base2_aln, 'w') as output:
    for sample_name in base2_df_filter.columns:
        sequence = ''.join(base2_df_filter[sample_name].values)
        output.write(f'>{sample_name}\n{sequence}\n')



