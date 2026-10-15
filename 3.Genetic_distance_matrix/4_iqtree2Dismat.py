import os
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="calculate k2p distance matrix from iqtree result")


parser.add_argument('-i',"--input", required = True,
                help="the name of input species")

parser.add_argument('--prefix', default='snv_mtx2base_filter',
                help="IQ-TREE run prefix (default: snv_mtx2base_filter)")

args = parser.parse_args()
species = args.input
prefix = args.prefix



iq_distance = pd.read_csv(f'../results/snv_distance/{species}/{prefix}.mldist', skiprows=1, header=None, sep=r"\s+", dtype=str)

sample_id = iq_distance.iloc[:, 0].tolist()
df = iq_distance.iloc[:, 1:]
df.index = sample_id
df.columns = sample_id

df_numeric = df.apply(pd.to_numeric, errors='coerce')
df_scaled = df_numeric * 100
mask_nan = df_scaled.isna()

max_val = np.nanmax(df_scaled.values)

df_scaled = df_scaled.fillna(max_val + 1)

dis = df_scaled

output_file = f'../results/snv_distance/{species}/{prefix}.mldist_pivotDis.xls'
dis.to_csv(output_file, sep='\t', index=True, header=True)

print(f"Distance matrix saved to {output_file}")





















