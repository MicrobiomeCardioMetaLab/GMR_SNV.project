import os
import pandas as pd
import argparse
import numpy as np
from scipy.stats import entropy

parser = argparse.ArgumentParser(description="calculate distance matrix from snippy result")


parser.add_argument('-i',"--input", required = True,
                help="the full names of input species")

args = parser.parse_args()

species = args.input


snv_mtx_file = f"{species}/snv_matrix2base_filter.xls.gz"
geneLen_file = f"{species}/{species}.geneLen.txt.gz"
scafLen_file = f"{species}_scafLen.xls"
snvGene_file = f"{species}.tsv.gz"

count_df = pd.read_table("GMR_rep_meta.xls")


count_N = count_df.loc[count_df['genome'] == species, 'N_genomePerSp']
total_mag = count_N.iloc[0]

gene_len = pd.read_table(geneLen_file, header=None, names=['locus_tag', 'length'])
scaf_len = pd.read_table(scafLen_file)
scaf_len = scaf_len.drop(['col1', 'col2'], axis=1)
scaf_len['scaffold'] = scaf_len['scaffold'].str.split('|', expand=True)[[1]]
scaf_dict = scaf_len.set_index('scaffold')['length'].to_dict()
gene_dict = gene_len.set_index('locus_tag')['length'].to_dict()


df = pd.read_table(snvGene_file)
df['scaf_len'] = df['scaffold'].map(scaf_dict)
df['locus_tag_len'] = df['locus_tag'].map(gene_dict)

df['eachScaf_N_site'] = df.groupby('scaffold')['scaffold'].transform('size')
df['eachGene_N_site'] = df.groupby('locus_tag')['locus_tag'].transform('size').where(df['locus_tag'] != '.')
cds_df = df[df['feature_type'] == 'CDS']
inter_df = df[df['feature_type'] != 'CDS']

df['gene_total_freq'] = df['eachGene_N_site'] / df['locus_tag_len']
df['gene_total_freq_perMAG'] = df['gene_total_freq'] / df['total_alleles']

df['scaf_total_freq'] = df['eachScaf_N_site'] / df['scaf_len']
df['scaf_total_freq_perMAG'] = df['scaf_total_freq'] / df['total_alleles']



grouped = cds_df.groupby('locus_tag')


results = []
for gene, subdf in grouped:
    length = subdf['locus_tag_len'].iloc[0]
    n_counts = subdf['total_alleles'].values
    total_mutations = n_counts.sum()



    p = n_counts / total_mutations 
    shannon_H = entropy(p, base=2) if len(p) > 1 else 0  




    max_entropy = np.log2(total_mag)
    normalized_H = shannon_H / max_entropy if max_entropy > 0 else 0
    burden_entropy = (total_mutations / length) * normalized_H


    results.append({
        'locus_tag': gene,
        'gene_N_mutations': total_mutations,
        'entropy_weighted_burden': burden_entropy,
    })

result_df = pd.DataFrame(results)

dat = pd.merge(df, result_df, on='locus_tag', how='outer')
dat = dat.sort_values(by='entropy_weighted_burden', ascending=False)
print(dat)


dat.to_csv(f"{species}_SE.xls", index=False, sep='\t')




