import scipy.stats as stats
import pandas as pd
import numpy as np
import argparse
import os
from statsmodels.stats.multitest import multipletests



pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


parser = argparse.ArgumentParser(
    description="KEGG pathway enrichment analysis of SNV-associated genes"
)

parser.add_argument(
    '-i',"--input", 
    required = False,
    help="the names of input species"
)


parser.add_argument(
    '-e',"--pheno",
    required = False,
    help="the names of pheno type"
)


parser.add_argument(
    '-p',"--parameter",
    default='greater',
    required = False,
    help="the parameter of fisher_exact test (two-sided /greater / less)"
)



args = parser.parse_args()



genome = args.input
pheno_id = args.pheno
alter_para = args.parameter



def py_pAdjust(pvalues, method="fdr_bh"):
    pvalues = pd.Series(pvalues, copy=True)

    valid = pvalues.notna()

    adjusted = pd.Series(np.nan, index=pvalues.index, dtype=float)

    if valid.sum() > 0:
        adjusted.loc[valid] = multipletests(
            pvalues.loc[valid],
            method=method
        )[1]
    return adjusted

pathway_df = pd.read_table("../data/ko_name_annotation.xls", header=None, names=["pathway_id", "pathway_name"])
pathway_df["pathway_name"] = (pathway_df["pathway_name"].str.replace(r" - .*?$", "", regex=True))
pathway_name_dict = dict(zip(pathway_df["pathway_id"], pathway_df["pathway_name"]))
#print(pathway_df)


ko_pathway = pd.read_table("../data/ko2pathway.txt", header=None, names=["ko", "pathway_id"])
ko_pathway = ko_pathway[ko_pathway['pathway_id'].str.contains('map')]
ko_pathway_dict = dict(zip(ko_pathway["ko"], ko_pathway["pathway_id"]))
#print(ko_pathway)


gene_ko_df = pd.read_table(f"../Example/{genome}/reference/{genome}.cdsKO.txt")
gene_ko_df.columns = ['gene_id', 'ko_id']
print(f"Total genes: {len(gene_ko_df):,}")

gene_ko_df['pathway_id'] = gene_ko_df['ko_id'].map(ko_pathway_dict)
gene_ko_df = gene_ko_df.dropna(subset=["pathway_id"])
#print(gene_ko_df)
print(f"Total genes with KO: {gene_ko_df['gene_id'].nunique():,}")




snv_df = pd.read_table(f'../Example/{genome}/snv_matrix/{genome}.xls')
snv_df = snv_df.drop_duplicates(['snv_site'])
#snv_df = snv_df[snv_df['fdr'] < 0.05]
snv_df = snv_df[snv_df['feature_type'] == 'CDS']
snv_df = snv_df[['locus_tag']]
#print(snv_df)
snv_genes = set(snv_df['locus_tag'])
print(f"Total SNV genes: {len(snv_genes)}")


snv_annotated = pd.merge(gene_ko_df, snv_df, left_on='gene_id', right_on='locus_tag', how='inner')
snv_annotated = snv_annotated.drop_duplicates(['gene_id', 'ko_id', 'pathway_id'])
#print(snv_annotated)
print(f"SNV genes with KO annotation: {len(set(snv_annotated['locus_tag'])):,}")


background_genes = set(gene_ko_df["gene_id"])
foreground_genes = set(snv_annotated['locus_tag'])


print("\nEnrichment universe:")
print(f"Background genes = {len(background_genes)}")
print(f"Foreground SNV genes = {len(foreground_genes)}")



gene_ko_count = pd.DataFrame(gene_ko_df['pathway_id'].value_counts())
snv_ko_count = pd.DataFrame(snv_annotated['pathway_id'].value_counts())

result = pd.merge(gene_ko_count, snv_ko_count, on='pathway_id',  how='outer', suffixes=("_sp_gene", "_snv_gene"))
result = result.fillna(0)
#print(result)


fisher_result = []
for _, row in result.iterrows():
    sp_path = int(row['count_sp_gene'])
    snv_path = int(row['count_snv_gene'])

    snv_nopath = len(foreground_genes) - snv_path
    nosnv_path = sp_path - snv_path
    nosnv_nopath = len(background_genes) - snv_path - snv_nopath - nosnv_path

    contingency_table = [
        [snv_path, snv_nopath],  
        [nosnv_path, nosnv_nopath]  
    ]
    ko_path = row.name
    n_bg = len(background_genes)
    n_snv = len(foreground_genes)
    odds_ratio, p_value = stats.fisher_exact(contingency_table, alternative=alter_para)
    fisher_result.append([ko_path, n_bg, n_snv, snv_path, snv_nopath, nosnv_path, nosnv_nopath, odds_ratio, p_value])

fisher_df = pd.DataFrame(fisher_result, columns=['ko_path', 'N_backGene', 'N_snvGene', 'snv_path', 'snv_nopath', 'nosnv_path', 'nosnv_nopath', 'odds_ratio', 'pvalue'])
fisher_df["FDR"] = py_pAdjust(fisher_df["pvalue"])


fisher_df = fisher_df.sort_values('pvalue')
fisher_df['pathway_name'] = fisher_df['ko_path'].map(pathway_name_dict)
fisher_df['genome'] = genome
print(fisher_df)


plt_df = fisher_df[fisher_df['pvalue'] < 0.05]
os.makedirs('../results/evolution', exist_ok=True)
fisher_df.to_csv(f'../results/evolution/{genome}_{alter_para}_KEGG.xls.gz', sep='\t', index=False, compression='gzip')




