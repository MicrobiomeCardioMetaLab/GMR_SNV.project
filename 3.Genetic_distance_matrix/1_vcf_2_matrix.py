#!/usr/bin/env python3
"""
SNV matrix builder + 2base filter
====================================================================

"""
import os
import csv
import argparse
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed


def build_sample_index(args):
    sample, sample_root, index_dir = args
    in_file = f'{sample_root}/{sample}_qr.snps.vcf'
    out_file = f'{index_dir}/{sample}.tsv'

    if not os.path.isfile(in_file):
        return sample, None

    sites = set()

    try:
        with open(in_file) as f, open(out_file, 'w') as out:
            for _ in range(4):
                f.readline()
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                site = f"{row['#CHROM']}__{row['POS']}__{row['REF']}"
                alt = row['ALT']
                out.write(site + '\t' + alt + '\n')
                sites.add(site)
        return sample, sites
    except Exception as e:
        print(f"[WARN] {sample}: {e}")
        return sample, None


def load_sample_chunk(sample, index_dir, sites_set):
    path = f'{index_dir}/{sample}.tsv'
    d = {}
    try:
        with open(path) as f:
            for line in f:
                site, alt = line.rstrip('\n').split('\t')
                if site in sites_set:
                    d[site] = alt
    except Exception:
        pass
    return sample, d


def process_chunk(ci, sites_chunk, samples, index_dir,
                  out_raw, out_b2, out_b234,
                  first_chunk, threads):
    sites_set = set(sites_chunk)
    results = {}
    with ProcessPoolExecutor(max_workers=threads) as exe:
        futures = [exe.submit(load_sample_chunk, s, index_dir, sites_set)
                   for s in samples]
        for f in as_completed(futures):
            s, d = f.result()
            results[s] = d

    raw_f = open(out_raw, 'a')
    b2_f  = open(out_b2, 'a')
    b234_f= open(out_b234, 'a')
    raw_w = csv.writer(raw_f, delimiter='\t')
    b2_w  = csv.writer(b2_f, delimiter='\t')
    b234_w= csv.writer(b234_f, delimiter='\t')

    if first_chunk:
        header = ['snv_site'] + samples
        raw_w.writerow(header)
        b2_w.writerow(header)
        b234_w.writerow(header)

    for site in sites_chunk:
        alts = []
        for s in samples:
            alt = results[s].get(site, '-')
            alts.append(alt)
        
        row = [site] + alts
        raw_w.writerow(row)

        n_unique = len(set(alts))
        if n_unique == 1 and alts[0] == '-':
            continue
        if n_unique == 3:
            b2_w.writerow(row)
            b234_w.writerow(row)
        elif n_unique in (4, 5):
            b234_w.writerow(row)

    raw_f.close()
    b2_f.close()
    b234_f.close()
    print(f"chunk {ci} done")


def filter_2base_matrix(input_path, output_path, chunksize=5000):
 
    first_chunk = True
    if os.path.exists(output_path):
        os.remove(output_path)

    for chunk in pd.read_table(input_path, index_col=0, chunksize=chunksize, low_memory=False):
        chunk = chunk.replace('-', np.nan)

        mask1 = chunk.count(axis=1) / chunk.shape[1] >= 0.1
        dat = chunk[mask1]
        if dat.shape[0] == 0:
            continue

        def freq_filter(row):
            vc = row.value_counts(dropna=True, normalize=True)
            if len(vc) == 0:
                return False
            if vc.iloc[0] > 0.95:
                return False
            if len(vc) > 1 and vc.iloc[1] >= 0.05:
                return True
            return False

        mask2 = dat.apply(freq_filter, axis=1)
        result = dat[mask2]
        if result.shape[0] == 0:
            continue

        result.to_csv(
            output_path,
            sep='\t',
            mode='w' if first_chunk else 'a',
            header=first_chunk
        )
        first_chunk = False




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--cluster', required=True)
    parser.add_argument('-t', '--threads', type=int, default=16)
    parser.add_argument('--chunk-size', type=int, default=5000, help='chunk size')
    args = parser.parse_args()

    cluster = args.cluster
    threads = args.threads
    chunk_size = args.chunk_size

    mummer_dir = f'../Example/{cluster}/MUMmer'
    out_dir = f'../results/snv_distance/{cluster}'
    os.makedirs(out_dir, exist_ok=True)
    index_dir = os.path.join(out_dir, 'index')
    os.makedirs(index_dir, exist_ok=True)

    sample_root = f'{mummer_dir}/'
    samples = []
    for fname in os.listdir(sample_root):
        if fname.endswith('_qr.snps.vcf'):
            sample = fname.replace('_qr.snps.vcf', '')
            samples.append(sample)
    samples = sorted(samples)



    all_sites = set()
    valid = []
    with ProcessPoolExecutor(max_workers=threads) as exe:
        futures = [exe.submit(build_sample_index, (s, sample_root, index_dir))
                   for s in samples]
        for f in as_completed(futures):
            s, sites = f.result()
            if sites:
                valid.append(s)
                all_sites.update(sites)
    samples = sorted(valid)
    all_sites = sorted(all_sites)


    chunks = [all_sites[i:i+chunk_size]
              for i in range(0, len(all_sites), chunk_size)]

    out_raw = os.path.join(out_dir, 'snv_matrix_raw.xls')
    out_b2  = os.path.join(out_dir, 'snv_matrix2base.xls')
    out_b234= os.path.join(out_dir, 'snv_matrix234base.xls')
    out_b2_filter = os.path.join(out_dir, 'snv_matrix2base_filter.xls')
    
    for f in [out_raw, out_b2, out_b234]:
        if os.path.exists(f):
            os.remove(f)


    for ci, chunk in enumerate(chunks):
        first_chunk = (ci == 0)
        process_chunk(ci, chunk, samples, index_dir,
                      out_raw, out_b2, out_b234,
                      first_chunk, threads)

    filter_2base_matrix(out_b2, out_b2_filter, chunksize=chunk_size)

    print("ALL DONE")


if __name__ == '__main__':
    main()
