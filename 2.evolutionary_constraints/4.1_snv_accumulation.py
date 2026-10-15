import os
import numpy as np
import pandas as pd
from Bio import SeqIO
import sys
import argparse
import gzip


def snv_accumulation(fasta_file, repeats=10, out_file=None, verbose=True):


    if fasta_file.endswith(".gz"):
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file, "r")

    with handle:
        seqs = [str(rec.seq) for rec in SeqIO.parse(handle, "fasta")]

    if len(seqs) < 2:
        raise ValueError("Need >= 2 sequences")


    L = len(seqs[0])

    if any(len(s) != L for s in seqs):
        raise ValueError("Sequences are not aligned")

    n = len(seqs)

    if verbose:
        print(f"{fasta_file}: {n} seq × {L} sites")


    matrix = np.frombuffer(
        "".join(seqs).encode("ascii"),
        dtype=np.uint8
    )

    matrix = matrix.reshape(n, L)


    results = np.zeros(
        (n, repeats),
        dtype=np.int32
    )


    for r in range(repeats):

        if verbose:
            print(f"  repeat {r + 1}/{repeats}")


        order = np.random.permutation(n)


        seen_min = matrix[order[0]].copy()
        seen_max = matrix[order[0]].copy()


        polymorphic = np.zeros(L, dtype=bool)

        results[0, r] = 0


        for k in range(1, n):

            seq = matrix[order[k]]

            seen_min = np.minimum(seen_min, seq)
            seen_max = np.maximum(seen_max, seq)

            polymorphic |= (seen_min != seen_max)

            results[k, r] = polymorphic.sum()


    df = pd.DataFrame(
        results,
        columns=[f"rep{i + 1}" for i in range(repeats)]
    )

    df.insert(
        0,
        "sample_size",
        np.arange(1, n + 1)
    )


    if out_file:
        df.to_csv(
            out_file,
            sep="\t",
            index=False
        )

    return df


def main():
    parser = argparse.ArgumentParser(description='SNV Accumulation Curve Calculation')
    parser.add_argument('-i',"--input", required = True, help='Input FASTA alignment file')
    parser.add_argument('--repeats', type=int, default=100, help='Number of random repeats (default: 100)')
    parser.add_argument('-d', '--dir', default='data', help='base directory holding <species>/ (default: data)')

    args = parser.parse_args()

    species = args.input
 
    input_fa = f'../Example/{species}/snv_matrix/snv_mtx2base_filter.fasta'

    if not os.path.exists(input_fa) and os.path.exists(input_fa + '.gz'):
        input_fa = input_fa + '.gz'

    os.makedirs('../results/evolution', exist_ok=True)
    out_file = f'../results/evolution/{species}_snv_accumulation.xls'

     
    df = snv_accumulation(input_fa, repeats=args.repeats, out_file=out_file)
    print("Calculation completed.")

if __name__ == "__main__":
    main()



