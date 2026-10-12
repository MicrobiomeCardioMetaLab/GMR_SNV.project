#!/usr/bin/env python3

import os
import gzip
import argparse
from bisect import bisect_right
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq



def open_text(filename, mode="rt"):

    if filename.endswith(".gz"):
        return gzip.open(filename, mode)

    return open(filename, mode)



def read_bakta_annotation(bakta_file):

    print("[INFO] Reading Bakta annotation...")

    temp = defaultdict(list)

    with open_text(bakta_file, "rt") as f:



        for _ in range(5):

            line = next(f, None)

            if line is None:

                raise ValueError(
                    "Bakta file contains fewer than 5 lines."
                )


        header_line = next(f, None)

        if header_line is None:

            raise ValueError(
                "Cannot find Bakta header."
            )

        header = (
            header_line
            .rstrip("\n")
            .split("\t")
        )

        header[0] = header[0].lstrip("\ufeff")

        col = {
            x.strip(): i
            for i, x in enumerate(header)
        }



        required = [
            "#Sequence Id",
            "Type",
            "Start",
            "Stop",
            "Strand",
            "Locus Tag",
        ]

        missing = [
            x for x in required
            if x not in col
        ]

        if missing:

            raise ValueError(
                "Missing Bakta columns:\n"
                + "\n".join(missing)
                + "\n\nAvailable columns:\n"
                + "\n".join(header)
            )

        gene_col = col.get("Gene")
        product_col = col.get("Product")
        dbxrefs_col = col.get("DbXrefs")



        for line in f:

            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            try:

                feature_type = fields[
                    col["Type"]
                ].strip()


                if feature_type.lower() != "cds":
                    continue

                scaffold = fields[
                    col["#Sequence Id"]
                ].strip()

                start = int(
                    fields[col["Start"]]
                )

                stop = int(
                    fields[col["Stop"]]
                )

                strand = fields[
                    col["Strand"]
                ].strip()

                locus_tag = fields[
                    col["Locus Tag"]
                ].strip()

                gene = (
                    fields[gene_col].strip()
                    if gene_col is not None
                    else ""
                )

                product = (
                    fields[product_col].strip()
                    if product_col is not None
                    else ""
                )

                dbxrefs = (
                    fields[dbxrefs_col].strip()
                    if dbxrefs_col is not None
                    else ""
                )

            except (
                IndexError,
                ValueError
            ):

                continue


            if start > stop:

                start, stop = stop, start

            temp[scaffold].append({

                "start": start,

                "stop": stop,

                "strand": strand,

                "locus_tag": locus_tag,

                "gene": gene,

                "product": product,

                "dbxrefs": dbxrefs,
            })



    annotations = {}

    for scaffold, records in temp.items():

        records.sort(
            key=lambda x: (
                x["start"],
                x["stop"]
            )
        )

        annotations[scaffold] = {

            "records": records,

            "starts": [
                x["start"]
                for x in records
            ]
        }

    n_cds = sum(
        len(x["records"])
        for x in annotations.values()
    )

    print(
        f"[INFO] CDS records: {n_cds:,}"
    )

    print(
        f"[INFO] Scaffolds containing CDS: "
        f"{len(annotations):,}"
    )

    return annotations




def read_ffn(ffn_file):

    print("[INFO] Reading Bakta .ffn...")

    cds_sequences = {}

    with open_text(ffn_file, "rt") as f:

        for record in SeqIO.parse(
            f,
            "fasta"
        ):


            locus_tag = record.id.strip()

            cds_sequences[locus_tag] = str(record.seq).upper()

    print(
        f"[INFO] FFN CDS sequences: "
        f"{len(cds_sequences):,}"
    )

    return cds_sequences



def check_annotation_ffn(
    annotations,
    cds_sequences
):

    bakta_ids = set()

    for scaffold_data in annotations.values():

        for record in scaffold_data["records"]:

            bakta_ids.add(
                record["locus_tag"]
            )

    ffn_ids = set(
        cds_sequences.keys()
    )

    missing = (
        bakta_ids - ffn_ids
    )

    extra = (
        ffn_ids - bakta_ids
    )

    print(
        "\n[INFO] Bakta <-> FFN consistency check"
    )

    print(
        f"[INFO] Bakta CDS : "
        f"{len(bakta_ids):,}"
    )

    print(
        f"[INFO] FFN CDS   : "
        f"{len(ffn_ids):,}"
    )

    print(
        f"[INFO] Missing FFN: "
        f"{len(missing):,}"
    )

    print(
        f"[INFO] Extra FFN: "
        f"{len(extra):,}"
    )

    if missing:

        print(
            "[WARNING] First 10 CDS missing from FFN:"
        )

        for x in list(missing)[:10]:

            print(
                "   ",
                x
            )




def parse_snv_site(snv_site):
 
    try:

        species, rest = snv_site.split(
            "|",
            1
        )

        # --------------------------------------------
        # 只需要 scaffold + position
        # --------------------------------------------

        scaffold, position, _ = rest.split(
            "__",
            2
        )

        return (
            species,
            scaffold,
            int(position)
        )

    except Exception:

        return None




def find_overlapping_cds(
    scaffold,
    position,
    annotations
):

    scaffold_data = annotations.get(
        scaffold
    )

    if scaffold_data is None:

        return []

    starts = scaffold_data["starts"]

    records = scaffold_data["records"]



    right = bisect_right(
        starts,
        position
    )

    if right == 0:

        return []

    hits = []



    for i in range(
        right - 1,
        -1,
        -1
    ):

        record = records[i]

        if record["stop"] >= position:

            hits.append(record)

    return hits



def get_major_alleles(row):

    values = row.dropna()

    if len(values) == 0:

        return None


    values = (
        values.astype(str)
        .str.upper()
    )


    values = values[
        ~values.isin([
            "-",
            "",
            "N",
            "NA",
            "NAN",
        ])
    ]


    values = values[
        values.isin([
            "A",
            "C",
            "G",
            "T",
        ])
    ]

    if len(values) == 0:

        return None

    counts = values.value_counts()


    if len(counts) < 2:

        return None

    base1 = counts.index[0]

    base2 = counts.index[1]

    count1 = int(
        counts.iloc[0]
    )

    count2 = int(
        counts.iloc[1]
    )

    total = int(
        counts.sum()
    )

    return (
        base1,
        base2,
        count1,
        count2,
        total
    )




def genomic_to_cds_position(
    position,
    start,
    stop,
    strand
):
    """
    genomic coordinate  1-based。

    CDS position  1-based。

    + strand:

        position - start + 1

    - strand:

        stop - position + 1
    """

    if strand == "+":

        return (
            position
            - start
            + 1
        )

    elif strand == "-":

        return (
            stop
            - position
            + 1
        )

    else:

        return None




def classify_two_alleles(
    cds_sequence,
    cds_position,
    base1,
    base2,
    strand
):


    if cds_position < 1:

        return None

    if cds_position > len(
        cds_sequence
    ):

        return None



    index = cds_position - 1



    codon_start = (
        index // 3
    ) * 3

    codon_position = (
        index % 3
    )


    codon_ref = cds_sequence[
        codon_start:
        codon_start + 3
    ]

    if len(codon_ref) != 3:

        return None



    if strand == "-":

        allele1_cds = str(
            Seq(base1)
            .reverse_complement()
        )

        allele2_cds = str(
            Seq(base2)
            .reverse_complement()
        )

    elif strand == "+":

        allele1_cds = base1

        allele2_cds = base2

    else:

        return None



    if allele1_cds not in "ACGT":

        return None

    if allele2_cds not in "ACGT":

        return None


    codon1 = (
        codon_ref[:codon_position]
        + allele1_cds
        + codon_ref[codon_position + 1:]
    )



    codon2 = (
        codon_ref[:codon_position]
        + allele2_cds
        + codon_ref[codon_position + 1:]
    )



    try:

        aa1 = str(
            Seq(codon1).translate()
        )

        aa2 = str(
            Seq(codon2).translate()
        )

    except Exception:

        return None

    if aa1 == aa2:

        mutation_type = "synonymous"

    else:

        mutation_type = "nonsynonymous"

    return {



        "reference_cds_base":
            codon_ref[codon_position],


        "base1":
            base1,

        "base2":
            base2,


        "base1_cds":
            allele1_cds,

        "base2_cds":
            allele2_cds,


        "codon_ref":
            codon_ref,

        "codon1":
            codon1,

        "codon2":
            codon2,



        "aa1":
            aa1,

        "aa2":
            aa2,



        "mutation_type":
            mutation_type,
    }



def process_chunk(
    chunk,
    annotations,
    cds_sequences,
    output_handle
):

    n_processed = 0

    n_cds_hit = 0

    n_classified = 0

    n_synonymous = 0

    n_nonsynonymous = 0

    for snv_site, row in chunk.iterrows():

        n_processed += 1



        parsed = parse_snv_site(
            str(snv_site)
        )

        if parsed is None:

            continue

        (
            species,
            scaffold,
            position
        ) = parsed



        alleles = get_major_alleles(
            row
        )

        if alleles is None:

            continue

        (
            base1,
            base2,
            count1,
            count2,
            total
        ) = alleles



        hits = find_overlapping_cds(
            scaffold,
            position,
            annotations
        )



        if not hits:

            output_handle.write(
                "\t".join([

                    str(snv_site),

                    species,

                    scaffold,

                    str(position),


                    ".",

                    base1,
                    base2,
                    str(count1),
                    str(count2),
                    str(total),

                    "intergenic",

                    ".",
                    ".",
                    ".",

                    ".",
                    ".",
                    ".",
                    ".",

                    ".",

                    ".",
                    ".",
                    ".",

                    ".",
                    ".",
                    ".",
                    ".",

                    ".",
                    ".",

                ])
                + "\n"
            )

            continue



        n_cds_hit += 1

        for cds in hits:

            start = cds["start"]

            stop = cds["stop"]

            strand = cds["strand"]

            locus_tag = cds["locus_tag"]

            gene = cds["gene"]

            product = cds["product"]

            dbxrefs = cds["dbxrefs"]



            cds_sequence = cds_sequences.get(
                locus_tag
            )

            if cds_sequence is None:

                output_handle.write(
                    "\t".join([

                        str(snv_site),
                        species,
                        scaffold,
                        str(position),
                        ".",

                        base1,
                        base2,
                        str(count1),
                        str(count2),
                        str(total),

                        "CDS",

                        str(start),
                        str(stop),
                        strand,

                        locus_tag,
                        gene,
                        product,
                        dbxrefs,

                        ".",

                        ".",
                        ".",
                        ".",

                        ".",
                        ".",

                        ".",
                        ".",

                        ".",

                        "missing_ffn",
                    ])
                    + "\n"
                )

                continue



            cds_position = genomic_to_cds_position(
                position,
                start,
                stop,
                strand
            )

            if cds_position is None:

                continue



            result = classify_two_alleles(
                cds_sequence=cds_sequence,

                cds_position=cds_position,

                base1=base1,

                base2=base2,

                strand=strand,
            )

            if result is None:

                output_handle.write(
                    "\t".join([

                        str(snv_site),
                        species,
                        scaffold,
                        str(position),
                        ".",

                        base1,
                        base2,
                        str(count1),
                        str(count2),
                        str(total),

                        "CDS",

                        str(start),
                        str(stop),
                        strand,

                        locus_tag,
                        gene,
                        product,
                        dbxrefs,

                        str(cds_position),

                        ".",
                        ".",
                        ".",

                        ".",
                        ".",

                        ".",
                        ".",

                        ".",

                        "invalid_position",
                    ])
                    + "\n"
                )

                continue


            n_classified += 1

            if result[
                "mutation_type"
            ] == "synonymous":

                n_synonymous += 1

            elif result[
                "mutation_type"
            ] == "nonsynonymous":

                n_nonsynonymous += 1



            output_handle.write(
                "\t".join([



                    str(snv_site),

                    species,

                    scaffold,

                    str(position),

                    ".",


                    base1,

                    base2,

                    str(count1),

                    str(count2),

                    str(total),


                    "CDS",

                    str(start),

                    str(stop),

                    strand,

                    locus_tag,

                    gene,

                    product,

                    dbxrefs,


                    str(cds_position),

                    result[
                        "reference_cds_base"
                    ],

                    result[
                        "base1_cds"
                    ],

                    result[
                        "base2_cds"
                    ],


                    result[
                        "codon_ref"
                    ],

                    result[
                        "codon1"
                    ],

                    result[
                        "codon2"
                    ],



                    result[
                        "aa1"
                    ],

                    result[
                        "aa2"
                    ],


                    result[
                        "mutation_type"
                    ],

                ])
                + "\n"
            )

    return (
        n_processed,
        n_cds_hit,
        n_classified,
        n_synonymous,
        n_nonsynonymous
    )



def run_pipeline(
    bakta_file,
    ffn_file,
    snv_matrix,
    output_file,
    chunk_size=100000
):

    print("=" * 70)

    print(
        "SNV -> CDS -> codon -> amino acid"
    )

    print("=" * 70)


    annotations = read_bakta_annotation(
        bakta_file
    )


    cds_sequences = read_ffn(
        ffn_file
    )



    check_annotation_ffn(
        annotations,
        cds_sequences
    )

    print(
        "\n[INFO] Processing SNV matrix..."
    )

    print(
        f"[INFO] Chunk size: "
        f"{chunk_size:,}"
    )


    total_processed = 0

    total_cds_hit = 0

    total_classified = 0

    total_synonymous = 0

    total_nonsynonymous = 0



    with gzip.open(
        output_file,
        "wt"
    ) as out:

        out.write(
            "\t".join([

                "snv_site",

                "species",

                "scaffold",

                "position",

                "snv_ref_ignored",

                "base1",

                "base2",

                "base1_count",

                "base2_count",

                "total_alleles",

                "feature_type",

                "start",

                "stop",

                "strand",

                "locus_tag",

                "gene",

                "product",

                "dbxrefs",

                "cds_position",

                "reference_cds_base",

                "base1_cds",

                "base2_cds",

                "codon_ref",

                "codon1",

                "codon2",

                "aa1",

                "aa2",

                "mutation_type",

            ])
            + "\n"
        )


        reader = pd.read_table(
            snv_matrix,

            index_col=0,

            chunksize=chunk_size,

            low_memory=False,
        )

        for chunk_id, chunk in enumerate(
            reader,
            start=1
        ):

            (
                n_processed,
                n_cds_hit,
                n_classified,
                n_synonymous,
                n_nonsynonymous
            ) = process_chunk(
                chunk=chunk,

                annotations=annotations,

                cds_sequences=cds_sequences,

                output_handle=out,
            )

            total_processed += (
                n_processed
            )

            total_cds_hit += (
                n_cds_hit
            )

            total_classified += (
                n_classified
            )

            total_synonymous += (
                n_synonymous
            )

            total_nonsynonymous += (
                n_nonsynonymous
            )


            print(
                f"\r[INFO] "
                f"Chunk={chunk_id:,} | "
                f"SNV={total_processed:,} | "
                f"CDS={total_cds_hit:,} | "
                f"classified={total_classified:,} | "
                f"S={total_synonymous:,} | "
                f"N={total_nonsynonymous:,}",
                end="",
                flush=True
            )

    print()



    if total_synonymous > 0:

        n_s_ratio = (
            total_nonsynonymous
            /
            total_synonymous
        )

    else:

        n_s_ratio = None

    print()
    print("=" * 70)

    print("Finished")

    print("=" * 70)

    print(
        f"Total SNVs processed : "
        f"{total_processed:,}"
    )

    print(
        f"CDS overlaps         : "
        f"{total_cds_hit:,}"
    )

    print(
        f"Classified SNVs      : "
        f"{total_classified:,}"
    )

    print(
        f"Synonymous           : "
        f"{total_synonymous:,}"
    )

    print(
        f"Nonsynonymous        : "
        f"{total_nonsynonymous:,}"
    )

    print(
        f"N/S ratio            : "
        f"{n_s_ratio}"
    )

    print()

    print(
        "[OUTPUT]"
    )

    print(
        output_file
    )

    print("=" * 70)

def calc_pnps(group):

    counts = group['mutation_type'].value_counts()
    n_nonsyn = counts.get('nonsynonymous', 0)
    n_syn    = counts.get('synonymous', 0)


    if len(group) == 1 or n_syn == 0:
        return pd.NA

    return round(n_nonsyn / n_syn, 4)




def main():

    parser = argparse.ArgumentParser(
        description=(
            "Annotate SNVs using Bakta CDS and FFN. "
            "SNV reference base is ignored."
        )
    )


    parser.add_argument('-i',"--sp_name", 
        required=True,
        help="Input species"
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="SNV chunk size, default=100000"
    )

    args = parser.parse_args()
    species = args.sp_name

    bakta_file = f"{species}/{species}.tsv.gz"
    ffn_file = f"{species}/{species}.ffn.gz"
    snv_matrix = f"{species}/snv_matrix2base_filter.xls.gz"

    run_pipeline(
        bakta_file=bakta_file,
        ffn_file=ffn_file,
        snv_matrix=snv_matrix,
        output_file=f"{species}_pnps.tsv.gz",
        chunk_size=10000
    )


    df = pd.read_table(f'/{species}_pnps.tsv.gz')
    pnps_series = df.groupby('locus_tag').apply(calc_pnps)
    pnps_series.name = 'gene_pNpS'

    df = df.merge(pnps_series, on='locus_tag', how='left')

    print(df)
    df.to_csv(f'{species}_pnps.tsv.gz', sep='\t', index=False, compression='gzip')




if __name__ == "__main__":

    main()





