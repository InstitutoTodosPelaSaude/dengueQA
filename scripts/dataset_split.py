#!/usr/bin/python

# Created by: Anderson Brito
# Release date: 2024-01-30
# Last update: 2024-01-30

import argparse
from Bio import SeqIO
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Splits sequence datasets according to a target virus type",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--query", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--references", required=True, help="FASTA file with reference genomes")
    parser.add_argument("--reftable", required=True, help="Tabulated file with list of virus types and their reference genomes")
    parser.add_argument("--refcol", required=True, help="Column in reference table that contains reference sequence names")
    parser.add_argument("--target", required=True, help="Name of the target virus type (e.g. DenV1) displayed in the type column.")
    parser.add_argument("--blast", required=True, help="TSV file with descriptions about the top hits of the blast search")
    parser.add_argument("--seqcol", required=True, help="Column in blast results that contains query sequence names")
    parser.add_argument("--typecol", required=True, help="Column in blast results that contains virus type information")
    parser.add_argument("--output1", required=True, help="FASTA file containing newly sequenced genomes belonging to the target type")
    parser.add_argument("--output2", required=True, help="FASTA file containing reference genomes belonging to the target type")
    args = parser.parse_args()

    genomes = args.query
    references = args.references
    reftable = args.reftable
    blast_result = args.blast
    seq_col = args.seqcol
    ref_col = args.refcol
    type_col = args.typecol
    target = args.target
    outfile1 = args.output1
    outfile2 = args.output2


    # path = "/Users/anderson/Documents/github/dengueQA/dengueQA_dev/"
    # genomes = path + "input/genomes_test.fasta"
    # references = path + "references/fasta/reference_genomes.fasta"
    # blast_result = path + 'output/blast_result.tsv'
    # seq_col = 'itps_id'
    # ref_col = 'reference'
    # type_col = "arbo_subtype"
    # target = 'DenV4'
    # outfile1 = path + 'output/' + target + "_sequences.fasta"
    # outfile2 = path + 'output/' + target + "_references.fasta"


    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    # Load sample metadata
    df1 = load_table(blast_result)
    df1 = df1[df1[type_col] == target] # filter sequence type set from dataframe
    df1.fillna('', inplace=True)

    # Load reference table
    df2 = load_table(reftable)
    df = pd.merge(df1, df2, on=type_col, how='left')

    with open(outfile1, "w") as output_handle:
        for record in SeqIO.parse(genomes, "fasta"):
            if record.description in df[seq_col].unique().tolist():
                SeqIO.write(record, output_handle, "fasta")

    with open(outfile2, "w") as output_handle:
        for record in SeqIO.parse(references, "fasta"):
            if record.description in df[ref_col].unique().tolist():
                SeqIO.write(record, output_handle, "fasta")
