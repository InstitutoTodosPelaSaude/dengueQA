#!/usr/bin/python

# Created by: Anderson Brito
# Release date: 2024-04-06
# Last update: 2024-04-06

import argparse
import pandas as pd
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate the final quality assurance output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--coverage", required=True, help="TSV file with genome coverage information")
    parser.add_argument("--outliers", required=True, help="TSV file with root-to-tip outliers")
    parser.add_argument("--index", required=True, help="Column with unique identifiers, found in both TSV files")
    parser.add_argument("--output", required=True, help="TSV file with quality assurance information")
    args = parser.parse_args()

    coverage = args.coverage
    outliers = args.outliers
    index = args.index
    output = args.output

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


    # Load coverage metadata
    dfC = load_table(coverage)
    dfC.fillna('', inplace=True)
    # print(dfC)

    # Load root-to-tip outliers
    dfR = load_table(outliers)
    dfR.fillna('', inplace=True)
    # print(dfR)

    # merge dataframes
    df = pd.merge(dfC, dfR, on=index, how='left')

    # fix column names
    newcolnames = {'strain': 'sample_id', 'serotype': 'arbo_subtype'}
    df = df.rename(columns=newcolnames)

    # Set quality for genomes not informed by root-to-tip file
    df['seq_quality'] = df['seq_quality'].fillna('Yes')
    # print(df)

    # Assess completeness
    df['comp_genome'] = np.where(df['genome_coverage'].astype(float) >= 0.70, 'Yes', 'No')
    df['comp_CprME'] = np.where(df['C-prM-E'].astype(float) >= 0.95, 'Yes', 'No')
    df['comp_E'] = np.where(df['E'].astype(float) >= 0.95, 'Yes', 'No')

    # assess what sequence should be sequenced
    df['submission'] = np.where(df['seq_quality'] == 'No', 'No',
                        np.where(df['comp_genome'] == 'Yes', 'genome',
                        np.where(df['comp_CprME'] == 'Yes', 'CprME',
                        np.where(df['comp_E'] == 'Yes', 'E', 'No'))))
    
    # reorder columns
    df = df[['sample_id', 'arbo_subtype', 'seq_quality', 'genome_coverage', 'C', 'prM', 'E', 'NS1', 'NS2A', 'NS2B', 'NS3', 'NS4A', '2K', 'NS4B', 'NS5', 'C-prM-E', 'comp_genome', 'comp_CprME', 'comp_E', 'submission']]
    # print(df)

    # Write quality assurance matrix
    df.to_csv(output, sep='\t', index=False)
