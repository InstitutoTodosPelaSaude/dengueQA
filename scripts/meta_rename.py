#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
# seqtree_handler.py -> This code renames (-action=rename) tree tips (-format=tree) and sequence
#                       headers (-format=fasta). given tab delimited file containing
#                       old and new taxa names, one pair per line (split
#                       by a \t character). It also removes (-action=remove) or keeps (-action=keep)
#                       sequences in sequence files or taxa in trees, given a simple .txt file
#                       matching the target sequences/taxa.
#
# Release date: 2020-05-24
# Last update: 2024-063-24

from Bio import Phylo
from Bio import SeqIO
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Rename taxa names in tree and sequence files according to columns from a given metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="FASTA or TREE input file")
    parser.add_argument("--format", required=True, nargs=1, type=str, default='fasta', choices=['fasta', 'tree'], help="File format")
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--columns", required=True, nargs='+', type=str,  help="List of column names, first element being the index")
    parser.add_argument("--separator", required=True, type=str,  help="Separator in new sequence/taxa names")
    parser.add_argument("--output", required=True, help="Renamed output file")
    args = parser.parse_args()
    
    input = args.input
    format = args.format[0]
    metadata = args.metadata
    list_columns = args.columns
    separator = args.separator
    output = args.output
    # output = input.split('.')[0] + '_ren.' + format


#     path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/pangolin/sub1_20210409/tree2/'
#     input = path + "B1526_alignment.fasta.tree"
#     metadata = path + 'metadata_nextstrain.tsv'
#     format = 'tree'
#     list_columns = ['strain', 'date']
#     separator = '|'
#     output = path + 'test.tree'
#     print(list_columns)

    # nextstrain metadata
    print('Loading metadata...')
    df = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)
    columns = [col.strip() for col in list_columns]
    df = df[columns]
    df.fillna('', inplace=True)

    ids = df[columns[0]].to_list()

    # get new names
    replacements = {}
    for id in ids:
        parts = []
        for col in columns:
            # print(col)
            field = df.loc[df[columns[0]] == id, col].values[0]
            parts.append(field)
        # print(parts)
        new_name = separator.join(parts)
        if id not in replacements.keys():
            replacements[id] = new_name
        # entry = id + '\t' + new_name
        # targets.append(entry)
    # print(replacements)

    # rename newick tree files
    if format == 'tree':
        c = 1

        tree = Phylo.read(input, 'newick')
        print('Starting tree file processing...')
        # rename clade names

        for clade in tree.find_clades():
            # for line in targets:
            #     oldName = line.split("\t")[0]
            #     newName = line.split("\t")[1].strip()

            if str(clade.name) in replacements.keys():
                print('Renaming ' + str(clade.name) + ' as ' + replacements[str(clade.name)])
                clade.name = replacements[str(clade.name)]
                c += 1


        Phylo.write([tree], output, 'newick')
        print('\nTree file successfully renamed: ' + output)


    # rename fasta files
    if format == 'fasta':
        print('Starting sequence file processing...')
        fasta_sequences = SeqIO.parse(open(input), 'fasta')

        # rename sequence names
        print('\n### Renaming sequences...')

        # perform simple renaming
        found = []
        not_found = []
        duplicates = []
        c = 1
        with open(output, 'w') as outfile:
            for fasta in fasta_sequences:
                id, seq = fasta.description, fasta.seq
                if id in replacements.keys():
                    if replacements[id] not in found:
                        entry = ">" + replacements[id] + "\n" + str(seq).upper() + "\n"
                        print(str(c) + '. Renamed - ' + replacements[id])
                        outfile.write(entry)
                        found.append(replacements[id])

                        c += 1
                    else:
                        print('Duplicate found: ' + replacements[id])
                        duplicates.append(replacements[id])
                else:
                    not_found.append(id)

        if len(not_found) > 0:
            print('\n### Total of sequences not found = ' + str(len(not_found)))
            for num, record in enumerate(not_found):
                print('\t* ' + str(num + 1) + ' - ' + record)

        if len(duplicates) > 0:
            print('\n### Total of duplicates = ' + str(len(duplicates)))
            for num, record in enumerate(duplicates):
                print('\t* ' + str(num + 1) + ' - ' + record)
        print('\n### A total of ' + str(len(found)) + ' sequences were renamed \n')
