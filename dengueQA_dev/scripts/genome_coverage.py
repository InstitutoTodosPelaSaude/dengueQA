#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
#
# Release date: 16/02/2024
# Last update: 16/02/2024
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from Bio import SeqIO
import os
import pandas as pd
import argparse

pd.set_option('display.max_columns', 500)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assess the level of genome completeness after sequencing, by CDS and full length",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="FASTA multiple sequence alingment (MSA)")
    parser.add_argument("--blast", required=True, help="TSV file with descriptions about the top hits of the blast search")
    parser.add_argument("--gff", required=True, help="Annotation file in GFF format")
    parser.add_argument("--target", required=True, help="Name of the target virus type (e.g. DenV1).")
    parser.add_argument("--typecol", required=True, help="Column in blast results that contains virus type information")
    parser.add_argument("--seqcol", required=True, help="Column in blast results that contains query sequence names")
    parser.add_argument("--refcol", required=True, help="Column in blast results that contains reference sequence names")
    parser.add_argument("--output", required=True, help="TSV file with genome completeness with respect to reference")
    args = parser.parse_args()

    input_file = args.input
    blast_results = args.blast
    target = args.target
    type_col = args.typecol
    ref_col = args.refcol
    seq_col = args.seqcol
    annotation = args.gff
    output_file = args.output


    # path = "/Users/Anderson/Documents/github/dengueQA/dengueQA_dev/"
    # target = 'DenV4'
    # input_file = path + "output/alignments/" + target + "_alignment.fasta"
    # seq_col = 'itps_id'
    # ref_col = 'reference'
    # type_col = "arbo_subtype"
    # annotation = path + "references/annotations/" + target + "_sequence.gff"
    # blast_results = path + 'output/blast/blast_result.tsv'
    # output_file = path + "output/coverage/" + target + "_coverage.tsv"



    """ BLAST RESULTS """
    blast_data = pd.read_csv(blast_results, sep='\t')
    outdf = blast_data[blast_data[type_col] == target]
    refname = outdf[ref_col].tolist()[0]
    outdf = outdf.drop(columns=[ref_col])
    outdf.insert(2, 'genome_coverage', '')


    """ ANNOTATIONS """
    # Read the GFF file into a pandas DataFrame
    gff_data = pd.read_csv(annotation, sep='\t', comment='#')
    gff_data = gff_data[['start', 'end', 'attribute']]
    gff_data['attribute'] = gff_data['attribute'].apply(lambda x: x.split('\"')[1]).replace('\"', '')

    """ GENOME PROCESSING """
    # open fasta alignment as dict
    record_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"), key_function=lambda x: x.description)
    for index, row in gff_data.iterrows():
        cds = row['attribute']
        start = row['start']
        end = row['end']
        outdf[cds] = ''

        # Check if the folder exists
        if not os.path.exists(os.getcwd() + '/' + target + '_alignments/'):
            # If it doesn't exist, create the folder
            os.makedirs(os.getcwd() + '/' + target + '_alignments/')

        cds_alignment = open(os.getcwd() + '/' + target + '_alignments/' + target + '_' + cds + '.fasta', 'w')
        # cds_alignment = open(input_file.replace('.fasta', '_' + cds + '.fasta'), 'w')

        for header in record_dict:
            sequence = str(record_dict[header].seq)
            seq = sequence[start - 1:end]
            cds_length = end - start + 1
            entry = '>' + header + '|' + cds + '\n' + seq + '\n'
            cds_alignment.write(entry)

            if header != refname:
                cds_size = len(str(seq).replace('N', '').replace('-', ''))
                cds_coverage = str(round(cds_size / cds_length, 2))
                outdf.loc[outdf[seq_col] == header, cds] = cds_coverage # add values to columns

                ref_size = len(str(record_dict[refname].seq))
                genome_coverage = str(round(len(sequence.replace('N', '').replace('-', '')) / ref_size, 2))
                outdf.loc[outdf[seq_col] == header, 'genome_coverage'] = genome_coverage # add values to columns


    outdf.to_csv(output_file, sep='\t', index=False)
    print('\nData successfully saved in:\n%s\n' % output_file)

