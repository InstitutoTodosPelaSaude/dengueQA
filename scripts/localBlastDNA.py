#!/usr/bin/python

# Created by: Anderson Brito
# Release date: 2024-01-30
# Last update: 2024-03-24

import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Performs blast searches of new sequences against a local database with reference sequences. As a result, descriptions about the top hits are exported in TSV format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--query", required=True, help="FASTA file with newly sequenced genomes")
    parser.add_argument("--id", required=True, default='id', type=str,  help="Query sequence identifier")
    parser.add_argument("--subjects", required=True, help="FASTA file with subject sequences")
    parser.add_argument("--fields", required=False, default=['match'], nargs='+', type=str, help="List of fields whose values are embedded in the subject's fasta headers, where the first field must be a sequence identifier (e.g. >OR167152|1V_E|Unassigned)")
    parser.add_argument("--separator", required=False, default='|', help="Field separator in subject's fasta headers")
    parser.add_argument("--evalue", required=False, default=1e-3, help="Maximum e-value threshold in scientific notation (e.g. 1e-10)")
    parser.add_argument("--output", required=True, help="TSV file with descriptions of Blast top hits")
    args = parser.parse_args()

    query = args.query
    query_dir, qFile = os.path.split(query)
    query_col = args.id
    subjects = args.subjects
    database_dir, dFile = os.path.split(subjects)
    fields = args.fields
    separator = args.separator
    evalue = float(args.evalue)
    outfile = args.output


    # path = "/Users/Anderson/Documents/github/dengueQA/dengueQA_dev/"
    # query_dir = path + "input/"
    # database_dir = path + "reference/fasta/"
    # qFile = "genomes_test.fasta"
    # dFile = "reference_genomes.fasta"
    # evalue = float(1e-10)

    # working directories
    root_dir = os.getcwd() + '/'
    query_dir = root_dir + query_dir + '/'
    database_dir = root_dir + database_dir + '/'

    if os.path.isfile(query_dir + qFile + '.nin'):
        print("Proceeding with pre-processed Blast databases\n")
        pass
    else:
        os.chdir(database_dir)
        os.system("makeblastdb -in %s -dbtype nucl" % (dFile))

        # os.chdir(query_dir)
        os.chdir(query_dir)
        os.system("makeblastdb -in %s -dbtype nucl" % (qFile))

    blastn_cline = NcbiblastnCommandline(query=query_dir + qFile, db=database_dir + dFile, evalue=evalue, outfmt=5, num_threads=6)

    result = os.popen(str(blastn_cline))

    # Parse the XML result
    blast_records = NCBIXML.parse(result)

    # Open a file for writing
    # output_file_path = path + "output/results.tsv"
    with open(root_dir + outfile, 'w') as output_file:
        # Write the header
        if fields not in ['', None]:
            list_fields = [query_col] + [col.strip() for col in fields]
        else:
            list_fields = [query_col]

        output_file.write('\t'.join(list_fields) + '\n')

        # Process each query result
        typing_data = ['None'] * len(fields)
        subject_name = 'None'
        for blast_record in blast_records:
            query_id = blast_record.query # Name of the query sequence

            if len(blast_record.alignments) > 0:
                top_hit = blast_record.alignments[0]  # Assuming the first alignment is the top hit
                hit_description = top_hit.hit_def
                # print(hit_description)
                # print(top_hit.title)
                subject_name = hit_description
                
                if float(top_hit.hsps[0].expect) < evalue:
                    # Extract arbo_subtype from hit_description
                    typing_data = hit_description.split(separator)
                    # print(typing_data)

            # Write to the file
            output_file.write('\t'.join([query_id] + typing_data) + '\n')


    os.system("rm %s/%s.*" % (query_dir, qFile))
    os.system("rm %s/%s.*" % (database_dir, dFile))

    print(f"Results written to {outfile}")



