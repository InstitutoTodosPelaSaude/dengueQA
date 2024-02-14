#!/usr/bin/python

# Created by: Anderson Brito
# Release date: 2024-01-30
# Last update: 2024-01-30

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
    parser.add_argument("--references", required=True, help="FASTA file with reference genomes")
    parser.add_argument("--evalue", required=False, default=1e-3, help="Maximum e-value threshold in scientific notation (e.g. 1e-10)")
    parser.add_argument("--output", required=True, help="TSV file with descriptions about the top hits")
    args = parser.parse_args()

    query = args.query
    references = args.references
    query_dir, qFile = os.path.split(query)
    database_dir, dFile = os.path.split(references)
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

    # # return to root directory
    # os.chdir(root_dir)
    blastn_cline = NcbiblastnCommandline(query=query_dir + qFile, db=database_dir + dFile, evalue=evalue, outfmt=5, num_threads=6)


    result = os.popen(str(blastn_cline))
    # result = subprocess.run(str(blastn_cline), shell=True, capture_output=True, text=True)

    # Parse the XML result
    blast_records = NCBIXML.parse(result)

    # Open a file for writing
    # output_file_path = path + "output/results.tsv"
    with open(root_dir + outfile, 'w') as output_file:
        # Write the header
        output_file.write("itps_id\tarbo_subtype\treference\n")

        # Process each query result
        for blast_record in blast_records:
            query_id = blast_record.query # Name of the query sequence

            if len(blast_record.alignments) > 0:
                top_hit = blast_record.alignments[0]  # Assuming the first alignment is the top hit
                hit_description = top_hit.hit_def
                # print(hit_description)
                # print(top_hit.title)
                seqname = hit_description

                if float(top_hit.hsps[0].expect) < evalue:
                    # Extract arbo_subtype from hit_description
                    arbo_subtype = hit_description.split('|')[1].replace("Dengue virus ", "DenV")
                else:
                    arbo_subtype = 'None'
                    seqname = 'None'
            else:
                arbo_subtype = 'None'
                seqname = 'None'

        # Write to the file
            output_file.write(f"{query_id}\t{arbo_subtype}\t{seqname}\n")

    os.system("rm %s/%s.*" % (query_dir, qFile))
    os.system("rm %s/%s.*" % (database_dir, dFile))

    print(f"Results written to {outfile}")



