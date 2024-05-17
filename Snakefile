import os

SAMPLES = [
	"DenV1",
	"DenV2",
	"DenV3",
	"DenV4"]

rule all:
	input:
		qa_matrix = "output/matrix_dengue_qa.tsv",
		coverage_matrix = "output/coverage/matrix_coverage.tsv",
		outliers_matrix = "output/root2tip/matrix_outliers.tsv",
		trees1 = expand("output/tree/{sample}_masked.tree", sample=SAMPLES),
		trees2 = expand("output/tree/{sample}_masked_renamed.tree", sample=SAMPLES)


# Define file names
rule files:
	params:
		new_genomes = "input/new_genomes.fasta",
		orig_metadata = "input/metadata.tsv",
		references = "references/fasta/reference_genomes.fasta",
		ref_typing = "references/fasta/reference_typing_s3.fasta",
		ref_table = "references/fasta/reference_table.tsv",
		blast_results = "output/blast/blast_results.tsv",
		metadata_typing = "references/metadata/table_ncbi_s3.tsv",
		ref_dataset = "references/fasta/DenV_genomes_ncbi.fasta",
files = rules.files.params


# Define parameters
rule parameters:
	params:
		evalue = 1e-10,
		bootstrap = 1, # default = 1, but ideally it should be >= 100
		model = "GTR",
		clock_filter = 10,
		clock_rate = 0.0007,
		clock_std_dev = 0.00004,
		root = "least-squares"
parameters = rules.parameters.params


rule options:
	params:
		threads = 4
options = rules.options.params



rule blast:
	message:
		"""
		Perform blast searches for serotyping, genotyping and assigning lineages
		"""
	input:
		new_genomes = files.new_genomes,
		ref_typing = files.ref_typing
	params:
		evalue = parameters.evalue,
		id = "strain",
		fields = "accno serotype genotype major_lineage minor_lineage"
	output:
		blast = files.blast_results
	shell:
		"""
		python scripts/localBlastDNA.py \
			--query {input.new_genomes} \
			--id {params.id} \
			--subjects {input.ref_typing} \
			--fields {params.fields} \
			--separator '|' \
			--evalue {params.evalue} \
			--output {output.blast}
		"""


rule split:
	message:
		"""
		Split datasets based on viral types (dengue serotypes)
		"""
	input:
		new_genomes = files.new_genomes,
		references = files.references,
		reftable = files.ref_table,
		blast = rules.blast.output.blast
	params:
		refcol = "reference",
		seqcol = "strain",
		typecol = "serotype",
		target = lambda wildcards: wildcards.sample
	output:
		typeseq = "output/sequences/{sample}_newsequences.fasta",
		refseq = "output/sequences/{sample}_reference.fasta"
	shell:
		"""
		python scripts/dataset_split.py \
			--query {input.new_genomes} \
			--references {input.references} \
			--reftable {input.reftable} \
			--refcol {params.refcol} \
			--target {params.target} \
			--blast {input.blast} \
			--seqcol {params.seqcol} \
			--typecol {params.typecol} \
			--output1 {output.typeseq} \
			--output2 {output.refseq}
		"""


rule align_reference:
	message:
		"""
		Align newly sequenced genomes with their reference genomes
		"""
	input:
		sequences = rules.split.output.typeseq,
		reference = rules.split.output.refseq
	params:
		threads = options.threads,
	output:
		alignment = "output/alignment/{sample}_alignref.fasta"
	shell:
		"""
		if [ -s {input.sequences} ]; then
			augur align \
				--sequences {input.sequences} \
				--reference-sequence {input.reference} \
				--nthreads {params.threads} \
				--output {output.alignment}
		else
			echo "The input sequences file {input.sequences} is empty, skipping alignment."
			touch {output.alignment}
		fi
		"""


rule coverage:
	message:
		"""
		Calculate the coverage per CDS sequences and whole genomes
		"""
	input:
		alignment = rules.align_reference.output.alignment,
		reftable = files.ref_table,
		gff_file = "references/annotations/{sample}_sequence.gff"
	params:
		refcol = "reference",
		seqcol = "strain",
		typecol = "serotype",
		target = lambda wildcards: wildcards.sample
	output:
		coverage = "output/coverage/{sample}_coverage.tsv"
	shell:
		"""
		if [ -s {input.alignment} ]; then
			python scripts/genome_coverage.py \
				--input {input.alignment} \
				--reftable {input.reftable} \
				--gff {input.gff_file} \
				--target {params.target} \
				--seqcol {params.seqcol} \
				--typecol {params.typecol} \
				--refcol {params.refcol} \
				--output {output.coverage}
		else
			echo "The input alignment file {input.alignment} is empty, skipping coverage assessment."
			touch {output.coverage}
		fi
		"""


rule selection:
	message:
		"""
		Select contextual genomes for root-to-tip analyses
		"""
	input:
		metadata = files.metadata_typing,
		sampling = "references/sampling/{sample}_sampling.tsv",
		coverage = rules.coverage.output.coverage,
	params:
		index = "strain",
		unit = "month"
	output:
		report = "output/sequences/{sample}_stats.txt",
		list = "output/sequences/{sample}_list.txt"
	shell:
		"""
		if [ -s {input.coverage} ]; then
			python scripts/genome_selector.py \
				--metadata {input.metadata} \
				--accno-column {params.index} \
				--time-unit {params.unit} \
				--scheme {input.sampling} \
				--report {output.report} \
				--output1 {output.list}
		else
			echo "The input coverage file {input.coverage} is empty, skipping contextual genome selection."
			touch {output.report}
			touch {output.list}
		fi
		"""


rule dataset:
	message:
		"""
		Joint newly sequenced genomes and contextual genomes in a single dataset
		"""
	input:
		genomes = files.ref_dataset,
		new_genomes = rules.split.output.typeseq,
		list = rules.selection.output.list
	output:
		dataset = "output/sequences/{sample}_dataset.fasta"
	shell:
		"""
		if [ -s {input.list} ]; then
			python scripts/add_new_sequences.py \
				--genomes {input.genomes} \
				--new-genomes {input.new_genomes} \
				--keep {input.list} \
				--output {output.dataset}
		else
			echo "The input list file {input.list} is empty, skipping contextual genome filtering."
			touch {output.dataset}
		fi
		"""



rule fix_metadata:
	message:
		"""
		Add missing columns in new sequence metadata
		"""
	input:
		metadata = files.orig_metadata,
		blast = files.blast_results
	params:
		index = "strain",
		targets = "serotype#11, genotype#12, major_lineage#13, minor_lineage#14",
	output:
		temp_metadata = "output/metadata/temp_metadata.tsv",
		new_metadata = "output/metadata/DenV_itps_metadata.tsv",
	shell:
		"""
		sed 's/sample_id/{params.index}/g' {input.metadata} > {output.temp_metadata}
		
		python scripts/reformat_dataframe.py \
			--input1 {output.temp_metadata} \
			--input2 {input.blast} \
			--index {params.index} \
			--action add \
			--mode columns \
			--targets \"{params.targets}\" \
			--output {output.new_metadata}
		"""


rule final_dataset:
	message:
		"""
		
		"""
	input:
		dataset = rules.dataset.output.dataset,
		metadata1 = files.metadata_typing,
		metadata2 = rules.fix_metadata.output.new_metadata
	output:
		final_metadata = "output/metadata/{sample}_metadata.tsv",
		final_genomes = "output/sequences/{sample}_genomes.tsv",
		rename_file = "output/tree/{sample}_rename.tsv"
	shell:
		"""
		if [ -s {input.dataset} ]; then
		python scripts/process_metadata.py \
			--sequences {input.dataset} \
			--metadata1 {input.metadata1} \
			--metadata2 {input.metadata2} \
			--output1 {output.final_metadata} \
			--output2 {output.final_genomes} \
			--output3 {output.rename_file}
		else
			echo "The input dataset file {input.dataset} is empty, skipping final dataset creation."
			touch {output.final_metadata}
			touch {output.final_genomes}
			touch {output.rename_file}
		fi
		"""


rule aligndb:
	message:
		"""
		Align new genomes and contextual genomes for phylogenetics
		"""
	input:
		genomes = rules.final_dataset.output.final_genomes,
		existing = rules.align_reference.output.alignment,
		reference = rules.split.output.refseq,
		reftable = files.ref_table,
	params:
		threads = options.threads,
		type = "{sample}"
	output:
		refname = "output/alignment/{sample}_reflist.txt",
		msa_noref = "output/alignment/{sample}_existing.fasta",
		msa_file = "output/alignment/{sample}_aligned.fasta"
	shell:
		"""
		if [ -s {input.existing} ]; then
			grep {params.type} {input.reftable} | cut -d$'\t' -f 2 > {output.refname}

			python scripts/masterkey.py \
				--input {input.existing} \
				--format fasta \
				--action remove \
				--list {output.refname} \
				--output {output.msa_noref}

			augur align \
				--sequences {input.genomes} \
				--existing-alignment {output.msa_noref} \
				--reference-sequence {input.reference} \
				--nthreads {params.threads} \
				--output {output.msa_file} \
				--remove-reference
		else
			echo "The input sequences file {input.existing} is empty, skipping alignment."
			touch {output.refname}
			touch {output.msa_noref}
			touch {output.msa_file}
		fi
		"""


def denv_utrs(sample):
	params_dict = {
		"DenV1": (94, 462),
		"DenV2": (96, 451),
		"DenV3": (94, 440),
		"DenV4": (101, 384)
	}
	return params_dict.get(sample, (None, None))

rule mask:
	message:
		"""
		Mask non-coding 5' and 3' regions
		"""
	input:
		sequences = rules.aligndb.output.msa_file
	params:
		start=lambda wildcards: denv_utrs(wildcards.sample)[0],
		end=lambda wildcards: denv_utrs(wildcards.sample)[1],
	output:
		masked = "output/alignment/{sample}_aligned_masked.fasta"
	shell:
		"""
		if [ -s {input.sequences} ]; then
			augur mask \
				--sequences {input.sequences} \
				--mask-from-beginning {params.start} \
				--mask-from-end {params.end} \
				--output {output.masked}
		else
			echo "The input sequences file {input.sequences} is empty, skipping masking step."
			touch {output.masked}
		fi
		"""


rule tree:
	message: "Building phylogenetic tree"
	input:
		alignment = rules.mask.output.masked
	params:
		threads = options.threads,
		model = parameters.model,
		bootstrap = parameters.bootstrap,
	output:
		tree = "output/tree/{sample}_masked.tree"
	shell:
		"""
		if [ -s {input.alignment} ]; then
			iqtree \
				-s {input.alignment} \
				-b {params.bootstrap} \
				-nt {params.threads} \
				-m {params.model}

			mv {input.alignment}.treefile {output.tree}
		else
			echo "The input file {input.alignment} is empty, skipping phylogenetic analysis."
			touch {output.tree}
		fi
		"""


rule root2tip:
	message:
		"""
		Perform root-to-tip analysis:
		  - detect molecular clock outliers
		  - generate root-to-tip regression plot
		"""
	input:
		tree = rules.tree.output.tree,
		alignment = rules.mask.output.masked,
		metadata = rules.final_dataset.output.final_metadata
	params:
		clock_filter = parameters.clock_filter,
		clock_rate = parameters.clock_rate,
		clock_std_dev = parameters.clock_std_dev,
		root = parameters.root,
		outdir = "./output/root2tip"
	output:
		rtt_plot = "output/root2tip/{sample}_rtt_plot.pdf",
		dates = "output/root2tip/{sample}_dates.txt",
		outliers = "output/root2tip/{sample}_outliers.tsv"
	shell:
		"""
		if [ -s {input.tree} ]; then
			treetime \
				--tree {input.tree} \
				--aln {input.alignment} \
				--dates {input.metadata} \
				--clock-filter {params.clock_filter} \
				--reroot {params.root} \
				--gtr JC69 \
				--clock-rate {params.clock_rate} \
				--clock-std-dev {params.clock_std_dev} \
				--max-iter 2 \
				--coalescent skyline \
				--plot-rtt rtt_plot.pdf \
				--tip-labels \
				--verbose 1 \
				--outdir {params.outdir}
			
			mv {params.outdir}/rtt_plot.pdf {output.rtt_plot}
			mv {params.outdir}/dates.tsv {output.dates}
			echo "--" >> {output.dates}
			grep "\-\-" {output.dates} | grep -v NODE | cut -d$'\t' -f 1 > {output.outliers}
			sed -i '' -e $'1s/^/strain\tseq_quality\\\n/' -e $'1,$s/$/\tNo/' {output.outliers}
		else
			echo "The input file {input.tree} is empty, skipping root-to-tip analysis."
			touch {output.dates}
			touch {output.outliers}
			touch {output.rtt_plot}
		fi
		"""


rule assurance:
	input:
		coverage_files = expand("output/coverage/{sample}_coverage.tsv", sample=SAMPLES),
		outlier_files = expand("output/root2tip/{sample}_outliers.tsv", sample=SAMPLES)
	params:
		index = "strain",
	output:
		coverage = "output/coverage/matrix_coverage.tsv",
		outliers = "output/root2tip/matrix_outliers.tsv",
		matrix = "output/matrix_dengue_qa.tsv"
	shell:
		"""
		# Concatenate coverage files
		(head -n 1 {input.coverage_files[0]} && tail -n +2 -q {input.coverage_files}) > {output.coverage}
		
		# Concatenate outlier files
		(head -n 1 {input.outlier_files[0]} && tail -n +2 -q {input.outlier_files}) > {output.outliers}

		# Assuming `qamatrix.py` takes the concatenated files to produce the final matrix
		python scripts/qamatrix.py \
			--coverage {output.coverage} \
			--outliers {output.outliers} \
			--index {params.index} \
			--output {output.matrix}
		"""


rule tempest:
	message:
		"""
		Inspect root-to-tip results using TemPest
		"""
	input:
		tree = rules.tree.output.tree,
		metadata = rules.final_dataset.output.final_metadata
	output:
		tree = "output/tree/{sample}_masked_renamed.tree"
	shell:
		"""
		if [ -s {input.tree} ]; then
		python scripts/meta_rename.py \
			--input {input.tree} \
			--format tree \
			--metadata {input.metadata} \
			--column strain country_code date \
			--separator "|" \
			--output {output.tree}
		else
			echo "The input file {input.tree} is empty, skipping TempEst analysis."
			touch {output.tree}
		fi
		"""

### Clearing the working directory (only executed when needed)

rule clean:
	message: "Removing directories: {params}"
	params:
		"output"

	shell:
		"""
		rm -rfv {params}
		"""



#rule xxx:
#	message:
#		"""
#		
#		"""
#	input:
#		metadata = arguments.
#	params:
#		index = arguments.,
#		date = arguments.
#	output:
#		matrix = "results/"
#	shell:
#		"""
#		python scripts/ \
#			--metadata {input.} \
#			--index-column {params.} \
#			--extra-columns {params.} \
#			--date-column {params.} \
#			--output {output.}
# 		
#		cp results/combined.tsv data/combined_cache.tsv
#		"""

