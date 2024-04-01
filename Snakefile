SAMPLES = [
	"DenV1",
	"DenV2",
	"DenV3",
	"DenV4"]


rule all:
	message:
		"""
		Execute all rules using terminal commands in a logical order.
		"""
	shell:
		"""
		snakemake --cores all test_results_go
		"""

# Define file names
rule files:
	params:
		new_genomes = "input/new_genomes.fasta",
		references = "references/fasta/reference_genomes.fasta",
		ref_typing = "references/fasta/reference_typing_s2.fasta",
		ref_table = "references/fasta/reference_table.tsv",
		blast_results = "output/blast/blast_results.tsv",
		metadata_typing = "references/metadata/table_ncbi_s2.tsv",
		orig_metadata = "input/pt_BR_all_itps_dengue-0-20240329000013.tsv",
		new_metadata = "input/DenV_itps_metadata.tsv",
		ref_dataset = "references/fasta/DenV_genomes_ncbi.fasta",


		# sequence_dataset = "input_files/gisaid_hcov-19.fasta",
		# metadata_gisaid = "input_files/metadata_nextstrain.tsv",
		# corelab_metadata = "input_files/metadata_corelab.tsv",
		# sample_metadata = "input_files/impacc-virology-clin-sample.csv",
		# patient_metadata = "input_files/impacc-virology-clin-individ.csv",
		# batch_layout = "input_files/batch_layout.csv",
		# reference = "config/reference_seq.fasta",
		# annotation = "config/sequence.gb",
		# genemap = "config/genemap.gff",
		# refgenome_size = "29903",
		# max_missing = "30"
files = rules.files.params


# Define parameters
rule parameters:
	params:
		evalue = 1e-10,
		mask_5prime = 142,
		mask_3prime = 548,
		bootstrap = 1, # default = 1, but ideally it should be >= 100
		model = "GTR",
		root = "JF912185",
		clock_rate = 0.0003,
		clock_std_dev = 0.0001,
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


rule split_go:
	input:
		expand("output/sequences/{sample}_newsequences.fasta", sample=SAMPLES),

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



rule alignref_go:
	input:
		expand("output/alignment/{sample}_alignref.fasta", sample=SAMPLES),

rule align_reference:
	message:
		"""
		Align newly sequenced genomes with there reference genomes
		"""
	input:
		sequences = rules.split.output.typeseq,
		reference = rules.split.output.refseq
	params:
		threads = options.threads
	output:
		alignment = "output/alignment/{sample}_alignref.fasta"
	shell:
		"""
		augur align \
			--sequences {input.sequences} \
			--reference-sequence {input.reference} \
			--nthreads {params.threads} \
			--output {output.alignment}
		"""



rule coverage_go:
	input:
		expand("output/coverage/{sample}_coverage.tsv", sample=SAMPLES),

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
		python scripts/genome_coverage.py \
			--input {input.alignment} \
			--reftable {input.reftable} \
			--gff {input.gff_file} \
			--target {params.target} \
			--seqcol {params.seqcol} \
			--typecol {params.typecol} \
			--refcol {params.refcol} \
			--output {output.coverage}
		"""



rule selection_go:
	input:
		expand("output/sequences/{sample}_list.txt", sample=SAMPLES),

rule selection:
	message:
		"""
		Select contextual genomes for root-to-tip analyses
		"""
	input:
		metadata = files.metadata_typing,
		sampling = "references/metadata/{sample}_sampling.tsv"
	params:
		index = "strain"
	output:
		report = "output/sequences/{sample}_stats.txt",
		list = "output/sequences/{sample}_list.txt"
	shell:
		"""		
		python scripts/genome_selector.py \
			--metadata {input.metadata} \
			--accno-column {params.index} \
			--scheme {input.sampling} \
			--report {output.report} \
			--output1 {output.list}
		"""



rule dataset_go:
	input:
		expand("output/sequences/{sample}_dataset.fasta", sample=SAMPLES),

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
		python add_new_sequences.py \
			--genomes {input.genomes} \
			--new-genomes {input.new_genomes} \
			--keep {input.list} \
			--output {output.dataset}
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
		targets = "serotype#21, genotype#22, major_lineage#23, minor_lineage#24",
	output:
		new_metadata = files.new_metadata
	shell:
		"""
		python reformat_dataframe.py \
			--input1 {input.metadata} \
			--input2 {input.blast} \
			--index {params.index} \
			--action add \
			--mode columns \
			--targets {params.targets} \
			--output {output.new_metadata}
		"""



# rule lineages_rep:
# 	message:
# 		"""
# 		Getting list of representative genomes per pangolin lineage:
# 		- Pick one representative genome per lineage
# 		- Export list with sequence names
# 		"""
# 	params:
# 		download_file = "yes",
# 	output:
# 		list = "output_files/sequences/rep_genomes.txt",
# 	shell:
# 		"""
# 		python.nextstrain scripts/get_lineage_reps.py \
# 			--download {params.download_file} \
# 			--output {output.list}
		
# 		echo Wuhan/Hu-1/2019 >> {output.list}
# 		echo Wuhan/WH01/2019 >> {output.list}
# 		"""



# rule base_dataset:
# 	message:
# 		"""
# 		Export base dataset with genomes and metadata
# 		"""
# 	input:
# 		genomes = "input_files/gisaid_hcov-19.fasta",
# 		metadata = "input_files/metadata_nextstrain.tsv",
# 		list = "output_files/sequences/rep_genomes.txt"
# 	output:
# 		base_dataset = "output_files/sequences/base_dataset.fasta",
# 		base_metadata = "output_files/metadata/base_metadata.tsv"
# 	shell:
# 		"""
	
# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.genomes} \
# 			--format fasta \
# 			--action keep \
# 			--list {input.list} \
# 			--output {output.base_dataset}
		
# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.metadata} \
# 			--format tsv \
# 			--action keep \
# 			--index strain \
# 			--list {input.list} \
# 			--output {output.base_metadata}
# 		"""


# rule filter_coverage:
# 	message:
# 		"""
# 		Filtering sequence files to:
# 		- Identify and remove poor quality genomes
# 		- Generate initial quality assurance matrix
# 		"""
# 	input:
# 		genomes = files.new_genomes
# 	params:
# 		size = files.refgenome_size,
# 		index = "sample_id",
# 		missing = files.max_missing
# 	output:
# 		matrix = "output_files/qa/qa_matrix1.tsv",
# 		rename = "output_files/sequences/rename.tsv",
# 		new_sequences = "output_files/sequences/renamed_genomes.fasta"
# 	shell:
# 		"""
# 		python.nextstrain scripts/filter_by_coverage.py \
# 			--genomes {input.genomes} \
# 			--index {params.index} \
# 			--refgenome-size {params.size} \
# 			--max-missing {params.missing} \
# 			--output {output.matrix} \
# 			--output2 {output.rename}
			
		
# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.genomes} \
# 			--format fasta \
# 			--action rename \
# 			--list {output.rename} \
# 			--output {output.new_sequences}	
# 		"""



# rule inspect_metadata:
# 	message:
# 		"""
# 		Inspecting metadata file to:
# 		- Detect genomes with missing metadata
# 		- Expand quality assurance matrix
# 		"""
# 	input:
# 		metadata1 = files.corelab_metadata,
# 		metadata2 = files.sample_metadata,
# 		metadata3 = files.patient_metadata,
# 		batch = files.batch_layout,
# 		matrix1 = rules.filter_coverage.output.matrix
# 	params:
# 		index = "sample_id"
# 	output:
# 		metadata = "output_files/metadata/metadata1.tsv",
# 		matrix = "output_files/qa/qa_matrix2.tsv",
# 		rename = "output_files/sequences/rename2.tsv"
# 	shell:
# 		"""
# 		python.nextstrain scripts/inspect_metadata.py \
# 			--metadata1 {input.metadata1} \
# 			--metadata2 {input.metadata2} \
# 			--metadata3 {input.metadata3} \
# 			--batch {input.batch} \
# 			--index {params.index} \
# 			--matrix {input.matrix1} \
# 			--output1 {output.metadata} \
# 			--output2 {output.matrix} \
# 			--output3 {output.rename}
# 		"""



# rule multifasta:
# 	message:
# 		"""
# 		Combining sequence files as a multifasta file
# 		"""
# 	input:
# 		base_dataset = "output_files/sequences/base_dataset.fasta",
# 		new_genomes = "output_files/sequences/renamed_genomes.fasta",
# 		qamatrix = "output_files/qa/qa_matrix2.tsv"
# 	params:
# 		format = "fasta"
# 	output:
# 		list_seqs = "output_files/sequences/full_genomes_list.txt",
# 		filtered_seqs = "output_files/sequences/filtered_seqs.fasta",
# 		combined_seqs = "output_files/sequences/quality_sequences.fasta"
# 	shell:
# 		"""
# 		grep -v FAIL {input.qamatrix} | cut -d$'\t' -f 1 | sed -e 1d > {output.list_seqs}

# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.new_genomes} \
# 			--format {params.format} \
# 			--action keep \
# 			--list {output.list_seqs} \
# 			--output {output.filtered_seqs}

# 		cat {input.base_dataset} {output.filtered_seqs} > {output.combined_seqs}
# 		"""



# rule combine_metadata:
# 	message:
# 		"""
# 		Combining metadata files
# 		"""
# 	input:
# 		base_metadata = "output_files/metadata/base_metadata.tsv",
# 		sample_metadata = rules.inspect_metadata.output.metadata
# 	params:
# 		index = "strain"
# 	output:
# 		combined_metadata = "output_files/metadata/combined_metadata.tsv"
# 	shell:
# 		"""
# 		python.nextstrain scripts/metadata_merger.py \
# 			--metadata1 {input.base_metadata} \
# 			--metadata2 {input.sample_metadata} \
# 			--index {params.index} \
# 			--output {output.combined_metadata}
# 		"""



# ### Aligning the sequences using nextalign
# rule align:
# 	message:
# 		"""
# 		Aligning sequences to {input.reference}
# 		    - gaps relative to reference are considered real
# 		"""
# 	input:
# 		sequences = rules.multifasta.output.combined_seqs,
# 		map = files.genemap,
# 		reference = files.reference
# 	params:
# 		threads = options.threads
# 	output:
# 		alignment = "output_files/sequences/aligned.fasta",
# 		insertions = "output_files/sequences/insertions.csv"
# 	shell:
# 		"""
# 		nextalign \
# 			--sequences {input.sequences} \
# 			--reference {input.reference} \
# 			--genemap {input.map} \
# 			--genes E,M,N,ORF10,ORF14,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
# 			--output-dir output_files/sequences/ \
# 			--output-basename nextalign
		
# 		mv output_files/sequences/nextalign.aligned.fasta {output.alignment}
# 		mv output_files/sequences/nextalign.insertions.csv {output.insertions}
# 		"""



# ### Detect genetic changes and potential sequencing errors
# rule mutations:
# 	message:
# 		"""
# 		Extract mutations from multiple sequence alignment
# 		"""
# 	input:
# 		alignment = rules.align.output.alignment,
# 		annotation = files.annotation,
# 		rename_file = rules.inspect_metadata.output.rename,
# 		qamatrix2 = rules.inspect_metadata.output.matrix,
# 		insertions = rules.align.output.insertions
# 	params:
# 		index = "sample_id",
# 		format = "fasta"
# 	output:
# 		list = "output_files/sequences/listseqs.txt",
# 		corelab_alignment = "output_files/sequences/corelab_aligned.fasta",
# 		mutations = "output_files/sequences/mutations.tsv",
# 		matrix = "output_files/qa/qa_matrix3.tsv",
# 		ren_sequences = "output_files/sequences/renamed_seqs.fasta"
# 	shell:
# 		"""
# 		cut -d$'\t' -f 1 {input.rename_file} > {output.list}
# 		echo Wuhan/Hu-1/2019 >> {output.list}
		
# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.alignment} \
# 			--format {params.format} \
# 			--action keep \
# 			--list {output.list} \
# 			--output {output.corelab_alignment}

# 		python.nextstrain scripts/get_mutations.py \
# 			--input {output.corelab_alignment} \
# 			--annotation {input.annotation} \
# 			--ref-name "Wuhan/Hu-1/2019" \
# 			--ignore ambiguity synonymous \
# 			--output {output.mutations}
			
# 		python.nextstrain scripts/inspect_mutations.py \
# 			--matrix {input.qamatrix2} \
# 			--index {params.index} \
# 			--mutations {output.mutations} \
# 			--insertions {input.insertions} \
# 			--output {output.matrix}
		
# 		python.nextstrain scripts/masterkey.py \
# 			--input {input.alignment} \
# 			--format {params.format} \
# 			--action rename \
# 			--list {input.rename_file} \
# 			--output {output.ren_sequences}	
# 		"""



# rule pangolin:
# 	message:
# 		"""
# 		Assign pango lineages to core lab sequences
# 		"""
# 	input:
# 		corelab_seq = rules.mutations.output.corelab_alignment,
# 		metadata = rules.inspect_metadata.output.metadata,
# 		qamatrix3 = rules.mutations.output.matrix
# 	params:
# 		threads = options.threads,
# 		index = "sample_id"
# 	output:
# 		report = "output_files/metadata/lineage_report.csv",
# 		matrix = "output_files/qa/qa_matrix4.tsv",
# 		metadata = "output_files/assured_data/metadata.tsv"
# 	shell:
# 		"""
# 		pangolin \
# 			{input.corelab_seq} \
# 			--threads {params.threads} \
# 			--outfile {output.report}
		
# 		python.nextstrain scripts/inspect_lineages.py \
# 			--lineages {output.report} \
# 			--metadata {input.metadata} \
# 			--matrix {input.qamatrix3} \
# 			--index {params.index} \
# 			--output1 {output.matrix} \
# 			--output2 {output.metadata}
# 		"""



# ### Masking alignment sites
# rule mask:
# 	message:
# 		"""
# 		Mask bases in alignment
# 		  - masking {params.mask_from_beginning} from beginning
# 		  - masking {params.mask_from_end} from end
# 		  - masking other sites: {params.mask_sites}
# 		"""
# 	input:
# 		alignment = rules.mutations.output.ren_sequences
# 	params:
# 		mask_from_beginning = 55,
# 		mask_from_end = 100,
# 		mask_sites = "150 153 635 1707 1895 2091 2094 2198 2604 3145 3564 3639 3778 4050 5011 5257 5736 5743 5744 6167 6255 6869 8022 8026 8790 8827 8828 9039 10129 10239 11074 11083 11535 13402 13408 13476 13571 14277 15435 15922 16290 16887 19298 19299 19484 19548 20056 20123 20465 21550 21551 21575 22335 22516 22521 22661 22802 24389 24390 24622 24933 25202 25381 26549 27760 27761 27784 28253 28985 29037 29039 29425 29553 29827 29830"
# 	output:
# 		alignment = "output_files/tree/masked_alignment.fasta"
# 	shell:
# 		"""
# 		python.nextstrain scripts/mask-alignment.py \
# 			--alignment {input.alignment} \
# 			--mask-from-beginning {params.mask_from_beginning} \
# 			--mask-from-end {params.mask_from_end} \
# 			--mask-sites {params.mask_sites} \
# 			--output {output.alignment}
# 		"""



# ### Inferring Maximum Likelihood tree using the default software (IQTree)

# rule tree:
# 	message: "Building tree"
# 	input:
# 		alignment = rules.mask.output.alignment
# 	params:
# 		threads = options.threads
# 	output:
# 		tree = "output_files/tree/tree_raw.nwk"
# 	shell:
# 		"""
# 		augur tree \
# 			--alignment {input.alignment} \
# 			--nthreads {params.threads} \
# 			--output {output.tree}
# 		"""



# ### Running TreeTime to estimate time for ancestral genomes

# rule root2tip:
# 	message:
# 		"""
# 		Perform root-to-tip analysis:
# 		  - detect molecular clock outliers
# 		  - generate root-to-tip regression plot
# 		"""
# 	input:
# 		tree = rules.tree.output.tree,
# 		alignment = rules.mask.output.alignment,
# 		metadata = rules.combine_metadata.output.combined_metadata
# 	params:
# 		clock_rate = 0.0008,
# 		clock_std_dev = 0.0004,
# 		root = "Wuhan/Hu-1/2019 Wuhan/WH01/2019",
# 		outdir = "./output_files/root2tip"
# 	output:
# 		outliers = "output_files/root2tip/outliers_list.txt"
# 	shell:
# 		"""
# 		treetime \
# 			--tree {input.tree} \
# 			--aln {input.alignment} \
# 			--dates {input.metadata} \
# 			--clock-filter 3 \
# 			--reroot {params.root} \
# 			--gtr JC69 \
# 			--clock-rate {params.clock_rate} \
# 			--clock-std-dev {params.clock_std_dev} \
# 			--max-iter 2 \
# 			--coalescent skyline \
# 			--plot-rtt root2tip_plot.pdf \
# 			--tip-labels \
# 			--verbose 1 \
# 			--outdir {params.outdir}
		
# 		echo "--" >> {params.outdir}/dates.tsv
# 		grep "\-\-" {params.outdir}/dates.tsv | cut -d$'\t' -f 1 > {output.outliers}
# 		"""



# rule assurance:
# 	message:
# 		"""
# 		Perform last step of quality assurance:
# 		  - generate final QA matrix
# 		  - move sequence and metadata files to final directory
# 		"""
# 	input:
# 		qamatrix4 = rules.pangolin.output.matrix,
# 		outliers = rules.root2tip.output.outliers,
# 		lineage_seqs = "output_files/sequences/rep_genomes.txt",
# 		genomes = rules.multifasta.output.filtered_seqs
# 	params:
# 		index = "sample_id",
# 		format = "fasta"
# 	output:
# 		qamatrix = "output_files/assured_data/qa_matrix.tsv",
# 		quality_seqs = "output_files/sequences/sequences.fasta"
# 	shell:
# 		"""
# 		python.nextstrain scripts/get_QAmatrix.py \
# 			--matrix {input.qamatrix4} \
# 			--genomes {input.genomes} \
# 			--outliers {input.outliers} \
# 			--index {params.index} \
# 			--output1 {output.qamatrix} \
# 			--output2 {output.quality_seqs}
# 		"""


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

