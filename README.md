# dengueQA - a pipeline for quality assurance of consensus dengue virus genomes

This pipeline processes dengue virus sequencing data, performing genotyping, sequence alignment, genome coverage calculation, and molecular clock analysis.

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/InstitutoTodosPelaSaude/dengueQA.git
    cd dengueQA
    ```

2. Install dependencies:
    ```bash
    conda env create -f config/dengueQA.yaml
    conda activate dengueQA
    ```

## Usage

1. Edit the Snakefile to specify your input files and parameters.
2. Run the pipeline:
    ```bash
    snakemake all --cores all
    ```

## Workflow

The pipeline consists of the following steps:

1. **BLAST Searches**: Perform BLAST searches for serotyping, genotyping, and lineage assignment according to 'dengue-lineages' nomenclature system.
2. **Dataset splitting**: Split datasets based on serotypes.
3. **Sequence alignment**: Align newly sequenced genomes with reference genomes.
4. **Coverage calculation**: Calculate coverage per coding region (CDS) and whole genomes.
5. **Contextual genome selection**: Select contextual genomes for root-to-tip analysis.
6. **Dataset compilation**: Join new sequences with contextual genomes.
7. **Metadata processing**: Add missing columns and process metadata.
8. **Phylogenetic analysis**: Align genomes, mask regions, build trees, and perform root-to-tip analysis.
9. **Quality assurance**: Generate quality assurance report containing information about genomic coverage (per coding region (CDS) and whole genome) and sequence quality (based on molecular clock analysis).


## Author

* **Anderson Brito, Instituto Todos pela Saúde (ITpS)** - [Website](https://www.itps.org.br/membros) - anderson.brito@itps.org.br

## License

This project is licensed under the MIT License.
