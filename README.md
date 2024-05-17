# dengueQA - a pipeline for quality assurance of consensus dengue virus genomes

This pipeline processes dengue virus sequencing data to perform quality assurance of genomic data. To do so it runs steps of genotyping, sequence alignment, genome coverage calculation, and molecular clock analysis.

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

1. **BLAST searches**: Perform BLAST searches for serotyping, genotyping, and lineage assignment according to reference files provided by ['dengue-lineages'](https://dengue-lineages.org/design.html) nomenclature system.
2. **Dataset splitting**: Split datasets based on dengue serotypes.
3. **Sequence alignment**: Align newly sequenced genomes with reference genomes using `augur align`.
4. **Coverage calculation**: Calculate coverage per coding region (CDS) and whole genomes.
5. **Contextual genome selection**: Select contextual genomes for root-to-tip analysis.
6. **Dataset compilation**: Join new sequences with contextual genomes.
7. **Metadata processing**: Add missing columns and process metadata.
8. **Phylogenetic analysis**: Align genomes (`augur align`), mask 5' and 3' untranslated regions (`augur mask`), build phylogenetic trees (`iqtree`), and perform root-to-tip analysis (`treetime`).
9. **Quality assurance**: Generate quality assurance report containing information about genomic coverage (flagging low coverage coding regions (CDS) and whole genomes) and sequence quality (flagging molecular clock outliers).

## Author

* **Anderson Brito, Instituto Todos pela Saúde (ITpS)** - [Website](https://www.itps.org.br/membros) - anderson.brito@itps.org.br

## License

This project is licensed under the MIT License.
