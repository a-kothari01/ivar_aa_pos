## Python script to add the amino acid position to ivar .tsv output file

iVar is a comprehensive tool that plays an important role in variant calling, especially for viral amplicon-based sequencing data. It can identify codons and translate variants into amino acids when supplied with a GFF file for the reference genome (in GFF3 format). The output is generated as a .tsv file that contains information like the position of a base on the reference sequence (POS), the reference base (REF), and the alternate base (ALT) among several others. 

This script, written in Python, has been optimised for use in a HIV-1 genome variant analysis pipeline. It intends to use the POS, REF and ALT columns from the iVar .tsv file to add 2 more columns to this output, for further downstream analyses of these generated variants. 
These columns are –   

1. **POS_AA** (the position of the amino acid where the mutation has occurred) and
2. **MUT_AA** (variant information in the format - {reference_aa}{position_aa}{alternate_aa}).

The latter is for ease of variant list extraction for input into the Stanford HIVDB for further analyses. 

### Workflow Overview

The script contains of 3 main segments –

1. **GFF file parsing function –** it parses the reference GFF to extract the desired CDS features of those regions, including the CDS start and end coordinates, and name of the corresponding gene.
2. **Amino acid annotation function –** It receives individual rows from the input .tsv files as a pandas (python library) Series
3. **TSV file processor –** iterates over every row of all input tsv files, concatenates them with the values received from the annotation function and generates a new tsv file as a pandas dataframe with the POS_AA and MUT_AA columns

### Directory Structure

.\
├── hiv1cref.fasta  *#Reference genome (FASTA)* \
├── hiv1cref.gff3 *#Reference annotation (GFF3)*\
├── ivar_pos_aa.py *#workflow definition*\
├── *variants.tsv *#input .tsv files for each sample*

### Dependencies

The entire pipeline has been executed on Ubuntu (v22.04) terminal.

Make sure the following tools are installed and available in your PATH:
- Python ≥ 3.8 

Python libraries/packages that need to be installed -
- Biopython 
- pandas

### Usage

1.	Place the input .tsv files (e.g., sample1_variants.tsv) in the working directory.
2.	Place the reference genome (.fasta) and annotation (.gff3) files in the working directory.
3.	Run the script with python:
   
bash\
python3 ivar_aa_pos.py

4.	Pipeline automatically finds all samples based on *.tsv patterns
5.	Outputs for each sample input will be generated in the respective directory (e.g., Sample1_variants_annotated.tsv).

### Configuration 

**Reference:** The program expects the reference files hiv1cref.fasta and hiv1cref.gff3 in the root directory.

**Input:** iVar output .tsv files named as: \
SAMPLE1_variants.tsv

The script may be edited directly for any adjustments in file paths or parameters.

### GFF3 Format

The code works on the assumption that the GFF file supplied follows the basic 9-column GFF3 format as described at the below link – \
[gff3 specifications](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

### .tsv File Format

The code works on the assumption that the columns in the .tsv file generated have the POS, REF, and ALT columns as described at the below link (Call variants with iVar section) – \
[iVar Manual](https://andersen-lab.github.io/ivar/html/manualpage.html)

### Citation

- [iVar](https://rdcu.be/eDoXw) 
- Python
- [biopython](http://dx.doi.org/10.1093/bioinformatics/btp163) 
- [pandas](https://doi.org/10.5281/zenodo.3509134)
- [GFF Basic Format](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/#basicformat)

### Notes

1. Can be optimised for specific reference genome as desired.
2. The gene names selected for cds_features can be updated and filtered as necessary.
3. Can be optimised for bases read from the negative strand as well.
4. File path and nomenclature for both input and output files can be updated as required.
5. Input files must have some pattern in nomenclature for the glob module to detect all files automatically.
6. GFF file must be in the typical GFF3 format for the code to work.

