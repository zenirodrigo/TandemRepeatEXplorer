# genome-manipulation-satDNA

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Zeni_repgen.sh - Genome Repeat Analysis Pipeline
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Description
This script processes genomic sequences using CD-HIT, allowing the user to define a variable threshold for filtering sequences. The script extracts relevant annotations, generates BED files, and creates FASTA files with extended sequences from the reference genome.

## Features
- Filters sequences using a user-defined threshold instead of a fixed value.
- Generates BED files with reference annotations.
- Extracts sequences from the reference genome, including 5000 bp before and after the annotated region.
- Ensures consistent processing across multiple libraries.

## Requirements
- **CD-HIT** (Cluster Database at High Identity with Tolerance)
- **Python 3.x**
- **BioPython** (for FASTA processing)
- **BEDTools** (for genomic coordinate manipulations)

## Installation
Before running the script, ensure that the required dependencies are installed:

```bash
sudo apt-get install bedtools
pip install biopython
```

Ensure that CD-HIT is installed and accessible in your system:
```bash
sudo apt-get install cd-hit
```

## Usage
Run the script by executing:

```bash
bash script.sh
```

During execution, the script will prompt the user to input the threshold for filtering sequences. This value will be used dynamically throughout the script.

## Input Files
- **FASTA files**: Contain genomic sequences to be processed.
- **BED files**: Store annotations with genomic coordinates.

## Output Files
- **Filtered FASTA files**: Contain sequences after applying the user-defined threshold.
- **BED annotation files**: Provide processed reference annotations.
- **Extracted FASTA files**: Contain reference genome sequences with 5000 bp padding on both ends.

## Example
Upon execution, the script will ask for a filtering threshold:
```bash
Enter the sequence filtering threshold: 10
```
This threshold will be applied throughout the analysis, influencing clustering and output files.




-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Zeni_annot_gffbed.sh -  BLAST Analysis and GFF Comparison in genomes assembled
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Description
This script automates the process of:
1. **BLAST Analysis**: Searches sequences in a library using BLAST and generates a `.bed` file with valid monomer intervals.
2. **Filtering the `.bed` File**: Retains only the largest overlapping regions to avoid redundancy.
3. **Sequence Extraction**: Uses `seqtk` to extract filtered sequences.
4. **Comparison with GFF**: Analyzes gene annotations for the identified regions, generating a final report.

---

## Dependencies
Before running the script, make sure you have installed:
- `BLAST+`
- `seqtk`
- `awk`
- `parallel`

To install these tools on Debian/Ubuntu-based systems:
```bash
sudo apt update && sudo apt install ncbi-blast+ seqtk parallel
```

---

## Usage
Run the script and follow the interactive instructions:
```bash
bash script.sh
```
The user will need to provide:
1. **Library file** (FASTA)
2. **Reference file** (FASTA)
3. **Reference sequence multiplier**
4. **BLAST output file**
5. **Output file for sequence extraction**
6. **Number of cores for parallel processing**
7. **`.gff` file for comparison**
8. **Final analysis GFF output file name**

---

## Outputs
- `valid_monomers.filtered.bed`: `.bed` file with filtered intervals.
- `<output_seqtk>`: FASTA file with extracted sequences.
- `<final_file>`: Tabular file containing annotated genes in the identified regions.

---

## Execution Example
```bash
bash Zeni_annot_gffbed.sh
```
And provide the required inputs as prompted.

---
## Acknowledgments
- Thanks to the CD-HIT, NCBI and BEDTools teams for developing essential bioinformatics tools.
- Support from Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).

## Contact
For questions or improvements, contact [rodrigo-zeni@outlook.com.br].


