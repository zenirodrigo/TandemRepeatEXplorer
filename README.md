# genome-manipulation-satDNA

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Zeni_repgen.sh - Genome Repeat Analysis Pipeline
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

A Bash/Python pipeline for identifying and visualizing satellite DNA arrays in genomes.

Features:
  Processes multiple genome files and reference sequences

  Automated BLAST analysis with configurable parameters

  Merges adjacent repeats within 2000bp (can be changed)

  Generates three visualization types:

      Chromosome-scale repeat distribution

      Heatmap of repeat prevalence

      Scatter plot of repeat sizes

  Handles common FASTA formats (.fa, .fna, .fasta)

Inputs:

  Genome assemblies

  Reference monomer sequences

  Parameters: sequence count, repeat multiplier

Outputs:

  BED files with merged repeat regions

  Publication-quality plots (PNG/PDF):

  chromosomes_with_annotations.*

  array_frequency_heatmap.png

  array_chromosome_vs_size_scatter.png

Dependencies:

        BLAST+ (makeblastdb, blastn)

        Python 3 with BioPython, pandas, matplotlib, seaborn
Usage:        
      ./Zeni_repgen.sh
      # Follow interactive prompts to input files and parameters




-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-
# Zeni_annot_gffbed.sh -  BLAST Analysis and GFF Comparison in genomes assembled
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
bash script.sh
```
And provide the required inputs as prompted.

---

## Contact
For questions or improvements, contact [rodrigo-zeni@outlook.com.br].


