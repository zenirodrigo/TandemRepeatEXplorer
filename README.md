[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)


# This repository is called T-REx and contains three scripts  🦖🧬:
-   ArrayScope.sh: A script for characterizing and locating tandem repeat arrays in assembled genomes.
-   SatFlank.sh: A script for studying the neighborhood of arrays using the assembled genome and its annotation.
-   MonoMiner.py: An automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ArrayScope.sh - Genome Repeat Analysis Pipeline
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Description
A Bash and Python pipeline for identifying and visualizing satellite DNA arrays in genomes.

---


## Features
- Processes multiple genome files and reference sequences
- Automated BLAST analysis with configurable parameters
- Merges adjacent repeats within 2000bp
- Generates three visualization types:
- Chromosome-scale repeat distribution
- Heatmap of repeat prevalence
- Scatter plot of repeat sizes
- Handles common FASTA formats (.fa, .fna, .fasta)

---

## 🛠️ Dependencies
- BLAST+ (makeblastdb, blastn)
- Python 3 with BioPython, pandas, matplotlib, seaborn

---

## 📥 Installation
Before running the script, ensure that the required dependencies are installed:

```bash
git clone https://github.com/zenirodrigo/TandemRepeatEXplorer.git
cd TandemRepeatEXplorer

sudo apt-get install bedtools
pip install biopython pandas matplotlib seaborn
```

---

## 📂 Inputs files:
- Genome assemblies
- Reference monomer sequences
- Parameters: sequence count, repeat multiplier

---

## 📂 Outputs files:
- BED files with merged repeat regions
- Publication-quality plots (PNG/PDF):


- chromosomes_with_annotations.*png and pdf

- array_frequency_heatmap.png

- array_chromosome_vs_size_scatter.png

---

## Usage
Run the script by executing:

```bash
bash ArrayScope.sh
```

During execution, the script will prompt the user to input the threshold for filtering sequences. This value will be used dynamically throughout the script.

## Example
Upon execution, the script will ask for a filtering threshold:
```bash
Enter the sequence filtering threshold: 10
```
This threshold will be applied throughout the analysis, influencing clustering and output files.




-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## SatFlank.sh -  BLAST Analysis and GFF Comparison in genomes assembled
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Description
This script automates the process of:
1. **BLAST Analysis**: Searches sequences in a library using BLAST and generates a `.bed` file with valid monomer intervals.
2. **Filtering the `.bed` File**: Retains only the largest overlapping regions to avoid redundancy.
3. **Sequence Extraction**: Uses `seqtk` to extract filtered sequences.
4. **Comparison with GFF**: Analyzes gene annotations for the identified regions, generating a final report.

---

## 🛠️ Dependencies
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

##📂 Outputs
- `valid_monomers.filtered.bed`: `.bed` file with filtered intervals.
- `<output_seqtk>`: FASTA file with extracted sequences.
- `<final_file>`: Tabular file containing annotated genes in the identified regions.

---

## Usage
Run the script and follow the interactive instructions:
```bash
bash SatFlank.sh
```
The user will need to provide:
1. **genome file** (FASTA)
2. **Reference file** (FASTA)
3. **Reference sequence multiplier**
4. **BLAST output file**
5. **Output file for sequence extraction**
6. **Number of cores for parallel processing**
7. **`.gff` file for comparison**
8. **Final analysis GFF output file name**

---

## Execution Example
```bash
bash SatFlank.sh
```
And provide the required inputs as prompted.

---

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MonoMiner.py - Genome Repeat Analysis Pipeline
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 📥 Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/zenirodrigo/TandemRepeatEXplorer.git
   cd TandemRepeatEXplorer
   
2. **Install CD-HIT (required)**:

Follow instructions at [CD-HIT](https://github.com/weizhongli/cdhit) Official Repository


Ensure cd_hit_filter_size.py is in your PATH.


Usage:

Step 1: Prepare Files
Reference: A .fasta file with the reference sequence.

Libraries: Place .fq files in the project folder.

Mapping: A .tsv file linking codes to species (example below).

Step 2: Run the Pipeline
```bash
python3 monoMiner.py
```
Follow the prompts:
```bash
> Enter mapping file path: species_mapping.tsv
> Enter reference sequence: examples/reference.fasta
> Minimum copies for CD-HIT: 5
```

## 🧩 Example Mapping File (species_mapping.tsv)
```bash
ame    Astyanax mexicanus
dre   Danio rerio
cpo    Catopsilia pomona

```

## 📂 Output Structure
.

├── final.fasta                    # Concatenated motifs

├── output.tsv                     # Annotated sequences

├── reference_motifs_*.fasta       # Per-library motif files

└── final.fasta.nr0.*.sel.fasta    # CD-HIT filtered results

## 🛠️ Key Features

- TAB autocompletion for file paths
- Parallel processing
- Adjustable similarity threshold (similarity_threshold)
- CD-HIT integration for duplicate removal
- TSV output for downstream analysis

## 📜 License

This project is licensed under the MIT License

## Acknowledgments
- Heng Li, Ole Tange, NCBI team, and BEDTools team for developing essential bioinformatics tools.
- Support from Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).
![Repository Icon](social-preview.png)
- This code was developed in my PhD with fellowship from Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP)

## ⁉️ Contact
For questions or improvements, contact rodrigo-zeni@outlook.com.br.


