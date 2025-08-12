[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)
<img align="right" src="social-preview.png" width="150" alt="Ícone do Repositório">

## This repository is also called T-REx and contains four scripts for analyzing tandem sequences using genome (the first two scripts) and short reads (monoMiner)🦖🧬

   - arrayScope.sh: A script for characterizing and locating tandem repeat arrays in assembled genomes.
   - satFlank.sh: A script for studying the neighborhood of arrays using the assembled genome and its annotation.
   - monoMiner.py: An automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.
   - gene_extractor.sh This script automates the retrieval and extraction of the genomic sequence of a specific gene across multiple species using NCBI tools 

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## arrayScope.sh - Genome Repeat Analysis Pipeline
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Description
ArrayScope is a Bash + Python pipeline designed to identify and visualize satellite DNA arrays in genome assemblies.
It processes multiple genomes and reference monomers, runs BLAST, merges nearby repeat regions (< 2000 bp apart), and generates high-quality plots for chromosome-scale visualization.

💡 If you have variants in your reference files — for example, two slightly different monomers or the same repeat from different species — use the arrayScope_variants.sh script and include an underscore _ in the filename.
Example:
ppfia1_Astyanax_mexicanus.fasta and ppfia1_Psalidodon_paranae.fasta will both be treated as ppfia1.

## ✅  Features
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
conda install -c conda-forge parallel

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


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## gene_extractor.sh – Extract Gene Sequences from NCBI by Species
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Description
This script automates the retrieval of genomic sequences for a given gene across multiple species using NCBI's E-utilities and datasets CLI. It fetches genomic coordinates when available, or falls back to transcriptomic data when necessary.

## ✅ Features
- Searches for Gene ID using gene symbol and species name

- Extracts genomic coordinates directly from NCBI

- Fallback to RNA FASTA when coordinates are missing

- Downloads sequences with efetch from GenBank

- Auto-detects and reverse-complements negative strands

- Outputs FASTA files named per gene and species

📊 Logs useful information: coordinates, exon count, sequence size


---


## 🛠️ Dependencies
Install the required tools via Conda:
```bash
conda install -c bioconda entrez-direct datasets-cli seqkit
```
Ensure the following are available in your environment:


Standard UNIX tools: unzip, cut, sed, find, head


---

## 📂 Inputs files:
A gene symbol (e.g., BRCA1)

A text file with one species per line:

Homo sapiens

Mus musculus

Danio rerio


## ▶️ Usage
```bash
bash gene_extractor.sh BRCA1 especies.txt
```
This will generate FASTA files like:

BRCA1_Homo_sapiens.fasta
BRCA1_Mus_musculus.fasta
BRCA1_Danio_rerio.fasta

🧪 Fallback Mode
If genomic coordinates cannot be extracted, the script automatically switches to RNA download mode using NCBI Datasets. It unpacks the archive and extracts the first available .fna file. This ensures recovery even for less-annotated organisms.


---

## 📂 Outputs files:
For each species, you’ll see messages like:
```bash
🔍 Searching for gene BRCA1 in Homo sapiens...
🧠 Gene ID found: 672
🗺️ Coordinates: NC_000017.11:43044294-43125483 (Exons: 24)
🌐 Downloading sequence...
✅ Saved as BRCA1_Homo_sapiens.fasta
🔄 Strand: reverse
```

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## satFlank.sh -  BLAST Analysis and GFF Comparison in genomes assembled
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

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## monoMiner.py - Genome Repeat Analysis Pipeline
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
- This code was developed in my PhD with fellowship from Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP)

## ⁉️ Contact
For questions or improvements, contact rodrigo-zeni@outlook.com.br.


