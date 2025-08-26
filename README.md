[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)

<img align="right" src="social-preview.png" width="150" alt="Repository Icon">
<br clear="right">
<h1>Tandem Repeat Explorer 🦖🧬</h1>


This repository, also called **T-REx**, contains four scripts for analyzing tandem sequences from genomes (first two scripts) and short reads (**monoMiner**):

- **arrayScope.sh** – Characterizing and locating tandem repeat arrays in assembled genomes.
- **satFlank.sh** – Studying the neighborhood of arrays using the assembled genome and its annotation.
- **monoMiner.py** – Automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.
- **gene_extractor.sh** – Automates the retrieval and extraction of the genomic sequence of a specific gene across multiple species using NCBI tools.

---

## arrayScope.sh – Genome Repeat Analysis Pipeline

### 📜 Description
**ArrayScope** is a **Bash + Python** pipeline designed to identify and visualize satellite DNA arrays in genome assemblies.  
It processes multiple genomes and reference monomers, runs **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**, merges nearby repeat regions (< 2000 bp apart), and generates high-quality plots for chromosome-scale visualization.

> 💡 For variants in reference files — e.g., two slightly different monomers or the same repeat from different species — use the `arrayScope_variants.sh` script and include an underscore `_` in the filename.

### ✅ Features
- Processes **multiple genome assemblies** and reference monomer files
- Fully automated **BLAST** analysis with configurable parameters
- Merges nearby repeat hits (< 2000 bp apart)
- Generates three publication-quality visualizations:
  1. Chromosome-scale repeat distribution (linear map)
  2. Heatmap of total base pairs per chromosome/reference
  3. Scatter plot of array sizes by chromosome
- Handles `.fa`, `.fna`, `.fasta` formats
- Parallel processing using **[GNU parallel](https://www.gnu.org/software/parallel/)**

### 🛠️ Dependencies
- **BLAST+** (`makeblastdb`, `blastn`)
- **[GNU parallel](https://www.gnu.org/software/parallel/)**
- **Python 3** with:
  - [BioPython](https://biopython.org/)
  - [pandas](https://pandas.pydata.org/)
  - [matplotlib](https://matplotlib.org/)
  - [seaborn](https://seaborn.pydata.org/)
- **[bedtools](https://bedtools.readthedocs.io/)**
- **awk**

**Quick install example:**
```bash
sudo apt-get install bedtools
conda install -c bioconda blast
conda install -c conda-forge parallel
pip install biopython pandas matplotlib seaborn
```

### 📂 Input Files
1. Genome assemblies (FASTA format)
2. Reference monomers (FASTA format)

### 📂 Output Files
- `valid_monomers.bed`
- `chromosomes_with_annotations.png`
- `array_frequency_heatmap.png`
- `array_chromosome_vs_size_scatter.png`

### 🚀 Usage
```bash
bash arrayScope.sh
```

---

## gene_extractor.sh – Extract Gene Sequences from NCBI by Species

### 📜 Description
Automates the retrieval of genomic sequences for a given gene across multiple species using **[NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)** and **datasets CLI**.

### ✅ Features
- Searches for Gene ID using gene symbol and species name
- Extracts genomic coordinates from NCBI
- Fallback to RNA FASTA if coordinates are missing
- Downloads sequences with `efetch`
- Reverse-complements negative strands automatically
- Outputs FASTA files named by gene and species

### 🛠️ Dependencies
```bash
conda install -c bioconda entrez-direct datasets-cli seqkit
```
Also requires: unzip, cut, sed, find, head

### 📂 Input Files
- Gene symbol (e.g., BRCA1)
- Text file with one species per line

### 📂 Output Files
- `{GENE_SYMBOL}_{SPECIES_NAME}.fasta`

### 🚀 Usage
```bash
bash gene_extractor.sh BRCA1 species.txt
```

---

## satFlank.sh – SatelliteDNA Flank Extraction and Gene Overlap Analysis Pipeline

### 📜 Description
**satFlank.sh** identifies satellite DNA flanking regions using **BLAST**, extracts sequences, and compares with GFF annotations to find overlapping genes.

### ✅ Features
- Multiplies reference sequence for improved BLAST sensitivity
- Filters `.bed` to keep largest non-overlapping regions
- Sequence extraction with **seqtk**
- Gene overlap analysis from GFF
- Parallel processing

### 🛠️ Dependencies
```bash
sudo apt-get install blast2 seqtk parallel
```

### 📂 Input Files
1. Library file (FASTA)
2. Reference sequence (FASTA)
3. GFF annotation file

### 📂 Output Files
- `monomeros_validos.filtrado.bed`
- Extracted sequences (FASTA)
- Gene overlap table (TSV)

### 🚀 Usage
```bash
bash satFlank.sh
```
You will be prompted to enter:
1. Library file name (FASTA)
2. Reference file name (FASTA)
3. Number of times to multiply the reference sequence
4. BLAST output file name
5. Output file name for extracted sequences
6. Number of CPU cores for parallel processing
7. GFF file name
8. Final output file name for GFF analysis


## 📌 Notes & Best Practices
- Use biologically relevant reference sequences for better BLAST detection
- Ensure the `.gff` file matches the genome assembly used in the library file


---

## monoMiner.py – Genome Repeat Analysis Pipeline

### 📜 Description
Automated Python pipeline for motif mining from sequencing libraries.

### 🛠️ Dependencies
- **Python 3**
- **cd_hit_filter_size.py**
- Standard Python libraries
- FASTQ sequencing files
- Mapping file (TSV: species_code → species_name)
- **Install CD-HIT (required)**:

Follow instructions at [CD-HIT](https://github.com/weizhongli/cdhit) Official Repository


Ensure cd_hit_filter_size.py is in your PATH.


### 🚀 Usage
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
> Minimum copies for CD-HIT: 2
```

##  Example Mapping File (species_mapping.tsv)
```bash
ame    Astyanax mexicanus
dre   Danio rerio
cpo    Catopsilia pomona

```

```bash
python3 monoMiner.py
```

## 📂 Output Files
.

├── final.fasta                    # Concatenated motifs

├── output.tsv                     # Annotated sequences

├── reference_motifs_*.fasta       # Per-library motif files

└── final.fasta.nr0.*.sel.fasta    # CD-HIT filtered results



---

## 📜 License

This project is licensed under the MIT License

## Acknowledgments
- Heng Li, Ole Tange, NCBI team, and BEDTools team for developing essential bioinformatics tools.
- Support from Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).
- This code was developed in my PhD with fellowship from Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP)

## ⁉️ Contact
For questions or improvements, contact:
**rodrigo.zeni@unesp.br**
or
**rodrigo-zeni@outlook.com.br**
