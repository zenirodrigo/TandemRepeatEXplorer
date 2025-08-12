[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)
<img align="right" src="social-preview.png" width="150" alt="Repository Icon">

# T-REx – Tandem Repeat Explorer 🦖🧬

**T-REx** is a toolkit containing four scripts for detecting, annotating, and analyzing tandem repeats from genome assemblies and short-read sequencing data.

- **arrayScope.sh** – Characterizes and locates tandem repeat arrays in assembled genomes.
- **satFlank.sh** – Studies the neighborhood of arrays using the assembled genome and its annotation.
- **monoMiner.py** – Automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.
- **gene_extractor.sh** – Automates the retrieval of a specific gene sequence across multiple species using NCBI tools.

---

## arrayScope.sh – Genome Repeat Analysis Pipeline

### 📜 Description
**ArrayScope** is a **Bash + Python** pipeline designed to identify and visualize satellite DNA arrays in genome assemblies. It processes multiple genomes and reference monomers, runs **BLAST**, merges nearby repeat regions (< 2000 bp apart), and generates high-quality plots.

> 💡 For variants in reference files (e.g., same repeat in different species), use `arrayScope_variants.sh` and include an underscore `_` in the filename.

### ✅ Features
- Processes multiple genomes and reference monomers
- Automated **BLAST** analysis with configurable parameters
- Merges nearby repeats within 2000 bp
- Generates:
  1. Chromosome-scale repeat map
  2. Heatmap of repeat abundance
  3. Scatter plot of array sizes
- Supports `.fa`, `.fna`, `.fasta`
- Parallel processing with **GNU parallel**

### 🛠️ Dependencies
- **BLAST+** (`makeblastdb`, `blastn`)
- **GNU parallel**
- **Python 3** with:
  - BioPython
  - pandas
  - matplotlib
  - seaborn
- **bedtools**
- **awk**

### 📂 Input Files
1. Genome assemblies (FASTA)
2. Reference monomers (FASTA)

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

## satFlank.sh – Satellite Flank Extraction and Gene Overlap Analysis

### 📜 Description
**satFlank.sh** identifies satellite DNA flanking regions with BLAST, extracts sequences, and cross-references them with genome annotations in GFF format.

### ✅ Features
- Multiplies reference sequence to improve BLAST sensitivity
- BLAST search for flanking regions
- Filters `.bed` to keep largest non-overlapping regions
- Sequence extraction with **seqtk**
- Gene overlap analysis with GFF annotations
- Parallel processing

### 🛠️ Dependencies
- **BLAST+**
- **seqtk**
- **GNU parallel**
- **awk**, **sed**

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

---

## monoMiner.py – Automated Motif Mining Pipeline

### 📜 Description
**monoMiner.py** scans sequencing libraries for motifs similar to a reference, saves results to FASTA, filters with **CD-HIT**, and annotates sequences by species.

### ✅ Features
- Reads mapping file (species code → name)
- Scans `.fq` files for motifs
- Saves motifs per library to FASTA
- Concatenates results to `final.fasta`
- Filters with **cd_hit_filter_size.py**
- Outputs `output.tsv` linking sequences to species
- Parallel processing

### 🛠️ Dependencies
- **Python 3**
- **cd_hit_filter_size.py**
- Standard Python libraries (`os`, `glob`, `re`, `subprocess`, `concurrent.futures`)

### 📂 Input Files
1. Mapping file (TSV)
2. Reference sequence (FASTA)
3. `.fq` sequencing files

### 📂 Output Files
- `{reference_name}_motifs_{library_name}.fasta`
- `final.fasta`
- CD-HIT filtered FASTA
- `output.tsv`

### 🚀 Usage
```bash
python monoMiner.py
```

---

## gene_extractor.sh – Automated Gene Retrieval from NCBI

### 📜 Description
**gene_extractor.sh** downloads genomic sequences for a given gene across multiple species from NCBI, with fallback to RNA data if needed.

### ✅ Features
- Search NCBI Gene by symbol and species
- Retrieve genomic coordinates and exon count
- Download sequence with **efetch**
- Reverse-complement negative strands
- RNA fallback using NCBI Datasets CLI

### 🛠️ Dependencies
- **Entrez Direct** (`esearch`, `efetch`, `esummary`, `xtract`)
- **seqkit**
- **datasets** CLI
- **unzip**

### 📂 Input Files
1. Gene symbol
2. Species list file (one species per line)

### 📂 Output Files
- `{GENE_SYMBOL}_{SPECIES_NAME}.fasta`

### 🚀 Usage
```bash
bash gene_extractor.sh BRCA1 species.txt
```

---

## 📜 License
This project is licensed under the MIT License.

## 🙌 Acknowledgments
- Heng Li, Ole Tange, NCBI team, and BEDTools team for essential tools.
- Universidade Estadual Paulista (UNESP).
- Developed during PhD with FAPESP fellowship.

## 📬 Contact
For questions or suggestions: **rodrigo-zeni@outlook.com.br**
