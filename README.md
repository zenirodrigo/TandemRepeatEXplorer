[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)

<img align="right" src="social-preview.png?raw=1" width="250" alt="Repository Icon">
<br clear="right">

# Tandem Repeat Explorer 🦖🧬

Tandem Repeat Explorer (T-REx) is a modular **Bash/Python toolkit** for the identification, characterization, and visualization of tandem arrays from genome assemblies and short-read sequencing data.

## Toolkit modules
- **arrayScope.sh** – Characterizing and locating tandem repeat arrays in assembled genomes.
- **satDNA_density.sh** – Circos-like visualization of satDNA density and GC content along chromosomes.
- **satFlank.sh** – Studying the neighborhood of arrays using the assembled genome and its annotation.
- **monoMiner.py** – Automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.
- **gene_extractor.sh** – Automates the retrieval and extraction of the genomic sequence of a specific gene across multiple species using NCBI tools.

---

## 📚 Table of Contents
- [Quick Start](#quick-start)
- [arrayScope.sh](#arrayscopesh--genome-repeat-analysis-pipeline)
- [satDNA_density.sh](#satdna_densitysh--circos-like-satelitome-density-visualization)
- [satFlank.sh](#satflanksh--satellitedna-flank-extraction-and-gene-overlap-analysis-pipeline)
- [monoMiner.py](#monominerpy--reference-guided-motif-mining-from-short-reads)
- [gene_extractor.sh](#gene_extractorsh--extract-gene-sequences-from-ncbi-by-species)
- [Notes & Best Practices](#-notes--best-practices)
- [License](#-license)
- [Acknowledgments](#acknowledgments)
- [Contact](#-contact)
- [Citation](#-citation)

---

## 🚀 Quick Start

```bash
git clone https://github.com/zenirodrigo/TandemRepeatEXplorer.git
cd T-REx

bash arrayScope.sh
bash satDNA_density.sh
bash satFlank.sh
python3 monoMiner.py
bash gene_extractor.sh BRCA1 species.txt
```

---

## arrayScope.sh – Genome Repeat Analysis Pipeline

### Description
arrayScope.sh is a **Bash + Python** pipeline designed to identify and visualize satellite DNA arrays in genome assemblies.  
It processes multiple genomes and reference monomers, runs **BLAST+**, merges nearby repeat regions (< 2000 bp apart), and generates high-quality plots for chromosome-scale visualization.

For variants in reference files — e.g., two slightly different monomers or the same repeat from different species — include an underscore `_` in the filename.

### Features
- Processes multiple genome assemblies and reference monomer files
- Fully automated BLAST analysis with configurable parameters
- Merges nearby repeat hits (< 2000 bp apart)
- Generates three publication-quality visualizations:
  1. Chromosome-scale repeat distribution (linear map)
  2. Heatmap of total base pairs per chromosome/reference
  3. Scatter plot of array sizes by chromosome
- Handles `.fa`, `.fna`, `.fasta` formats
- Parallel processing using GNU parallel

### Dependencies
- BLAST+ (`makeblastdb`, `blastn`)
- GNU parallel
- Python 3 with:
  - biopython
  - pandas
  - matplotlib
  - seaborn
- bedtools
- awk

### Usage
```bash
bash arrayScope.sh
```

---

## satDNA_density.sh – Circos-like satelitome density visualization

### Description
satDNA_density.sh generates a **circos-like overview** of satellite DNA distribution across chromosomes or contigs.

The pipeline runs BLASTn searches of satDNA reference sequences against a genome assembly, merges nearby hits (≤ 2000 bp), and produces a circular plot including:

- Outer ring: chromosomes/contigs (preserving the original FASTA order)
- Ruler: tick marks every 10 Mb and labels every 50 Mb
- Main ring: stacked satDNA coverage per 100 kb bin
- Inner track: GC fraction per 100 kb bin
- Top-10 highlighting: either the first 10 satDNAs in FASTA order or the 10 most abundant satDNAs (by total bp covered)

Chromosome identifiers are never re-ordered or re-indexed.  
Only a safe label abbreviation is applied: `ChromosomeN` → `Chr N`.

### Features
- Accepts:
  - One multi-FASTA satelitome (automatically split into one FASTA per satDNA), or
  - Multiple single-FASTA reference files
- Efficient parallel BLAST execution using GNU parallel
- Two Top-10 coloring strategies:
  1. First 10 satDNAs in reference FASTA order  
  2. Top 10 most abundant satDNAs (by bp covered)
- Automatic GC content profiling
- Publication-ready figures and tables

### Dependencies
- BLAST+
- GNU parallel
- Python 3.8+ with numpy, pandas, matplotlib

### Usage
```bash
bash satDNA_density.sh
```

---

## satFlank.sh – SatelliteDNA Flank Extraction and Gene Overlap Analysis Pipeline

### Description
satFlank.sh identifies satellite DNA flanking regions using BLAST, extracts sequences, and compares them with GFF annotations to detect overlapping genes.

### Usage
```bash
bash satFlank.sh
```

---

## monoMiner.py – Reference-guided Motif Mining from Short Reads

### Description
monoMiner.py performs reference-guided motif discovery from sequencing libraries with clustering and filtering.

### Usage
```bash
python3 monoMiner.py
```

---

## gene_extractor.sh – Extract Gene Sequences from NCBI by Species

### Description
Automates the retrieval of genomic sequences for a given gene across multiple species.

### Usage
```bash
bash gene_extractor.sh BRCA1 species.txt
```

---

## 📌 Notes & Best Practices
- Use biologically meaningful reference sequences.
- Ensure genome FASTA and annotations are consistent.
- Prefer high-core machines for large satelitomes.

---

## 📜 License
This project is licensed under the MIT License.

---

## Acknowledgments
- Ole Tange (GNU Parallel), NCBI, BEDTools community.
- Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).
- Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP).

---

## ⁉️ Contact
rodrigo.zeni@unesp.br  
rodrigo-zeni@outlook.com.br  

---

## 📖 Citation
If you use T-REx in your research, please cite the associated publication (in preparation).
