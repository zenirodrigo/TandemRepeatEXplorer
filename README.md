[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT) 
[![BLAST](https://img.shields.io/badge/BLAST-Required-0077CC?logo=NCBI&logoColor=white)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Shell Script](https://img.shields.io/badge/Shell_Script-Compatible-brightgreen?logo=gnu-bash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python 3](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python&logoColor=white)](https://www.python.org/)

<img align="right" src="social-preview.png?raw=1" width="250" alt="Repository Icon">
<br clear="right">
<h1>Tandem Repeat Explorer 🦖🧬</h1>


Tandem Repeat Explorer (T-REx) is a modular Bash/Python toolkit for the identification, characterization, and visualization of tandem arrays from genome assemblies and short-read sequencing data.

- **arrayScope.sh** – Characterizing and locating tandem repeat arrays in assembled genomes.
- **satDNA_density.sh** - Generating a circus overview of satDNAs distributions across a  genome assembly.
- **satDNA_similarity.py** – Biological clustering of satellite DNA monomers into families, accounting for circularity, reverse-complement, and indels.
- **satFlank.sh** – Studying the neighborhood of arrays using the assembled genome and its annotation.
- **monoMiner.py** – Automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and intelligent filtering.
- **gene_extractor.sh** – Automates the retrieval and extraction of the genomic sequence of a specific gene across multiple species using NCBI tools.

---
## 📚 Table of Contents
- [arrayScope.sh](#arrayscopesh--genome-repeat-analysis-pipeline)
- [satDNA_density.sh](#satdna_densitysh--circos-like-satellitome-density-plot)
- [satDNA_similarity.py](#satdna_similaritypy--automated-alignment-of-satellite-dna-monomers-for-variant-and-superfamily-analysis)
- [satFlank.sh](#satflanksh--satellitedna-flank-extraction-and-gene-overlap-analysis-pipeline)
- [gene_extractor.sh](#gene_extractorsh--extract-gene-sequences-from-ncbi-by-species)
- [monoMiner.py](#monominerpy--reference-guided-motif-mining-from-short-reads)
    
## Quick Start
```bash
git clone https://github.com/zenirodrigo/TandemRepeatEXplorer.git
cd T-REx
# T-REx provides a complete Conda environment containing all required dependencies for every module in the toolkit.

conda env create -f environment.yml
conda activate trex_env

#Alternatively, advanced users may install dependencies manually.

```

## arrayScope.sh – Genome Repeat Analysis Pipeline

### 📜 Description
**ArrayScope** is a **Bash + Python** pipeline designed to identify and visualize satellite DNA arrays in genome assemblies.  
It processes multiple genomes and reference monomers, runs **[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**, merges nearby repeat regions (< 2000 bp apart), and generates high-quality plots for chromosome-scale visualization.

> 💡 For variants in reference files — e.g., two slightly different monomers or the same repeat from different species, include an underscore `_` in the filename.

###  Features
- Processes **multiple genome assemblies** and reference monomer files
- Fully automated **BLAST** analysis with configurable parameters
- Merges nearby repeat hits (< 2000 bp apart)
- Generates three publication-quality visualizations:
  1. Chromosome-scale repeat distribution (linear map)
  2. Heatmap of total base pairs per chromosome/reference
  3. Scatter plot of array sizes by chromosome
- Handles `.fa`, `.fna`, `.fasta` formats
- Parallel processing using **[GNU parallel](https://www.gnu.org/software/parallel/)**

###  Dependencies
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

###  Usage
```bash
bash arrayScope.sh

```
### Exemple of results
<img width="1077" height="717" alt="chromosomes_with_annotations" src="https://github.com/user-attachments/assets/96fc8939-23a9-46fa-b957-52c160d372ae" />



---

---


## satDNA_density.sh – Circos-like satellitome density plot

### 📜 Description
**satDNA_density.sh** is a **Bash + Python** module for generating a **circos-like** overview of satDNA distributions across chromosomes/contigs from a genome assembly.

It runs **BLASTn** searches of satDNA references against a genome FASTA, merges nearby hits (≤ 2000 bp), and produces a circular plot with:

- **Outer ring:** chromosomes/contigs (preserves the original FASTA order)
- **Ruler:** tick marks every **10 Mb** and labels every **50 Mb**
- **Main ring:** stacked satDNA coverage per **100 kb** window (Top 10 colored + Other)
- **Inner track:** **GC fraction** per **100 kb** window
- **Top 10 selection mode (interactive):**
  1) **First 10 sequences in the reference FASTA (FASTA order)**  
  2) **Top 10 most abundant (by total bp covered in genome hits)**

>  IMPORTANT: contig/chromosome IDs are **not re-ordered and not re-indexed**.  
> The plot only applies a **safe abbreviation** for labels: `ChromosomeN` → `Chr N`.  
> All other IDs remain unchanged (e.g., `B1`, `Sex_chromosome`, `scaffold_1`, `unplaced` and any other name.).

###  Features
- Works with **multiple genome FASTA** files
- Uses only the **first N contigs/chromosomes** from each genome (preserving FASTA order)
- BLAST hits are merged within **2000 bp** only when on the **same chromosome AND same satDNA**
- Parallel BLAST execution using **GNU parallel**
- Produces both figures and tabular outputs (TSV) for downstream analysis

###  Dependencies
- **BLAST+** (`makeblastdb`, `blastn`)
- **GNU parallel**
- **Python 3.8+** with:
  - `numpy`
  - `pandas`
  - `matplotlib`

### 📂 Inputs
1. Genome FASTA file(s) (space-separated)
2. Number of chromosome sequences (contigs) to use from each genome FASTA
3. Reference FASTA file(s) (satDNA monomers OR a satellitome multi-FASTA)
4. Repeat multiplier (number of monomers to build arrays for BLAST sensitivity)
5. Number of parallel jobs
6. Top-10 selection mode (FASTA order or abundance)

### 📂 Outputs (per genome)
Inside a folder named after the genome FASTA (basename):
- `valid_monomers.bed`
- `satellitome_density_100kb_by_ref_long.tsv`
- `satellitome_density_100kb_total.tsv`
- `gc_track_100kb_full.tsv`
- `satellitome_density_circos_like.png`

###  Usage
```bash
bash satDNA_density.sh
```
### Exemple of results
<img width="3009" height="2835" alt="satelitome_density_circos_like" src="https://github.com/user-attachments/assets/4fe4dcf5-4540-44b4-bbda-8e11be003352" />


---

## satDNA_similarity.py – Automated alignment of satellite DNA monomers for variant and superfamily analysis

### 📜 Description
**satDNA_similarity.py** performs automated, biologically informed alignments of **satellite DNA (satDNA) monomers** to classify them at two hierarchical levels:

- **satDNA families** (true sequence variants of the same satellite)
- **satDNA superfamilies** (homologous but distinct satellites)

---

### Features
- Alignment-based comparison of satDNA monomers
- Explicit handling of **circular sequences** (no fixed start position)
- **Reverse-complement equivalence**
- **Insertions and deletions (gaps)** included in similarity estimates
- Automatic distinction between **families (variants)** and **superfamilies (homology only)**
- Supports **very long monomers** (including >5 kb)
- Generates alignment **proof files** for full biological traceability
- Produces strict FASTA outputs compatible with downstream tools

---

###  Dependencies
- **Python 3.8+**
- Standard Python libraries only  
(no external bioinformatics dependencies required)

---

### 📂 Input
- A FASTA file containing **satDNA monomers**  
  (strict FASTA format: `>ID` followed by sequence)

---

### 📂 Output Files
For an input file named `satellitome.fasta`, the script generates:

- `satellitome.id80.family_reps.fasta`  
  FASTA file containing **one representative sequence per satDNA family**

- `satellitome.id80.families.tsv`  
  Tabular file describing family membership

- `satellitome.id80.proof.txt`  
  Alignment proofs documenting **why sequences were grouped as variants**

- `satellitome.id80.superfamilies.tsv`  
  Pairwise superfamily relationships (homology without collapsing sequences)

---

### Usage
The script is fully interactive and requires only two inputs:
- A FASTA file containing satDNA monomers
- A minimum identity threshold (default: 0.80)
```bash
python3 satDNA_similarity.py
```


---

## gene_extractor.sh – Extract Gene Sequences from NCBI by Species

### 📜 Description
Automates the retrieval of genomic sequences for a given gene across multiple species using **[NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)** and **datasets CLI**.

###  Features
- Searches for Gene ID using gene symbol and species name
- Extracts genomic coordinates from NCBI
- Fallback to RNA FASTA if coordinates are missing
- Downloads sequences with `efetch`
- Reverse-complements negative strands automatically
- Outputs FASTA files named by gene and species

###  Dependencies
```bash
conda install -c bioconda entrez-direct datasets-cli seqkit
```
Also requires: unzip, cut, sed, find, head

### 📂 Input Files
- Gene symbol (e.g., BRCA1)
- Text file with one species per line

### 📂 Output Files
- `{GENE_SYMBOL}_{SPECIES_NAME}.fasta`

###  Usage
```bash
bash gene_extractor.sh BRCA1 species.txt
```

---

## satFlank.sh – SatelliteDNA Flank Extraction and Gene Overlap Analysis Pipeline

### 📜 Description
**satFlank.sh** identifies satellite DNA flanking regions using **BLAST**, extracts sequences, and compares with GFF annotations to find overlapping genes.

###  Features
- Multiplies reference sequence for improved BLAST sensitivity
- Filters `.bed` to keep largest non-overlapping regions
- Sequence extraction with **seqtk**
- Gene overlap analysis from GFF
- Parallel processing

###  Dependencies
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

###  Usage
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



---

## monoMiner.py – Reference-guided Motif Mining from Short Reads

### 📜 Description
monoMiner performs reference-guided motif discovery from sequencing libraries, enabling rapid detection and clustering of tandem repeat motifs across multiple datasets with parallel processing and automated filtering.

###  Dependencies
- **Python 3**
- **cd_hit_filter_size.py**
- Standard Python libraries
- FASTQ sequencing files
- Mapping file (TSV: species_code → species_name)
- **Install CD-HIT (required)**:

Follow instructions at [CD-HIT](https://github.com/weizhongli/cdhit) Official Repository


Ensure cd_hit_filter_size.py is in your PATH.


###  Usage
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
- Heng Li (for existing and inspiring generations of biologists and bioinformaticians), Ole Tange, the NCBI  and the BEDTools teams for creating  tools that keep bioinformatics alive and reproducible.
- Support from Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).
- These codes were developed during my PhD, supported by a scholarship from Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP).


## ⁉️ Contact
For questions or improvements, contact:
**rodrigo.zeni@unesp.br**
or
**rodrigo-zeni@outlook.com.br**


## 📖 Citation
If you use T-REx in your research, please cite the associated publication (in preparation).


