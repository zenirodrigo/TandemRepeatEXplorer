[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Conda environment](https://img.shields.io/badge/Conda-environment.yml-44A833?logo=anaconda&logoColor=white)](environment.yml)
[![Python ≥3.10](https://img.shields.io/badge/Python-3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![R ≥4.2](https://img.shields.io/badge/R-%E2%89%A54.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bash](https://img.shields.io/badge/Bash-pipelines-4EAA25?logo=gnubash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey?logo=linux&logoColor=white)](#installation)

<img align="right" src="social-preview.png?raw=1" width="250" alt="Repository Icon">
<br clear="right">
<h1>Tandem Repeat Explorer</h1>


Tandem Repeat Explorer (T-REx) is a modular Bash/Python toolkit for the identification, characterization, and visualization of tandem arrays from genome assemblies and short-read sequencing data.

- **arrayScope.sh** – Characterizing and locating tandem repeat arrays in assembled genomes.
- **satDNA_density.sh** - Generating a circus overview of satDNAs density and distributions across a  genome assembly.
- **satDNA_similarity.py** – Biological clustering of satellite DNA monomers into families and superfamilies, accounting for circularity, reverse-complement, and indels.
- **satFlank.sh** – Studying the neighborhood of arrays using the assembled genome and its annotation.
- **monoMiner.py** – Automated pipeline for identifying biological motifs in sequencing libraries, with parallel processing and filtering.
- **gene_extractor.sh** – Automates the retrieval and extraction of the genomic sequence of a specific gene across multiple species using NCBI tools.
- **run_satDNA_synteny.sh** – Investigates the relationship between satellite DNA (satDNA) arrays and conserved syntenic blocks identified by MCScan/JCVI.
  
---
## Table of Contents
- [arrayScope.sh](#arrayscopesh--genome-repeat-analysis-pipeline)
- [satDNA_density.sh](#satdna_densitysh--circos-like-satellitome-density-plot)
- [satDNA_similarity.py](#satdna_similaritypy--automated-alignment-of-satellite-dna-monomers-for-variant-and-superfamily-analysis)
- [gene_extractor.sh](#gene_extractorsh--extract-gene-sequences-from-ncbi-by-species)
- [monoMiner.py](#monominerpy--reference-guided-motif-mining-from-short-reads)
- [run_satDNA_synteny.sh](#run_satdna_syntenysh--relationship-between-collinear-blocks-and-satdna-arrays)
## Quick Start
```bash
git clone https://github.com/zenirodrigo/TandemRepeatEXplorer.git
cd TandemRepeatEXplorer
# T-REx provides a complete Conda environment containing all required dependencies for every module in the toolkit.

conda env create -f environment.yml
conda activate trex_env

#Alternatively, if you will use only one tool, you can install dependencies manually.

```

## arrayScope.sh – Genome Repeat Analysis Pipeline

### Description
**ArrayScope** is a **Bash + Python** pipeline designed to identify, characterize and visualize tandem repeat arrays in genome assemblies.

The pipeline processes multiple genome assemblies and reference monomers, performs **BLAST+** searches, automatically groups neighboring hits into arrays using an adaptive distance approach, and generates chromosome-scale visualizations and quantitative summaries.

ArrayScope was designed to work with both individual satDNA monomers and complete satellitome datasets.

### Features

- Processes **multiple genome assemblies** and reference monomer files
- Accepts single satDNA monomers or complete satellitome datasets
- Fully automated **BLAST** analysis with configurable parameters
- Adaptive array construction based on observed genomic spacing between hits
- Automatic estimation of array boundaries
- Generates publication-quality chromosome-scale visualizations
- Produces detailed tabular summaries of arrays
- Handles `.fa`, `.fna`, and `.fasta` formats
- Parallel processing using **GNU Parallel**
- Suitable for large genome assemblies
  
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
ArrayScope requires the following software.

#### Bioinformatics software

```bash
conda install -c bioconda blast bedtools
conda install -c conda-forge parallel
```

#### Python packages

```bash
conda install -c conda-forge \
    python=3.11 \
    biopython \
    numpy \
    pandas \
    matplotlib \
    openpyxl
```

or, alternatively,

```bash
pip install biopython numpy pandas matplotlib openpyxl
```

### Input Files
1. Genome assemblies (FASTA format)
2. Reference monomers (FASTA format)

OBS: If you have a fasta with multi satDNAs, you can run the split_satdna.sh before running the arrayScope, to separate each sequence in a diferent fasta.
### Output Files

Main outputs include:

#### BLAST and array tables

- raw_blast_hits.tsv
- raw_blast_hits_mapped.tsv
- valid_monomers.bed
- arrays_by_reference.tsv
- merged_regions_multi_satdna.tsv
- adaptive_merge_distance_by_reference.tsv

#### Quantitative summaries

- summary_by_reference.tsv
- heatmap_total_array_bp_by_chromosome_reference.tsv
- heatmap_array_count_by_chromosome_reference.tsv
- satdna_metrics_by_chromosome.tsv
- satdna_metrics_by_chromosome_reference.tsv
- satellitome_composition.tsv
- Excel workbook containing all summary tables

#### Figures

- chromosomes_with_annotations_exact_scale.png
- chromosomes_with_annotations.png
- chromosomes_with_annotations_plus_top_N.png
- array_chromosome_vs_size_scatter.png
- heatmap_chromosome_vs_satdna_total_bp.png
- satdna_percentage_composition_by_chromosome.png

#### Auxiliary files

- reference_colors.json
- sequence_name_mapping.tsv

PNG and PDF versions are generated for all figures.


###  Usage
```bash
bash arrayScope.sh

```
### Exemple of results
<img width="22759" height="8626" alt="prochilodus_lineatus_assembly_final_10x_chromosomes_with_annotations" src="https://github.com/user-attachments/assets/cd5d8cda-8aa7-4a8c-bff7-fc87b635e31b" />
<img width="5375" height="2365" alt="prochilodus_lineatus_assembly_final_10x_array_chromosome_vs_size_scatter" src="https://github.com/user-attachments/assets/2130f82c-834b-421c-b52a-888bb5556a2e" />
![Uploading prochilodus_lineatus_assembly_final_10x_chromosomes_with_annotations.png…]()
<img width="17394" height="7239" alt="prochilodus_lineatus_assembly_final_10x_satdna_percentage_composition_by_chromosome" src="https://github.com/user-attachments/assets/8d6be067-a8af-4a52-b2e5-ecb43bd39cb7" />
<img width="2831" height="2924" alt="prochilodus_lineatus_assembly_final_10x_heatmap_chromosome_vs_satdna_total_bp" src="https://github.com/user-attachments/assets/1c243184-2d61-4d23-9527-e5e46344b10e" />



---

---


## satDNA_density.sh – Circos-like satellitome density plot

### Description
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

### Inputs
1. Genome FASTA file(s) (space-separated)
2. Number of chromosome sequences (contigs) to use from each genome FASTA
3. Reference FASTA file(s) (satDNA monomers OR a satellitome multi-FASTA)
4. Repeat multiplier (number of monomers to build arrays for BLAST sensitivity)
5. Number of parallel jobs
6. Top-10 selection mode (FASTA order or abundance)

### Outputs (per genome)
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

### Description

**satDNA_similarity.py** is a Python-based toolkit for clustering satellite DNA monomers into biologically meaningful families.

The program performs circular-aware sequence comparisons, automatically handling:

- Circular phase shifts
- Reverse-complement relationships
- Insertions and deletions (indels)
- Monomer length variation

Families are built using a complete-linkage strategy, preventing unrelated monomers from being merged through transitive similarity chains.

The program was specifically designed for comparative satellitome analyses involving hundreds to thousands of satDNA monomers from multiple species.
---

### Features

- Circular sequence comparison
- Reverse-complement equivalence
- Gap-aware similarity estimation
- Complete-linkage family clustering
- Family representative selection
- Automatic generation of family FASTA files
- Circular phase correction relative to family representatives
- Pairwise representative-member alignments
- Family-wide representative-anchored alignments
- Scalable to large comparative satellitome datasets
---

###  Dependencies
- **Python 3.8+**
- Standard Python libraries only  
(no external bioinformatics dependencies required)

---

### Input
- A FASTA file containing **satDNA monomers**  
  (strict FASTA format: `>ID` followed by sequence)

---

### Output Files
### Family-level Outputs
Each family generates a set of files including:

- Family_000001_AmeSat01.members.fasta
 
Contains all original monomers assigned to the family.

- Family_000001_AmeSat01.frame_corrected_monomers.fasta

All family members are rotated and oriented relative to the family representative while preserving their original monomer length.

- Family_000001_AmeSat01.pairwise_to_first_representative.fasta

Contains pairwise alignments between the family representative and each family member.

- Family_000001_AmeSat01.firstseq_reference_anchored_alignment.fasta

Family-wide alignment generated after circular phase correction, suitable for visualization in Geneious, AliView, Jalview and related software.


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

### Description
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

### Input Files
- Gene symbol (e.g., BRCA1)
- Text file with one species per line

### Output Files
- `{GENE_SYMBOL}_{SPECIES_NAME}.fasta`

###  Usage
```bash
bash gene_extractor.sh BRCA1 species.txt
```

---


## monoMiner.py – Reference-guided Motif Mining from Short Reads

### Description
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

## Output Files
.

├── final.fasta                    # Concatenated motifs

├── output.tsv                     # Annotated sequences

├── reference_motifs_*.fasta       # Per-library motif files

└── final.fasta.nr0.*.sel.fasta    # CD-HIT filtered results




---
## run_satDNA_synteny.sh – Relationship between collinear blocks and satDNA arrays

**satDNA-Synteny** is a lightweight pipeline for investigating the relationship between satellite DNA (satDNA) arrays and conserved collinear (syntenic) blocks identified from whole-genome comparisons.

The pipeline converts gene-based collinear blocks into genomic coordinates, identifies satDNA loci located inside or outside conserved blocks, calculates the distance from each satDNA array to the nearest collinear block, and summarizes the results using comparative heatmaps.

---

## Requirements

- Python ≥ 3.8
- BEDTools ≥ 2.30
- R ≥ 4.2

### Required R packages

```r
install.packages(c("ggplot2","dplyr"))
```

---

## Input files

For each genome, the pipeline requires three input files.

### 1. Genome annotation (BED)

A BED file containing the genomic coordinates of annotated genes used during the MCScan/JCVI analysis.

**Example**

```text
Chr1    12000    13500    Gene0001
Chr1    18200    19450    Gene0002
...
```

---

### 2. MCScan/JCVI anchors.simple

The `anchors.simple` file generated by MCScan/JCVI describing each collinear block.

**Example**

```text
Gene0001    Gene0258
Gene0002    Gene0259
...
```

The script `blocks_with_id.py` converts these gene pairs into genomic coordinates.

> **Important**
>
> Depending on the orientation of the MCScan/JCVI comparison, the genes corresponding to the first and last positions of each block may be located in different columns of the `anchors.simple` file.
>
> Users should verify which columns correspond to the block boundaries and modify the following lines in `blocks_with_id.py` if necessary:
>
> ```python
> g_start = cols[0]
> g_end   = cols[1]
> ```
>
> This adjustment is only required when the column order differs from the default expected by the script.

---

### 3. satDNA BED

A BED file containing genomic coordinates of satellite DNA arrays.

This file can be generated by **arrayScope.sh**, although **any standard four-column BED file** is accepted.

**Required format**

```text
chromosome    start    end    satDNA_family
```

**Example**

```text
Chr1    420000    424300    Sat01
Chr1    812000    814500    Sat02
...
```

---

## Running the pipeline

Edit the species list inside `run_satDNA_synteny.sh`.

```bash
species=(
Species1
Species2
Species3
Species4
)
```

Then execute

```bash
chmod +x run_satDNA_synteny.sh

./run_satDNA_synteny.sh
```

---

## Pipeline steps

### Step 1 — Build genomic coordinates of collinear blocks

Convert MCScan/JCVI collinear blocks into genomic coordinates.

**Input**

```
Species.bed
Species.anchors.simple
```

**Output**

```
Species.blocks.id.bed
```

---

### Step 2 — Identify satDNA loci inside and outside conserved blocks

Uses BEDTools to determine whether each satDNA locus overlaps a conserved collinear block.

Commands

```bash
bedtools intersect
```

Outputs

```
Species_inside.bed
Species_outside.bed
```

---

### Step 3 — Calculate distances to the nearest collinear block

Computes the distance between each satDNA locus and the nearest conserved collinear block.

Command

```bash
bedtools closest
```

Output

```
Species.dist.txt
```

The plotting script automatically supports both standard BEDTools outputs (8 or 9 columns).

---

### Step 4 — Generate comparative heatmaps

Produces two summary figures:

- Percentage of satDNA loci located **outside** conserved collinear blocks.
- Percentage of satDNA loci located **more than 100 kb** from the nearest conserved block.

Outputs

```
heatmap_outside.pdf
heatmap_far.pdf
```

---

## Output files

For each genome

```
Species1.Species2.blocks.id.bed
Species_inside.bed
Species_outside.bed
Species.dist.txt
```

Final figures

```
heatmap_outside.pdf
heatmap_far.pdf
```

---

## Notes

- This pipeline accepts **any four-column satDNA BED file** and is **not restricted to ArrayScope outputs**.
- Collinear blocks may originate from any pairwise whole-genome comparison performed with MCScan/JCVI.
- The only assumption made by the pipeline concerns the columns used to define block boundaries in `anchors.simple`. Depending on the orientation of the synteny comparison, users may need to modify the corresponding lines in `blocks_with_id.py`.


---

## License

This project is licensed under the MIT License

## Acknowledgments
- Heng Li (for existing and inspiring generations of biologists and bioinformaticians), Ole Tange, the NCBI  and the BEDTools teams for creating  tools that keep bioinformatics alive and reproducible.
- Support from Universidade Estadual Paulista Júlio de Mesquita Filho (UNESP).
- These codes were developed during my PhD, supported by a scholarship from Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP).


##  Contact
For questions or improvements, contact:
**rodrigo.zeni@unesp.br**
or
**rodrigo-zeni@outlook.com.br**


##  Citation
If you use T-REx in your research, please cite:
dos Santos, R. Z. Tandem Repeat Explorer (T-REx) [Computer software]. https://github.com/zenirodrigo/TandemRepeatEXplorer


