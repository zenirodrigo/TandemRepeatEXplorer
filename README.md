# genome-manipulation-satDNA

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Zeni_repgen.sh - Genome Repeat Analysis Pipeline

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
