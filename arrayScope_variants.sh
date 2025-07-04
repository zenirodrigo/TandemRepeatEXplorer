#!/bin/bash

# Enable filename autocompletion
shopt -s progcomp

# Autocomplete function
_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}

# Assign autocomplete to read command
complete -F _autocomplete read

# ------------------------------------------------
# Function to remove .fa, .fna, .fasta extensions
# ------------------------------------------------
remove_extensions() {
    local filename="$1"
    filename=$(basename "$filename")
    filename="${filename%.fasta}"
    filename="${filename%.fa}"
    filename="${filename%.fna}"
    echo "$filename"
}

# Get input genomes
read -e -p "Enter genome file names (space-separated): " input_biblios

# Get number of sequences to use
read -p "How many chromosomes sequences will be used? " num_sequences

# Get reference files and multiplier
read -e -p "Enter reference (satDNA or another tandem repeat MONOMER) files (space-separated): " refs_in
read -p "How many monomers will be used to create a array (here is the minimum monomers to form a array in this study? " multiplier

# -------------------------------------------------
# 1) Remove extensions from reference files
# -------------------------------------------------
expanded_refs=()
for r in $refs_in; do
    # Expand wildcards (e.g., C*)
    matches=( $(compgen -G "$r") )
    if [ ${#matches[@]} -eq 0 ]; then
        expanded_refs+=( "$r" )
    else
        expanded_refs+=( "${matches[@]}" )
    fi
done

# Remove extensions from all reference names
expanded_refs_no_ext=()
for ref in "${expanded_refs[@]}"; do
    expanded_refs_no_ext+=( "$(remove_extensions "$ref")" )
done

# Display cleaned references
echo "Final references (no extension): ${expanded_refs_no_ext[@]}"

# -------------------------------------------------
# 2) Process each genome
# -------------------------------------------------
for input_biblio in $input_biblios; do
    # Create genome-named directory
    genome_name=$(remove_extensions "$input_biblio")
    if [ -d "$genome_name" ]; then
        echo "Error: Directory $genome_name already exists."
        exit 1
    fi
    mkdir -p "$genome_name" || { echo "Failed to create directory $genome_name"; exit 1; }

    # Create temp file with first N sequences
    temp_genome="$genome_name/temp_genome.fasta"
    awk -v num_seq="$num_sequences" '
    BEGIN { count = 0 }
    /^>/ { if (count >= num_seq) exit; count++ }
    { print }
    ' "$input_biblio" > "$temp_genome"

    if [ ! -f "$temp_genome" ]; then
        echo "Error: Temporary file $temp_genome not created."
        exit 1
    fi

    # Define BLAST output name
    blast_output="$genome_name/${genome_name}_blast"

    # Create temp BED file
    temp_bed="$genome_name/valid_monomers_temp.bed"
    > "$temp_bed"

    # Process each reference
    for ref_no_ext in "${expanded_refs_no_ext[@]}"; do
        # Find actual FASTA file
        if [ -f "${ref_no_ext}.fasta" ]; then
            ref_fasta="${ref_no_ext}.fasta"
        elif [ -f "${ref_no_ext}.fa" ]; then
            ref_fasta="${ref_no_ext}.fa"
        elif [ -f "${ref_no_ext}.fna" ]; then
            ref_fasta="${ref_no_ext}.fna"
        else
            echo "Warning: No .fasta, .fa or .fna file found for $ref_no_ext"
            continue
        fi

        # Multiply reference sequence
        awk 'NR % 2 == 0 { for (i=0;i<'"$multiplier"';i++) printf $0 } NR % 2 != 0' "$ref_fasta" > multiplied_reference.fasta

        # Prepare BLAST database
        makeblastdb -in "$temp_genome" -dbtype nucl

        # Run BLASTn
        blastn -task blastn -outfmt "6" -db "$temp_genome" -query "multiplied_reference.fasta" -out "$blast_output" -evalue 1e-10 -qcov_hsp_perc 50

        # Check BLAST results
        if [ ! -s "$blast_output" ]; then
            echo "BLAST found no matches for $ref_no_ext in $genome_name."
            continue
        fi

        # Process valid monomers
        awk '{start=($9 < $10) ? $9 : $10; end=($9 < $10) ? $10 : $9; print $2, start, end}' "$blast_output" \
        | sort -k1,1 -k2,2n \
        | awk -v OFS='\t' -v dist=2000 '
          {
              if (NR == 1) {
                  chr=$1; start=$2; end=$3
              } else {
                  if ($1 == chr && ($2 <= end + dist)) {
                      end = ($3 > end) ? $3 : end
                  } else {
                      print chr, start, end
                      chr=$1; start=$2; end=$3
                  }
              }
          }
          END { print chr, start, end }
        ' \
        | awk -v ref="$ref_no_ext" '{split(ref, a, "_"); print $0"\t"a[1]}' >> "$temp_bed"
    done

    # Finalize BED file
    mv "$temp_bed" "$genome_name/valid_monomers.bed"

    if [ ! -f "$genome_name/valid_monomers.bed" ]; then
        echo "Error: BED file $genome_name/valid_monomers.bed not created."
        exit 1
    fi

    echo "File valid_monomers.bed generated for genome $genome_name."

    # -------------------------------------------------------
    # Python visualization script (MODIFICADO PARA CORES CONSISTENTES)
    # -------------------------------------------------------
python3 - <<EOF
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import re
import hashlib
import random

# 1) Read BED file
def read_bed(file):
    df = pd.read_csv(file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Reference'])
    df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
    df['End'] = pd.to_numeric(df['End'], errors='coerce')
    return df

# 2) Natural sort function
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

# 3) Get chromosome lengths
def get_chromosome_lengths(fasta_file, relevant_chromosomes):
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in relevant_chromosomes:
            lengths[record.id] = len(record.seq)
    return lengths

# 4) Merge intervals
def merge_intervals(df):
    merged = []
    for chrom, group in df.groupby('Chromosome'):
        sorted_group = group.sort_values(by='Start')
        start, end = sorted_group.iloc[0]['Start'], sorted_group.iloc[0]['End']
        refs = {sorted_group.iloc[0]['Reference']}
        for _, row in sorted_group.iloc[1:].iterrows():
            if row['Start'] <= end + 2000:
                end = max(end, row['End'])
                refs.add(row['Reference'])
            else:
                merged.append({'Chromosome': chrom, 'Start': start, 'End': end, 'References': list(refs)})
                start, end = row['Start'], row['End']
                refs = {row['Reference']}
        merged.append({'Chromosome': chrom, 'Start': start, 'End': end, 'References': list(refs)})
    return pd.DataFrame(merged)

# 5) Color configuration using hash for all references
distinct_colors = [
    (0.90, 0.12, 0.12),
    (0.12, 0.55, 0.12),
    (0.12, 0.12, 0.90),
    (0.90, 0.75, 0.12),
    (0.54, 0.17, 0.89),
    (0.89, 0.47, 0.20),
    (0.20, 0.89, 0.67),
    (0.20, 0.67, 0.89),
    (0.67, 0.20, 0.89),
    (0.89, 0.20, 0.67),
    (0.0, 0.5, 0.0),
    (0.5, 0.0, 0.5),
    (0.0, 0.8, 0.8),
    (0.8, 0.0, 0.0),
    (0.3, 0.7, 0.9),
    (0.9, 0.6, 0.0),
    (0.4, 0.2, 0.6),
    (0.7, 0.3, 0.3),
    (0.1, 0.9, 0.1),
    (0.6, 0.4, 0.8),
    (0.95, 0.9, 0.1),
    (0.8, 0.5, 0.9),
    (0.5, 0.5, 0.5),
    (0.0, 0.3, 0.6),
    (0.9, 0.7, 0.4),
    (0.2, 0.8, 0.2),
    (0.7, 0.0, 0.7),
    (0.4, 0.6, 0.2),
    (0.9, 0.2, 0.5),
    (0.3, 0.4, 0.9),
    (0.6, 0.1, 0.1),
    (0.1, 0.6, 0.6),
    (0.9, 0.4, 0.6),
    (0.5, 0.3, 0.0),
    (0.8, 0.8, 0.0),
    (0.0, 0.7, 0.3),
    (0.7, 0.5, 0.0),
    (0.29, 0.0, 0.51),
    (1.0, 0.41, 0.71),
    (0.0, 0.0, 0.4),
    (0.8, 0.2, 0.8),
    (0.4, 0.0, 0.0),
    (0.6, 0.8, 0.2),
    (0.2, 0.2, 0.2),
    (0.94, 0.94, 0.86),
    (0.5, 0.0, 0.0),
    (0.82, 0.71, 0.55),
    (0.9, 0.0, 0.9),
    (0.0, 0.5, 0.5)
]

#num_extra_colors = 40
#extra_colors = sns.color_palette("husl", num_extra_colors)

big_palette = distinct_colors #+ extra_colors
#random.shuffle(big_palette[len(distinct_colors):])
num_colors = len(big_palette)

def get_color_for_ref(ref):
    # Generate consistent color based on reference name hash
    h = hashlib.md5(ref.encode('utf-8')).hexdigest()
    idx = int(h, 16) % num_colors
    return big_palette[idx]

# Input paths
fasta_file = "$temp_genome"
bed_file   = "$genome_name/valid_monomers.bed"
df = read_bed(bed_file)

# Chromosome processing
relevant_chromosomes = df['Chromosome'].unique()
chromosome_lengths   = get_chromosome_lengths(fasta_file, relevant_chromosomes)

df_filtered = df[df['Chromosome'].isin(relevant_chromosomes)]
filtered_chromosomes = {chrom: length for chrom, length in chromosome_lengths.items() if chrom in relevant_chromosomes}

df_filtered = df_filtered[df_filtered['Chromosome'].isin(filtered_chromosomes.keys())]
df_merged   = merge_intervals(df_filtered)

# Reference processing
bed_refs = df_merged['References'].explode().unique()
all_refs = sorted(bed_refs, key=natural_sort_key)
color_map = {r: get_color_for_ref(r) for r in all_refs}

# Size calculations
reference_sizes = {}
for ref in all_refs:
    subset = df_merged[df_merged['References'].apply(lambda x: ref in x)]
    if subset.empty:
        reference_sizes[ref] = (0, 0)
    else:
        sizes_ = subset['End'] - subset['Start']
        reference_sizes[ref] = (sizes_.min(), sizes_.max())

sorted_chromosomes = sorted(filtered_chromosomes.keys(), key=natural_sort_key)

# ---------------------------
# FIGURE 1: Chromosome plot
# ---------------------------
plt.figure(figsize=(36, 24))
spacing = 1.5

for idx, chrom in enumerate(sorted_chromosomes):
    length = filtered_chromosomes[chrom]
    length_mb = length / 1e6
    y_position = idx * spacing
    
    # Chromosome baseline
    plt.plot([0, length_mb], [y_position, y_position], color='lightgray', linewidth=3)
    plt.fill_between([0, length_mb], y_position - 0.4, y_position + 0.4, color='lightgray', alpha=0.6)
    
    chrom_data = df_merged[df_merged['Chromosome'] == chrom]
    for _, row in chrom_data.iterrows():
        start_mb = row['Start'] / 1e6
        end_mb   = row['End']   / 1e6
        refs = row['References']
        segment_width = (end_mb - start_mb) / len(refs) if len(refs) else 0
        for i, ref in enumerate(refs):
            c = color_map.get(ref, "#cccccc")
            plt.fill_between(
                [start_mb + i*segment_width, start_mb + (i+1)*segment_width],
                y_position - 0.4,
                y_position + 0.4,
                color=c,
                alpha=0.8
            )

# Legend
handles = [plt.Line2D([0], [0], color=color_map[ref], lw=5) for ref in all_refs]
labels  = [f"{ref} (min: {reference_sizes[ref][0]:,.0f} bp, max: {reference_sizes[ref][1]:,.0f} bp)" for ref in all_refs]

plt.legend(handles, labels, title="References", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.yticks([i * spacing for i in range(len(sorted_chromosomes))], sorted_chromosomes, fontsize=12)
plt.xlabel('Base Pairs (Mb)', fontsize=14)
plt.ylabel('Chromosome', fontsize=14)
plt.title('Satellite Array Distribution Across Chromosomes', fontsize=16)
if filtered_chromosomes:
    plt.xlim(0, max(filtered_chromosomes.values())/1e6)
plt.grid(False)
plt.tight_layout()
plt.savefig("$genome_name/chromosomes_with_annotations.png", dpi=300, bbox_inches='tight')
plt.savefig("$genome_name/chromosomes_with_annotations.pdf", format='pdf', dpi=300, bbox_inches='tight')

# --------------------------
# FIGURE 2: Heatmap
# --------------------------
df_exploded = df_merged.explode('References')
df_exploded['Size'] = df_exploded['End'] - df_exploded['Start']

heatmap_data = df_exploded.groupby(['Chromosome', 'References'])['Size'].sum().unstack()
heatmap_data = heatmap_data.reindex(index=sorted_chromosomes, columns=all_refs)

plt.figure(figsize=(12, 8))
sns.heatmap(
    heatmap_data,
    cmap='viridis',
    cbar_kws={'label': 'Total Base Pairs (bp)'}
)
plt.title('Total Base Pairs per Chromosome and Reference', fontsize=16)
plt.xlabel('Reference', fontsize=14)
plt.ylabel('Chromosome', fontsize=14)
plt.xticks(rotation=60)
plt.tight_layout()
plt.savefig("$genome_name/array_frequency_heatmap.png", dpi=300, bbox_inches='tight')

# --------------------------
# FIGURE 3: Scatter plot
# --------------------------
df_exploded['Chromosome'] = pd.Categorical(
    df_exploded['Chromosome'],
    categories=sorted_chromosomes,
    ordered=True
)

plt.figure(figsize=(12, 8))
sns.scatterplot(
    data=df_exploded,
    x='Chromosome',
    y='Size',
    hue='References',
    hue_order=all_refs,
    palette=color_map,
    alpha=0.7,
    s=100
)
plt.title('Relationship Between Chromosome and Array Size', fontsize=16)
plt.xlabel('Chromosome', fontsize=14)
plt.ylabel('Array Size (bp)', fontsize=14)
plt.xticks(rotation=60)
plt.legend(title="References", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig("$genome_name/array_chromosome_vs_size_scatter.png", dpi=300, bbox_inches='tight')

print("Visualization plots saved.")
EOF
done
