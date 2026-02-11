#!/bin/bash

# Enable filename autocompletion
shopt -s progcomp

_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}
complete -F _autocomplete read

set -euo pipefail

remove_extensions() {
    local filename="$1"
    filename=$(basename "$filename")
    filename="${filename%.fasta}"
    filename="${filename%.fa}"
    filename="${filename%.fna}"
    echo "$filename"
}

make_multiplied_fasta() {
    local in_fa="$1"
    local multiplier="$2"
    local out_fa="$3"

    # Robust FASTA reader:
    # - Works with multi-line sequences
    # - Cleans sequence to ACGTN only (and removes hyphens)
    # - Writes one record per entry, with proper newlines
    awk -v m="$multiplier" '
        BEGIN { RS=">"; ORS=""; }
        NR>1 {
            n = split($0, a, "\n");
            header = a[1];
            seq = "";
            for (i=2; i<=n; i++) seq = seq a[i];

            gsub(/[ \t\r]/, "", seq);
            seq = toupper(seq);
            gsub(/-/, "", seq);
            gsub(/[^ACGTN]/, "", seq);

            printf(">%s\n", header);
            for (j=0; j<m; j++) printf("%s", seq);
            printf("\n");
        }
    ' "$in_fa" > "$out_fa"
}

run_blast_for_ref() {
    local ref_no_ext="$1"
    local temp_genome="$2"
    local multiplier="$3"
    local temp_bed="$4"

    local ref_fasta=""
    if [ -f "${ref_no_ext}.fasta" ]; then
        ref_fasta="${ref_no_ext}.fasta"
    elif [ -f "${ref_no_ext}.fa" ]; then
        ref_fasta="${ref_no_ext}.fa"
    elif [ -f "${ref_no_ext}.fna" ]; then
        ref_fasta="${ref_no_ext}.fna"
    else
        echo "Warning: No .fasta, .fa or .fna file found for ${ref_no_ext}"
        return 0
    fi

    local multiplied="multiplied_reference_${ref_no_ext}.fasta"
    make_multiplied_fasta "$ref_fasta" "$multiplier" "$multiplied"

    local blast_output="blast_${ref_no_ext}.out"

    # IMPORTANT:
    # - Each parallel job runs its own BLAST -> use -num_threads 1
    # - Control total CPU via GNU parallel --jobs
    blastn -task blastn \
        -outfmt "6" \
        -db "$temp_genome" \
        -query "$multiplied" \
        -out "$blast_output" \
        -evalue 1e-10 \
        -qcov_hsp_perc 70 \
        -num_threads 1

    if [ ! -s "$blast_output" ]; then
        echo "BLAST found no matches for $ref_no_ext"
        return 0
    fi

    # Merge nearby hits per chromosome (dist=2000)
    awk '{start=($9 < $10) ? $9 : $10; end=($9 < $10) ? $10 : $9; print $2, start, end}' "$blast_output" \
    | sort -k1,1 -k2,2n \
    | awk -v OFS='\t' -v dist=2000 '
        {
            if (NR == 1) { chr=$1; start=$2; end=$3 }
            else {
                if ($1 == chr && ($2 <= end + dist)) {
                    end = ($3 > end) ? $3 : end
                } else {
                    print chr, start, end
                    chr=$1; start=$2; end=$3
                }
            }
        }
        END { if (NR>0) print chr, start, end }
    ' \
    | awk -v ref="$ref_no_ext" '{split(ref, a, "_"); print $0"\t"a[1]}' >> "$temp_bed"
}
export -f run_blast_for_ref remove_extensions make_multiplied_fasta

read -e -p "Enter genome file names (space-separated): " input_biblios
read -p "How many chromosomes sequences will be used? " num_sequences
read -e -p "Enter reference (satDNA or another tandem repeat MONOMER) files (space-separated): " refs_in
read -p "How many monomers will be used to create a array (here is the minimum monomers to form a array in this study)? " multiplier
read -p "Quantas threads deseja utilizar? (ex: 4, 8, etc.): " NUM_THREADS

if ! command -v parallel &> /dev/null; then
    echo "Erro: GNU parallel não está instalado. Use: conda install -c conda-forge parallel"
    exit 1
fi

expanded_refs=()
for r in $refs_in; do
    matches=( $(compgen -G "$r") )
    if [ ${#matches[@]} -eq 0 ]; then
        expanded_refs+=( "$r" )
    else
        expanded_refs+=( "${matches[@]}" )
    fi
done

expanded_refs_no_ext=()
for ref in "${expanded_refs[@]}"; do
    expanded_refs_no_ext+=( "$(remove_extensions "$ref")" )
done

echo "Final references (no extension): ${expanded_refs_no_ext[@]}"

for input_biblio in $input_biblios; do
    genome_name=$(remove_extensions "$input_biblio")
    mkdir -p "$genome_name"

    temp_genome="$genome_name/temp_genome.fasta"
    awk -v num_seq="$num_sequences" '
        BEGIN { count = 0 }
        /^>/ { if (count >= num_seq) exit; count++ }
        { print }
    ' "$input_biblio" > "$temp_genome"

    makeblastdb -in "$temp_genome" -dbtype nucl -out "$temp_genome" -parse_seqids

    temp_bed="$genome_name/valid_monomers_temp.bed"
    > "$temp_bed"

    parallel --jobs "$NUM_THREADS" run_blast_for_ref {} "$temp_genome" "$multiplier" "$temp_bed" ::: "${expanded_refs_no_ext[@]}"

    mv "$temp_bed" "$genome_name/valid_monomers.bed"
    echo "Arquivo valid_monomers.bed criado para $genome_name."

python3 - <<EOF
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import re

def read_bed(filepath):
    if (not os.path.exists(filepath)) or os.path.getsize(filepath) == 0:
        return pd.DataFrame(columns=['Chromosome','Start','End','Reference','References'])
    df = pd.read_csv(filepath, sep='\\t', header=None, names=['Chromosome', 'Start', 'End', 'Reference'])
    df['References'] = df['Reference'].apply(lambda x: [x] if pd.notna(x) else [])
    return df

def get_chromosome_lengths(fasta_path, chromosomes):
    lengths = {}
    chrom_set = set(chromosomes)
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in chrom_set:
            lengths[record.id] = len(record.seq)
    return lengths

def natural_sort_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', str(text))]

def merge_intervals(df, dist=2000):
    merged = []
    df_sorted = df.sort_values(by=['Chromosome', 'Start', 'End']).reset_index(drop=True)
    for chrom in df_sorted['Chromosome'].unique():
        chrom_data = df_sorted[df_sorted['Chromosome'] == chrom]
        current_start, current_end = None, None
        refs = set()
        for _, row in chrom_data.iterrows():
            s, e, rlist = int(row['Start']), int(row['End']), row['References']
            if current_start is None:
                current_start, current_end = s, e
                refs.update(rlist)
            elif s <= current_end + dist:
                current_end = max(current_end, e)
                refs.update(rlist)
            else:
                merged.append([chrom, current_start, current_end, list(refs)])
                current_start, current_end = s, e
                refs = set(rlist)
        if current_start is not None:
            merged.append([chrom, current_start, current_end, list(refs)])
    return pd.DataFrame(merged, columns=['Chromosome', 'Start', 'End', 'References'])

fasta_file = "$temp_genome"
bed_file = "$genome_name/valid_monomers.bed"

df = read_bed(bed_file)
if df.empty:
    print("No BLAST hits found (valid_monomers.bed vazio). Pulando plots.")
    raise SystemExit(0)

relevant_chromosomes = df['Chromosome'].unique()
chromosome_lengths = get_chromosome_lengths(fasta_file, relevant_chromosomes)
df_filtered = df[df['Chromosome'].isin(chromosome_lengths)]

if df_filtered.empty or not chromosome_lengths:
    print("Sem cromossomos relevantes após filtragem. Pulando plots.")
    raise SystemExit(0)

df_merged = merge_intervals(df_filtered)

bed_refs = df_merged['References'].explode().dropna().unique()
all_refs = sorted(bed_refs, key=natural_sort_key)

if len(all_refs) == 0:
    print("Sem referências após merge. Pulando plots.")
    raise SystemExit(0)

color_file = os.path.join("$genome_name", "reference_colors.json")

if os.path.exists(color_file):
    with open(color_file, 'r') as f:
        color_map = json.load(f)
else:
    color_map = {}

missing_refs = [r for r in all_refs if r not in color_map]
if missing_refs:
    # Paleta suave/pastel por padrão
    new_palette = sns.color_palette("pastel", len(missing_refs))
    for i, ref in enumerate(missing_refs):
        rgb = new_palette[i]
        color_map[ref] = tuple(rgb)

with open(color_file, 'w') as f:
    json.dump(color_map, f)

color_map = {k: tuple(v) for k, v in color_map.items()}

reference_sizes = {}
for ref in all_refs:
    subset = df_merged[df_merged['References'].apply(lambda x: ref in x)]
    sizes_ = subset['End'] - subset['Start']
    reference_sizes[ref] = (int(sizes_.min()) if not sizes_.empty else 0, int(sizes_.max()) if not sizes_.empty else 0)

sorted_chromosomes = sorted(chromosome_lengths.keys(), key=natural_sort_key)

# Plot 1: Chromosomes with annotations
plt.figure(figsize=(36, 24))
spacing = 1.5
for idx, chrom in enumerate(sorted_chromosomes):
    length = chromosome_lengths[chrom]
    length_mb = length / 1e6
    y_pos = idx * spacing
    plt.plot([0, length_mb], [y_pos, y_pos], color='lightgray', linewidth=3)
    plt.fill_between([0, length_mb], y_pos - 0.4, y_pos + 0.4, color='lightgray', alpha=0.6)

    chrom_data = df_merged[df_merged['Chromosome'] == chrom]
    for _, row in chrom_data.iterrows():
        start_mb = row['Start'] / 1e6
        end_mb = row['End'] / 1e6
        refs = row['References']
        if not refs:
            continue
        segment_width = (end_mb - start_mb) / len(refs)
        for i, ref in enumerate(refs):
            c = color_map.get(ref, "#cccccc")
            plt.fill_between(
                [start_mb + i*segment_width, start_mb + (i+1)*segment_width],
                y_pos - 0.4,
                y_pos + 0.4,
                color=c,
                alpha=0.8
            )

handles = [plt.Line2D([0], [0], color=color_map[ref], lw=5) for ref in all_refs]
labels  = [f"{ref} (min: {reference_sizes[ref][0]:,.0f} bp, max: {reference_sizes[ref][1]:,.0f} bp)" for ref in all_refs]
plt.legend(handles, labels, title="References", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
plt.yticks([i * spacing for i in range(len(sorted_chromosomes))], sorted_chromosomes, fontsize=12)
plt.xlabel('Base Pairs (Mb)', fontsize=14)
plt.ylabel('Chromosome', fontsize=14)
plt.title('Satellite Array Distribution Across Chromosomes', fontsize=16)
plt.xlim(0, max(chromosome_lengths.values())/1e6)
plt.grid(False)
plt.tight_layout()
plt.savefig("$genome_name/chromosomes_with_annotations.png", dpi=300, bbox_inches='tight')
plt.close()

# Prepare exploded data
df_exploded = df_merged.explode('References')
df_exploded = df_exploded[df_exploded['References'].notna()].copy()
df_exploded['Size'] = df_exploded['End'] - df_exploded['Start']

if df_exploded.empty:
    print("Dados explodidos vazios. Pulando heatmap e scatter.")
    raise SystemExit(0)

# Plot 2: Heatmap
heatmap_data = df_exploded.groupby(['Chromosome', 'References'])['Size'].sum().unstack()
heatmap_data = heatmap_data.reindex(index=sorted_chromosomes, columns=all_refs).fillna(0)

if heatmap_data.size == 0 or heatmap_data.to_numpy().sum() == 0:
    print("Heatmap sem dados (tudo zero). Pulando heatmap.")
else:
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'Total Base Pairs (bp)'})
    plt.title('Total Base Pairs per Chromosome and Reference', fontsize=16)
    plt.xlabel('Reference', fontsize=14)
    plt.ylabel('Chromosome', fontsize=14)
    plt.xticks(rotation=60)
    plt.tight_layout()
    plt.savefig("$genome_name/array_frequency_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()

# Plot 3: Scatter
df_exploded['Chromosome'] = pd.Categorical(df_exploded['Chromosome'], categories=sorted_chromosomes, ordered=True)
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
plt.close()

print("Visualization plots saved.")
EOF

done
