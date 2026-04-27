#!/bin/bash

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

    awk -v m="$multiplier" '
        BEGIN { RS=">"; ORS=""; }
        NR>1 {
            n = split($0, a, "\n")
            header = a[1]

            seq = ""
            for (i=2; i<=n; i++) seq = seq a[i]

            gsub(/[ \t\r]/, "", seq)
            seq = toupper(seq)
            gsub(/-/, "", seq)
            gsub(/[^ACGTN]/, "", seq)

            printf(">%s\n", header)
            for (j=0; j<m; j++) printf("%s", seq)
            printf("\n")
        }
    ' "$in_fa" > "$out_fa"
}

run_blast_for_ref() {
    local ref_no_ext="$1"
    local temp_genome="$2"
    local multiplier="$3"
    local out_bed="$4"

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

    blastn -task blastn \
        -outfmt "6 qseqid sseqid sstart send" \
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

    awk -v OFS='\t' '
        {
            q=$1
            chr=$2
            s=($3<$4)?$3:$4
            e=($3<$4)?$4:$3
            print q, chr, s, e
        }
    ' "$blast_output" \
    | sort -k2,2 -k1,1 -k3,3n \
    | awk -v OFS='\t' -v dist=2000 '
        {
            q=$1
            chr=$2
            s=$3
            e=$4

            if (NR==1) {
                cq=q
                cchr=chr
                cs=s
                ce=e
            } else {
                if (q==cq && chr==cchr && s <= ce + dist) {
                    if (e > ce) ce=e
                } else {
                    print cchr, cs, ce, cq
                    cq=q
                    cchr=chr
                    cs=s
                    ce=e
                }
            }
        }
        END {
            if (NR>0) print cchr, cs, ce, cq
        }
    ' > "$out_bed"
}

export -f run_blast_for_ref remove_extensions make_multiplied_fasta

read -e -p "Enter genome file names (space-separated): " input_biblios
read -p "How many chromosome/scaffold sequences will be used? " num_sequences
read -e -p "Enter reference (satDNA or another tandem repeat MONOMER) files (space-separated): " refs_in
read -p "How many monomers will be used to create an array (minimum monomers in this study)? " multiplier
read -p "How many threads do you want to use? (e.g., 4, 8, etc.): " NUM_THREADS

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed. Use: conda install -c conda-forge parallel"
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
        /^>/ {
            if (count >= num_seq) exit
            count++
        }
        { print }
    ' "$input_biblio" > "$temp_genome"

    makeblastdb -in "$temp_genome" -dbtype nucl -out "$temp_genome" -parse_seqids

    tmp_parallel_dir="$genome_name/tmp_parallel_beds"
    rm -rf "$tmp_parallel_dir"
    mkdir -p "$tmp_parallel_dir"

    parallel --jobs "$NUM_THREADS" \
        run_blast_for_ref {} "$temp_genome" "$multiplier" "$tmp_parallel_dir/{}.bed" \
        ::: "${expanded_refs_no_ext[@]}"

    cat "$tmp_parallel_dir"/*.bed 2>/dev/null | sort -k1,1 -k2,2n > "$genome_name/valid_monomers.bed" || true

    if [ ! -s "$genome_name/valid_monomers.bed" ]; then
        echo "The file valid_monomers.bed is empty for $genome_name."
    else
        echo "The file valid_monomers.bed was created for $genome_name."
    fi

python3 - <<EOF
import os
import json
import re
from collections import OrderedDict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

fasta_file = "$temp_genome"
bed_file = "$genome_name/valid_monomers.bed"
genome_name = "$genome_name"

def natural_sort_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'([0-9]+)', str(text))]

def sanitize_accession(acc):
    acc = str(acc).strip()
    acc = acc.replace("ref|", "").replace("|", "")
    acc = acc.split()[0]
    return acc

def infer_pretty_name_from_header(header, fallback_index):
    h = str(header).strip()
    h_low = h.lower()

    m_chr = re.search(r'chromosome\\s+([0-9]+)', h_low)
    if m_chr:
        return f"Chromosome{m_chr.group(1)}"

    m_scaf = re.search(r'scaffold[_\\s]+([0-9]+)', h_low)
    if m_scaf:
        return f"scaffold{m_scaf.group(1)}"

    return f"Chromosome{fallback_index}"

def build_header_mapping(fasta_path):
    records_info = []
    accession_to_pretty = OrderedDict()
    pretty_to_length = OrderedDict()

    idx = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        idx += 1
        raw_header = record.description.strip()
        accession = sanitize_accession(record.id)

        pretty = infer_pretty_name_from_header(raw_header, idx)

        if pretty in pretty_to_length:
            pretty = f"{pretty}_{idx}"

        accession_to_pretty[accession] = pretty
        pretty_to_length[pretty] = len(record.seq)

        records_info.append({
            "index": idx,
            "accession": accession,
            "raw_header": raw_header,
            "pretty_name": pretty,
            "length": len(record.seq)
        })

    return accession_to_pretty, pretty_to_length, records_info

def read_bed(filepath):
    if (not os.path.exists(filepath)) or os.path.getsize(filepath) == 0:
        return pd.DataFrame(columns=["Chromosome","Start","End","Reference","References"])

    df = pd.read_csv(
        filepath,
        sep="\\t",
        header=None,
        names=["Chromosome", "Start", "End", "Reference"]
    )

    df["Chromosome"] = df["Chromosome"].astype(str).map(sanitize_accession)
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Chromosome", "Start", "End", "Reference"]).copy()
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["References"] = df["Reference"].apply(lambda x: [x] if pd.notna(x) else [])
    return df

def merge_intervals(df, dist=2000):
    merged = []
    df_sorted = df.sort_values(by=["Chromosome", "Start", "End"]).reset_index(drop=True)

    for chrom in df_sorted["Chromosome"].unique():
        chrom_data = df_sorted[df_sorted["Chromosome"] == chrom]
        current_start, current_end = None, None
        refs = set()

        for _, row in chrom_data.iterrows():
            s, e, rlist = int(row["Start"]), int(row["End"]), row["References"]

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

    return pd.DataFrame(merged, columns=["Chromosome", "Start", "End", "References"])

def save_header_mapping(records_info, out_tsv):
    pd.DataFrame(records_info).to_csv(out_tsv, sep="\\t", index=False)

accession_to_pretty, pretty_to_length, records_info = build_header_mapping(fasta_file)
save_header_mapping(records_info, os.path.join(genome_name, "sequence_name_mapping.tsv"))

df = read_bed(bed_file)

if df.empty:
    print("No BLAST hits found (valid_monomers.bed is empty). Skipping plots.")
    raise SystemExit(0)

df["Chromosome_original_accession"] = df["Chromosome"]
df["Chromosome"] = df["Chromosome"].map(lambda x: accession_to_pretty.get(x, x))

valid_names = set(pretty_to_length.keys())
df_filtered = df[df["Chromosome"].isin(valid_names)].copy()

if df_filtered.empty or len(pretty_to_length) == 0:
    print("No relevant chromosomes/scaffolds after standardization. Skipping plots.")
    raise SystemExit(0)

df_merged = merge_intervals(df_filtered)

if df_merged.empty:
    print("No arrays found after merge. Skipping plots.")
    raise SystemExit(0)

bed_refs = df_merged["References"].explode().dropna().unique()
all_refs = sorted(bed_refs, key=natural_sort_key)

if len(all_refs) == 0:
    print("No references found after merge. Skipping plots.")
    raise SystemExit(0)

color_file = os.path.join(genome_name, "reference_colors.json")

# High-contrast colors for the chromosome plot.
# This intentionally overwrites old pastel colors, because pastel colors made the arrays hard to see.
base_palette = []
base_palette.extend(sns.color_palette("tab10", 10))
base_palette.extend(sns.color_palette("Dark2", 8))
base_palette.extend(sns.color_palette("Set1", 9))

color_map = {}
for i, ref in enumerate(all_refs):
    rgb = base_palette[i % len(base_palette)]
    color_map[ref] = [float(rgb[0]), float(rgb[1]), float(rgb[2])]

with open(color_file, "w") as f:
    json.dump(color_map, f, indent=2)

color_map = {k: tuple(v) for k, v in color_map.items()}

df_merged.to_csv(os.path.join(genome_name, "merged_arrays.tsv"), sep="\\t", index=False)
df_filtered.to_csv(os.path.join(genome_name, "valid_monomers_mapped.tsv"), sep="\\t", index=False)

sorted_chromosomes = sorted(pretty_to_length.keys(), key=natural_sort_key)

reference_sizes = {}
for ref in all_refs:
    subset = df_merged[df_merged["References"].apply(lambda x: ref in x)]
    sizes_ = subset["End"] - subset["Start"]
    reference_sizes[ref] = (
        int(sizes_.min()) if not sizes_.empty else 0,
        int(sizes_.max()) if not sizes_.empty else 0
    )

# =========================================================
# PLOT 1: CHROMOSOME ANNOTATION PLOT - X-SCALED ARRAYS
# =========================================================
# Key interpretation:
#   - X position = genomic coordinate on the chromosome/scaffold.
#   - Rectangle width = array genomic interval, from Start to End.
#   - Rectangle height = only visual thickness; it does NOT encode array size.
#
# Output logic:
#   1) chromosomes_with_annotations_exact_scale.png/pdf
#      Every array has its exact genomic width. This is the biologically faithful version.
#      Very small arrays can be almost invisible on whole-chromosome scale.
#
#   2) chromosomes_with_annotations.png/pdf
#      Same genomic positions, but very small arrays receive a minimum visible display width.
#      This is the main readable version. Tiny arrays are visually exaggerated, while large
#      arrays remain true-scale when they are above the minimum display width.
#
# I am not saving a separate chromosomes_with_annotations_visible.* anymore, because it was
# identical to chromosomes_with_annotations.* and only created confusing duplicate outputs.
# =========================================================
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MultipleLocator


def merge_intervals_by_reference(df, dist=2000):
    rows = []
    tmp = df.copy()
    tmp["Reference"] = tmp["Reference"].astype(str)
    tmp = tmp.sort_values(["Chromosome", "Reference", "Start", "End"]).reset_index(drop=True)

    for (chrom, ref), sub in tmp.groupby(["Chromosome", "Reference"], sort=False):
        current_start = None
        current_end = None

        for _, row in sub.iterrows():
            s = int(row["Start"])
            e = int(row["End"])
            if s > e:
                s, e = e, s

            if current_start is None:
                current_start = s
                current_end = e
            elif s <= current_end + dist:
                current_end = max(current_end, e)
            else:
                rows.append([chrom, current_start, current_end, ref])
                current_start = s
                current_end = e

        if current_start is not None:
            rows.append([chrom, current_start, current_end, ref])

    out = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Reference"])
    if not out.empty:
        out["ArraySize"] = out["End"] - out["Start"] + 1
        out = out.sort_values(["Chromosome", "Reference", "Start", "End"]).reset_index(drop=True)
    return out


def draw_chromosome_array_plot(plot_arrays, out_prefix, enhanced_visibility=False, min_visible_kb=80):
    sorted_chromosomes_for_axis = list(sorted_chromosomes)[::-1]
    max_len_mb = max(pretty_to_length.values()) / 1e6

    chrom_height = 0.58
    array_height = 0.50
    spacing = 1.08
    fig_width = 38
    fig_height = max(14, len(sorted_chromosomes_for_axis) * 0.48)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    y_positions = {}

    for idx, chrom in enumerate(sorted_chromosomes_for_axis):
        y = idx * spacing
        y_positions[chrom] = y
        length_mb = pretty_to_length[chrom] / 1e6

        body = FancyBboxPatch(
            (0, y - chrom_height / 2),
            length_mb,
            chrom_height,
            boxstyle=f"round,pad=0.00,rounding_size={chrom_height * 0.22}",
            linewidth=0.40,
            edgecolor="#bdbdbd",
            facecolor="#f1f1f1",
            alpha=1.00,
            zorder=1,
        )
        ax.add_patch(body)

    plot_arrays_sorted = plot_arrays.sort_values("ArraySize", ascending=True).reset_index(drop=True)
    min_visible_mb = float(min_visible_kb) / 1000.0

    for _, row in plot_arrays_sorted.iterrows():
        chrom = row["Chromosome"]
        if chrom not in y_positions:
            continue

        ref = row["Reference"]
        start_bp = int(row["Start"])
        end_bp = int(row["End"])
        if start_bp > end_bp:
            start_bp, end_bp = end_bp, start_bp

        true_start_mb = start_bp / 1e6
        true_end_mb = end_bp / 1e6
        true_width_mb = max(true_end_mb - true_start_mb, 1e-9)

        if enhanced_visibility:
            display_width_mb = max(true_width_mb, min_visible_mb)
            midpoint_mb = (true_start_mb + true_end_mb) / 2.0
            start_mb = max(0.0, midpoint_mb - display_width_mb / 2.0)
            end_mb = min(pretty_to_length[chrom] / 1e6, midpoint_mb + display_width_mb / 2.0)
            width_mb = max(end_mb - start_mb, 1e-9)
        else:
            start_mb = true_start_mb
            width_mb = true_width_mb

        y = y_positions[chrom]
        color = color_map.get(ref, "#777777")

        rect = Rectangle(
            (start_mb, y - array_height / 2),
            width_mb,
            array_height,
            linewidth=0.18,
            edgecolor=color,
            facecolor=color,
            alpha=0.98,
            zorder=5,
        )
        ax.add_patch(rect)

    reference_sizes_for_plot = {}
    for ref in all_refs:
        subset = plot_arrays[plot_arrays["Reference"] == ref]
        sizes_ = subset["ArraySize"]
        reference_sizes_for_plot[ref] = (
            int(sizes_.min()) if not sizes_.empty else 0,
            int(sizes_.max()) if not sizes_.empty else 0,
        )

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=color_map[ref], edgecolor=color_map[ref], alpha=0.98)
        for ref in all_refs
    ]
    labels = [
        f"{ref} (min: {reference_sizes_for_plot[ref][0]:,} bp; max: {reference_sizes_for_plot[ref][1]:,} bp)"
        for ref in all_refs
    ]

    ax.set_yticks([y_positions[c] for c in sorted_chromosomes_for_axis])
    ax.set_yticklabels(sorted_chromosomes_for_axis, fontsize=10)
    ax.set_xlabel("Position on chromosome/scaffold (Mb)", fontsize=12)
    ax.set_ylabel("Chromosomes/scaffolds", fontsize=12)

    if enhanced_visibility:
        title_suffix = f"visible mode; arrays < {min_visible_kb} kb widened for display"
    else:
        title_suffix = "exact scale"
    ax.set_title(f"Distribution of tandem repeat arrays across chromosomes/scaffolds ({title_suffix})", fontsize=14, pad=12)

    ax.set_xlim(0, max_len_mb * 1.015)
    ax.set_ylim(-spacing, (len(sorted_chromosomes_for_axis) - 1) * spacing + spacing)

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.grid(axis="x", which="major", color="#d0d0d0", linewidth=0.45, alpha=0.65, zorder=0)
    ax.grid(axis="x", which="minor", color="#eeeeee", linewidth=0.25, alpha=0.55, zorder=0)
    ax.grid(axis="y", visible=False)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")

    ax.legend(
        handles,
        labels,
        title="References",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=10,
        frameon=True,
        borderaxespad=0.0,
    )

    fig.tight_layout()
    fig.savefig(os.path.join(genome_name, f"{out_prefix}.png"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(genome_name, f"{out_prefix}.pdf"), bbox_inches="tight")
    plt.close(fig)


plot_arrays = merge_intervals_by_reference(df_filtered, dist=2000)
plot_arrays.to_csv(os.path.join(genome_name, "arrays_by_reference_for_plot.tsv"), sep="\t", index=False)

if plot_arrays.empty:
    print("No arrays available for chromosome plot. Skipping chromosome plot.")
else:
    draw_chromosome_array_plot(
        plot_arrays=plot_arrays,
        out_prefix="chromosomes_with_annotations_exact_scale",
        enhanced_visibility=False,
        min_visible_kb=80,
    )
    draw_chromosome_array_plot(
        plot_arrays=plot_arrays,
        out_prefix="chromosomes_with_annotations",
        enhanced_visibility=True,
        min_visible_kb=80,
    )

# =========================================================
# PLOT 2: ARRAY SIZE VS CHROMOSOME/ SCAFFOLD SCATTER
# =========================================================
df_exploded = df_merged.explode("References")
df_exploded = df_exploded[df_exploded["References"].notna()].copy()
df_exploded["ArraySize"] = df_exploded["End"] - df_exploded["Start"] + 1
df_exploded["Chromosome"] = pd.Categorical(
    df_exploded["Chromosome"],
    categories=sorted_chromosomes,
    ordered=True
)

if df_exploded.empty:
    print("No data available for the scatter plot. Skipping scatter.")
    raise SystemExit(0)

plt.figure(figsize=(18, 8))
sns.stripplot(
    data=df_exploded,
    x="Chromosome",
    y="ArraySize",
    hue="References",
    hue_order=all_refs,
    palette=color_map,
    dodge=True,
    alpha=0.85,
    size=5
)
plt.xticks(rotation=90)
plt.xlabel("Sequences")
plt.ylabel("Array size (bp)")
plt.title("Array size distribution by chromosome/scaffold")
plt.legend(title="References", bbox_to_anchor=(1.01, 1), loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(genome_name, "array_chromosome_vs_size_scatter.png"), dpi=300, bbox_inches="tight")
plt.close()

print("Plots successfully generated in:", genome_name)
print("Saved image:", os.path.join(genome_name, "chromosomes_with_annotations.png"))
print("Saved image:", os.path.join(genome_name, "array_chromosome_vs_size_scatter.png"))
print("Saved table:", os.path.join(genome_name, "valid_monomers.bed"))
print("Saved table:", os.path.join(genome_name, "valid_monomers_mapped.tsv"))
print("Saved table:", os.path.join(genome_name, "merged_arrays.tsv"))
print("Saved table:", os.path.join(genome_name, "sequence_name_mapping.tsv"))
EOF

done


