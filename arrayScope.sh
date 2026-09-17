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
    local out_hits_tsv="$4"
    local blast_threads="$5"
    local min_qcov="$6"

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

    blastn \
        -task blastn \
        -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
        -db "$temp_genome" \
        -query "$multiplied" \
        -out "$blast_output" \
        -evalue 1e-10 \
        -qcov_hsp_perc "$min_qcov" \
        -max_target_seqs 10000 \
        -dust no \
        -soft_masking false \
        -num_threads "$blast_threads"

    if [ ! -s "$blast_output" ]; then
        echo "BLAST found no matches for $ref_no_ext"
        return 0
    fi

    awk -v OFS='\t' -v ref="$ref_no_ext" '
        BEGIN {
            print "Reference","Query","Chromosome","Start","End","Strand","Pident","HitLength","QueryLength","Qstart","Qend","Evalue","Bitscore"
        }
        {
            q=$1
            chr=$2
            pident=$3
            hitlen=$4
            qlen=$5
            qstart=$6
            qend=$7
            sstart=$8
            send=$9
            evalue=$10
            bitscore=$11

            if (sstart <= send) {
                s=sstart
                e=send
                strand="+"
            } else {
                s=send
                e=sstart
                strand="-"
            }

            print ref,q,chr,s,e,strand,pident,hitlen,qlen,qstart,qend,evalue,bitscore
        }
    ' "$blast_output" > "$out_hits_tsv"
}

export -f run_blast_for_ref remove_extensions make_multiplied_fasta

# Enhanced user interface with better descriptions
echo "=== SatDNA Analysis Pipeline ==="
echo "This pipeline detects and characterizes satellite DNA arrays in genome assemblies."
echo ""

read -e -p "Enter genome file names (space-separated): " input_biblios
read -p "How many chromosome/scaffold sequences will be used? " num_sequences
read -e -p "Enter reference (satDNA or tandem repeat MONOMER) files (space-separated): " refs_in
read -p "How many monomers will be used to create an array (minimum monomers)? " multiplier
read -p "Maximum gap allowed for adaptive merging, in bp [default: 10000]: " adaptive_gap_cap
read -p "Fallback gap if adaptive distance cannot be estimated, in bp [default: 2000]: " fallback_gap
read -p "Minimum BLAST percent identity to keep a hit [default: 80]: " min_pident
read -p "Minimum BLAST query coverage to keep a hit [default: 70]: " min_qcov
read -p "How many largest arrays per satDNA should be highlighted in the TOP plot? [default: 2]: " top_n_arrays
read -p "Comma-separated chromosome/scaffold names to highlight (example: ChrB,B,microB): " highlight_chromosomes
read -p "How many threads do you want to use? (e.g., 4, 8, etc.): " NUM_THREADS

adaptive_gap_cap="${adaptive_gap_cap:-10000}"
fallback_gap="${fallback_gap:-2000}"
min_pident="${min_pident:-80}"
min_qcov="${min_qcov:-70}"
top_n_arrays="${top_n_arrays:-2}"

if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed. Use: conda install -c conda-forge parallel"
    exit 1
fi

if ! command -v blastn &> /dev/null; then
    echo "Error: blastn is not available in PATH."
    exit 1
fi

if ! command -v makeblastdb &> /dev/null; then
    echo "Error: makeblastdb is not available in PATH."
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
    output_prefix="${genome_name}_${multiplier}x"
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

    tmp_parallel_dir="$genome_name/tmp_parallel_hits"
    rm -rf "$tmp_parallel_dir"
    mkdir -p "$tmp_parallel_dir"

    # Run BLAST with additional coverage filter
    parallel --jobs "$NUM_THREADS" \
        run_blast_for_ref {} "$temp_genome" "$multiplier" "$tmp_parallel_dir/{}.raw_hits.tsv" 1 "$min_qcov" \
        ::: "${expanded_refs_no_ext[@]}"

    {
        first=1
        for f in "$tmp_parallel_dir"/*.raw_hits.tsv; do
            [ -e "$f" ] || continue
            if [ "$first" -eq 1 ]; then
                cat "$f"
                first=0
            else
                tail -n +2 "$f"
            fi
        done
    } > "$genome_name/${output_prefix}_raw_blast_hits.tsv"

    if [ ! -s "$genome_name/${output_prefix}_raw_blast_hits.tsv" ]; then
        echo "The file ${output_prefix}_raw_blast_hits.tsv is empty for $genome_name."
    else
        echo "The file ${output_prefix}_raw_blast_hits.tsv was created for $genome_name."
    fi

python3 - <<EOF
import os
import json
import re
import math
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# Improved color palettes
COLOR_PALETTES = {
    'vibrant': ['#E63946', '#457B9D', '#2A9D8F', '#F4A261', '#E76F51', '#6D597A', '#B5838D', '#4A7C59'],
    'pastel': ['#FFB3BA', '#FFDFBA', '#FFFFBA', '#BAFFC9', '#BAE1FF', '#E8BAFF', '#FFB3E6', '#B3FFE5'],
    'nature': ['#2E4057', '#6B8C5E', '#A8C69F', '#D4A373', '#E8D5B7', '#8CB0C9', '#B5C9D9', '#C9B8A8'],
    'scientific': ['#1B4965', '#5FA8D3', '#90C9E2', '#C2E3E0', '#F7D6E0', '#F2B5D4', '#E5989B', '#B56576'],
}

FASTA_FILE = "$temp_genome"
GENOME_NAME = "$genome_name"
OUTPUT_PREFIX = "$output_prefix"
RAW_HITS_FILE = os.path.join(GENOME_NAME, f"{OUTPUT_PREFIX}_raw_blast_hits.tsv")

MIN_MONOMERS_PER_ARRAY = int("$multiplier")
ADAPTIVE_GAP_CAP = int("$adaptive_gap_cap")
FALLBACK_GAP = int("$fallback_gap")
MIN_PIDENT = float("$min_pident")
TOP_N_ARRAYS = int("$top_n_arrays")
HIGHLIGHT_CHROMOSOMES_RAW = "$highlight_chromosomes"
SELECTED_PALETTE = 'scientific'  # Change to 'vibrant', 'pastel', 'nature', or 'scientific'

os.makedirs(GENOME_NAME, exist_ok=True)

def output_path(filename):
    return os.path.join(GENOME_NAME, f"{OUTPUT_PREFIX}_{filename}")

def natural_sort_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'([0-9]+)', str(text))]

def sanitize_accession(acc):
    acc = str(acc).strip().split()[0]
    # BLAST+ can report a FASTA ID as gb|CM038724.1| while Biopython
    # reads the same FASTA record as CM038724.1. Strip only the namespace.
    match = re.fullmatch(r'(?:gb|ref|emb|dbj|lcl)\|([^|]+)\|?', acc)
    return match.group(1) if match else acc

def infer_pretty_name_from_header(header, fallback_index):
    h = str(header).strip()
    h_low = h.lower()

    patterns = [
        r'chromosome\s+([0-9A-Za-z_.-]+)',
        r'chr(?:omosome)?[_\s-]*([0-9A-Za-z_.-]+)',
        r'scaffold[_\s-]+([0-9A-Za-z_.-]+)',
    ]

    for pat in patterns:
        m = re.search(pat, h_low)
        if m:
            label = m.group(1)
            if pat.startswith("scaffold"):
                return f"scaffold{label}"
            return f"Chromosome{label}"

    if re.fullmatch(r'[A-Za-z0-9_.-]+', h):
        return h

    return f"Chromosome{fallback_index}"

def build_header_mapping(fasta_path):
    records_info = []
    accession_to_pretty = OrderedDict()
    pretty_to_length = OrderedDict()
    accession_to_length = OrderedDict()

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
        accession_to_length[accession] = len(record.seq)

        records_info.append({
            "index": idx,
            "accession": accession,
            "raw_header": raw_header,
            "pretty_name": pretty,
            "length": len(record.seq)
        })

    return accession_to_pretty, pretty_to_length, accession_to_length, records_info

def read_raw_hits(path):
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return df

    expected = [
        "Reference","Query","Chromosome","Start","End","Strand","Pident",
        "HitLength","QueryLength","Qstart","Qend","Evalue","Bitscore"
    ]
    for col in expected:
        if col not in df.columns:
            raise ValueError(f"Missing expected column in raw BLAST table: {col}")

    df["Chromosome"] = df["Chromosome"].astype(str).map(sanitize_accession)
    df["Reference"] = df["Reference"].astype(str)
    df["Query"] = df["Query"].astype(str)

    numeric_cols = ["Start","End","Pident","HitLength","QueryLength","Qstart","Qend","Evalue","Bitscore"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["Chromosome","Start","End","Reference","Pident","HitLength"]).copy()
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["HitLength"] = df["HitLength"].astype(int)
    df["QueryLength"] = df["QueryLength"].astype(int)

    df = df[df["End"] >= df["Start"]].copy()

    if MIN_PIDENT > 0:
        df = df[df["Pident"] >= MIN_PIDENT].copy()

    return df.reset_index(drop=True)

def save_header_mapping(records_info, out_tsv):
    pd.DataFrame(records_info).to_csv(out_tsv, sep="\t", index=False)

def compute_gap_stats_for_reference(sub):
    gaps = []
    tmp = sub.sort_values(["Chromosome", "Start", "End"]).reset_index(drop=True)

    for chrom, chrom_data in tmp.groupby("Chromosome", sort=False):
        chrom_data = chrom_data.sort_values(["Start", "End"]).reset_index(drop=True)
        prev_end = None

        for _, row in chrom_data.iterrows():
            s = int(row["Start"])
            e = int(row["End"])
            if prev_end is not None:
                gap = s - prev_end - 1
                if gap > 0:
                    gaps.append(gap)
            prev_end = max(prev_end, e) if prev_end is not None else e

    gaps = np.array(gaps, dtype=float)
    if len(gaps) == 0:
        return {
            "n_gaps": 0,
            "adaptive_gap": FALLBACK_GAP,
            "gap_median": np.nan,
            "gap_p75": np.nan,
            "gap_p90": np.nan,
            "gap_p95": np.nan,
            "method": "fallback_no_positive_gaps"
        }

    local_gaps = gaps[gaps <= ADAPTIVE_GAP_CAP]
    if len(local_gaps) < 5:
        return {
            "n_gaps": int(len(gaps)),
            "adaptive_gap": FALLBACK_GAP,
            "gap_median": float(np.median(gaps)),
            "gap_p75": float(np.percentile(gaps, 75)),
            "gap_p90": float(np.percentile(gaps, 90)),
            "gap_p95": float(np.percentile(gaps, 95)),
            "method": "fallback_few_local_gaps"
        }

    gap_p95 = float(np.percentile(local_gaps, 95))
    adaptive_gap = int(math.ceil(min(max(gap_p95, 1), ADAPTIVE_GAP_CAP)))

    return {
        "n_gaps": int(len(gaps)),
        "adaptive_gap": adaptive_gap,
        "gap_median": float(np.median(local_gaps)),
        "gap_p75": float(np.percentile(local_gaps, 75)),
        "gap_p90": float(np.percentile(local_gaps, 90)),
        "gap_p95": gap_p95,
        "method": "p95_local_gaps_capped"
    }

def create_arrays_by_reference(df):
    arrays = []
    factor_rows = []

    for ref, ref_data in df.groupby("Reference", sort=False):
        stats = compute_gap_stats_for_reference(ref_data)
        merge_gap = int(stats["adaptive_gap"])
        factor_rows.append({
            "Reference": ref,
            "AdaptiveMergeGap_bp": merge_gap,
            "FallbackGap_bp": FALLBACK_GAP,
            "GapCap_bp": ADAPTIVE_GAP_CAP,
            "NumPositiveGaps": stats["n_gaps"],
            "GapMedian_local_bp": stats["gap_median"],
            "GapP75_local_bp": stats["gap_p75"],
            "GapP90_local_bp": stats["gap_p90"],
            "GapP95_local_bp": stats["gap_p95"],
            "Method": stats["method"]
        })

        ref_data = ref_data.sort_values(["Chromosome", "Start", "End"]).reset_index(drop=True)

        for chrom, chrom_data in ref_data.groupby("Chromosome", sort=False):
            chrom_data = chrom_data.sort_values(["Start", "End"]).reset_index(drop=True)

            current = None
            for _, row in chrom_data.iterrows():
                s = int(row["Start"])
                e = int(row["End"])

                if current is None:
                    current = {
                        "Chromosome": chrom,
                        "Start": s,
                        "End": e,
                        "Reference": ref,
                        "NumMonomers": 1,
                        "HitLengths": [int(row["HitLength"])],
                        "Pidents": [float(row["Pident"])],
                        "Bitscores": [float(row["Bitscore"])],
                        "Strands": [str(row["Strand"])],
                        "Gaps": [],
                        "AdaptiveMergeGap_bp": merge_gap
                    }
                    continue

                gap = s - current["End"] - 1

                if gap <= merge_gap:
                    if gap > 0:
                        current["Gaps"].append(int(gap))
                    current["End"] = max(current["End"], e)
                    current["NumMonomers"] += 1
                    current["HitLengths"].append(int(row["HitLength"]))
                    current["Pidents"].append(float(row["Pident"]))
                    current["Bitscores"].append(float(row["Bitscore"]))
                    current["Strands"].append(str(row["Strand"]))
                else:
                    arrays.append(current)
                    current = {
                        "Chromosome": chrom,
                        "Start": s,
                        "End": e,
                        "Reference": ref,
                        "NumMonomers": 1,
                        "HitLengths": [int(row["HitLength"])],
                        "Pidents": [float(row["Pident"])],
                        "Bitscores": [float(row["Bitscore"])],
                        "Strands": [str(row["Strand"])],
                        "Gaps": [],
                        "AdaptiveMergeGap_bp": merge_gap
                    }

            if current is not None:
                arrays.append(current)

    rows = []
    for arr in arrays:
        num_plus = arr["Strands"].count("+")
        num_minus = arr["Strands"].count("-")
        strand_mix = "mixed" if (num_plus > 0 and num_minus > 0) else ("+" if num_plus > 0 else "-")

        # Calculate orientation consistency
        orientation_consistency = max(num_plus, num_minus) / len(arr["Strands"]) if arr["Strands"] else 0

        gaps = arr["Gaps"]
        rows.append({
            "Chromosome": arr["Chromosome"],
            "Start": int(arr["Start"]),
            "End": int(arr["End"]),
            "Reference": arr["Reference"],
            "ArraySize": int(arr["End"] - arr["Start"] + 1),
            "NumMonomers": int(arr["NumMonomers"]),
            "MeanMonomerHitLength": float(np.mean(arr["HitLengths"])) if arr["HitLengths"] else np.nan,
            "MedianMonomerHitLength": float(np.median(arr["HitLengths"])) if arr["HitLengths"] else np.nan,
            "MeanPident": float(np.mean(arr["Pidents"])) if arr["Pidents"] else np.nan,
            "MedianPident": float(np.median(arr["Pidents"])) if arr["Pidents"] else np.nan,
            "MeanBitscore": float(np.mean(arr["Bitscores"])) if arr["Bitscores"] else np.nan,
            "MedianGap": float(np.median(gaps)) if gaps else 0,
            "MaxGap": int(max(gaps)) if gaps else 0,
            "AdaptiveMergeGap_bp": int(arr["AdaptiveMergeGap_bp"]),
            "PlusStrandHits": int(num_plus),
            "MinusStrandHits": int(num_minus),
            "StrandMix": strand_mix,
            "OrientationConsistency": orientation_consistency
        })

    out = pd.DataFrame(rows)
    factors = pd.DataFrame(factor_rows)

    if not out.empty:
        out = out[out["NumMonomers"] >= MIN_MONOMERS_PER_ARRAY].copy()
        out = out.sort_values(["Chromosome","Reference","Start","End"]).reset_index(drop=True)

    return out, factors

def create_merged_regions_multi_satdna(arrays_by_reference, overlap_dist=0):
    if arrays_by_reference.empty:
        return pd.DataFrame(columns=["Chromosome","Start","End","References","NumReferences","ArraySize"])

    rows = []
    tmp = arrays_by_reference.sort_values(["Chromosome","Start","End"]).reset_index(drop=True)

    for chrom, sub in tmp.groupby("Chromosome", sort=False):
        current_start = None
        current_end = None
        refs = set()

        for _, row in sub.iterrows():
            s = int(row["Start"])
            e = int(row["End"])
            ref = str(row["Reference"])

            if current_start is None:
                current_start = s
                current_end = e
                refs = {ref}
            elif s <= current_end + overlap_dist:
                current_end = max(current_end, e)
                refs.add(ref)
            else:
                rows.append({
                    "Chromosome": chrom,
                    "Start": current_start,
                    "End": current_end,
                    "References": ",".join(sorted(refs, key=natural_sort_key)),
                    "NumReferences": len(refs),
                    "ArraySize": current_end - current_start + 1
                })
                current_start = s
                current_end = e
                refs = {ref}

        if current_start is not None:
            rows.append({
                "Chromosome": chrom,
                "Start": current_start,
                "End": current_end,
                "References": ",".join(sorted(refs, key=natural_sort_key)),
                "NumReferences": len(refs),
                "ArraySize": current_end - current_start + 1
            })

    return pd.DataFrame(rows)

def make_valid_monomers_bed(raw_hits_mapped, out_path):
    if raw_hits_mapped.empty:
        pd.DataFrame(columns=["Chromosome","Start","End","Reference"]).to_csv(out_path, sep="\t", header=False, index=False)
        return
    raw_hits_mapped[["Chromosome","Start","End","Reference"]].sort_values(
        ["Chromosome","Start","End","Reference"]
    ).to_csv(out_path, sep="\t", header=False, index=False)

def get_color_map(refs):
    """Get color map using professional palettes with better contrast."""
    colors = []
    palette = COLOR_PALETTES.get(SELECTED_PALETTE, COLOR_PALETTES['scientific'])
    
    # Expand palette if more references than colors
    if len(refs) <= len(palette):
        colors = palette[:len(refs)]
    else:
        # Cycle through palette with slight variations
        extended = []
        for i in range(len(refs)):
            extended.append(palette[i % len(palette)])
        colors = extended
    
    # Add some variation for readability
    color_map = {}
    for i, ref in enumerate(refs):
        if i < len(colors):
            color_map[ref] = colors[i]
        else:
            # Fallback
            color_map[ref] = '#888888'
    
    return color_map

def parse_highlight_chromosomes(raw, sorted_chromosomes):
    if raw.strip() == "":
        return set()

    requested = [x.strip() for x in raw.split(",") if x.strip()]
    requested_lower = {x.lower() for x in requested}
    result = set()

    for chrom in sorted_chromosomes:
        if chrom.lower() in requested_lower:
            result.add(chrom)
        for req in requested_lower:
            if req and req in chrom.lower():
                result.add(chrom)

    return result

def draw_chromosome_array_plot(plot_arrays, pretty_to_length, sorted_chromosomes, all_refs, color_map, out_prefix, enhanced_visibility=False, min_visible_kb=80, title_extra="", top_arrays=None, combine_with_top=False, highlight_chromosomes=None):
    """Improved chromosome array plot with better aesthetics."""
    highlight_chromosomes = highlight_chromosomes or set()

    sorted_chromosomes_for_axis = list(sorted_chromosomes)[::-1]
    max_len_mb = max(pretty_to_length.values()) / 1e6

    # Enhanced aesthetics: better sizes and spacing
    chrom_height = 0.62
    array_height = 0.48
    spacing = 1.10
    fig_width = 38
    fig_height = max(14, len(sorted_chromosomes_for_axis) * 0.55)

    if combine_with_top:
        fig, axes = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(fig_width, fig_height * 1.65),
            sharex=True,
            gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.15}
        )
        ax_list = list(axes)
        datasets = [
            (plot_arrays, "All arrays"),
            (top_arrays if top_arrays is not None else plot_arrays, f"Top {TOP_N_ARRAYS} largest arrays per satDNA")
        ]
    else:
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax_list = [ax]
        datasets = [(plot_arrays, "All arrays")]

    min_visible_mb = float(min_visible_kb) / 1000.0

    for ax, (dataset, panel_title) in zip(ax_list, datasets):
        y_positions = {}

        for idx, chrom in enumerate(sorted_chromosomes_for_axis):
            y = idx * spacing
            y_positions[chrom] = y
            length_mb = pretty_to_length[chrom] / 1e6

            if chrom in highlight_chromosomes:
                facecolor = "#FFF8E1"  # Soft yellow highlight
                edgecolor = "#D4A017"
                linewidth = 2.0
            else:
                # Gradient based on chromosome length for visual interest
                intensity = min(0.95, 0.75 + 0.2 * (length_mb / max_len_mb))
                facecolor = f"#{int(245*intensity):02x}{int(245*intensity):02x}{int(245*intensity):02x}"
                edgecolor = "#B0B0B0"
                linewidth = 0.6

            body = FancyBboxPatch(
                (0, y - chrom_height / 2),
                length_mb,
                chrom_height,
                boxstyle=f"round,pad=0.00,rounding_size={chrom_height * 0.20}",
                linewidth=linewidth,
                edgecolor=edgecolor,
                facecolor=facecolor,
                alpha=1.00,
                zorder=1,
            )
            ax.add_patch(body)
            
            # Add chromosome length label
            ax.text(length_mb + max_len_mb * 0.005, y, 
                   f'{length_mb:.1f} Mb', 
                   fontsize=7, va='center', ha='left', color='#666666')

        dataset_sorted = dataset.sort_values("ArraySize", ascending=True).reset_index(drop=True)

        # Add subtle shadow for arrays
        for _, row in dataset_sorted.iterrows():
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
            
            # Parse color for transparency
            try:
                rgb = to_rgb(color)
                face_alpha = 0.85
                edge_alpha = 0.95
            except:
                rgb = (0.5, 0.5, 0.5)
                face_alpha = 0.85
                edge_alpha = 0.95

            # Main array rectangle
            rect = Rectangle(
                (start_mb, y - array_height / 2),
                width_mb,
                array_height,
                linewidth=0.5,
                edgecolor=color,
                facecolor=color,
                alpha=face_alpha,
                zorder=5,
            )
            ax.add_patch(rect)
            
            # Add small array size label for large arrays
            if width_mb > 1.0:  # Only label if visible
                size_kb = (end_bp - start_bp + 1) / 1000
                if size_kb > 50:
                    ax.text(start_mb + width_mb/2, y, 
                           f'{size_kb:.0f}kb', 
                           fontsize=6, va='center', ha='center', 
                           color='white', weight='bold', alpha=0.9)

        ax.set_yticks([y_positions[c] for c in sorted_chromosomes_for_axis])
        ax.set_yticklabels(sorted_chromosomes_for_axis, fontsize=10)
        ax.set_ylabel("Chromosomes/scaffolds", fontsize=12, weight='bold')

        if enhanced_visibility:
            title_suffix = f"visible mode; arrays < {min_visible_kb} kb widened for display"
        else:
            title_suffix = "exact scale"

        ax.set_title(f"{panel_title} ({title_suffix}){title_extra}", fontsize=14, pad=15, weight='bold')
        ax.set_xlim(0, max_len_mb * 1.08)
        ax.set_ylim(-spacing * 0.7, (len(sorted_chromosomes_for_axis) - 1) * spacing + spacing * 0.7)

        # Better grid styling
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(2))
        ax.grid(axis="x", which="major", color="#D0D0D0", linewidth=0.5, alpha=0.7, zorder=0)
        ax.grid(axis="x", which="minor", color="#EEEEEE", linewidth=0.2, alpha=0.5, zorder=0)
        ax.grid(axis="y", visible=False)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.spines["left"].set_color("#888888")
        ax.spines["bottom"].set_color("#888888")
        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)

    ax_list[-1].set_xlabel("Position on chromosome/scaffold (Mb)", fontsize=12, weight='bold')

    # Enhanced legend
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=color_map[ref], edgecolor=color_map[ref], alpha=0.9)
        for ref in all_refs
    ]
    labels = [str(ref) for ref in all_refs]

    ax_list[0].legend(
        handles,
        labels,
        title="SatDNA References",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=10,
        frameon=True,
        borderaxespad=0.0,
        ncol=1,
        framealpha=0.9,
        edgecolor='#888888'
    )

    fig.tight_layout()
    fig.savefig(output_path(f"{out_prefix}.png"), dpi=600, bbox_inches="tight", facecolor='white')
    fig.savefig(output_path(f"{out_prefix}.pdf"), bbox_inches="tight", facecolor='white')
    plt.close(fig)

def draw_scatter(arrays_by_reference, sorted_chromosomes, all_refs, color_map):
    """Improved scatter plot with better aesthetics."""
    df = arrays_by_reference.copy()
    if df.empty:
        print("No data available for the scatter plot. Skipping scatter.")
        return

    df["Chromosome"] = pd.Categorical(df["Chromosome"], categories=sorted_chromosomes, ordered=True)
    df = df.sort_values(["Chromosome","Reference","ArraySize"]).reset_index(drop=True)

    fig_width = max(18, len(sorted_chromosomes) * 0.45)
    fig_height = 8
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ref_offsets = {}
    if len(all_refs) == 1:
        ref_offsets[all_refs[0]] = 0
    else:
        offsets = np.linspace(-0.35, 0.35, len(all_refs))
        ref_offsets = {ref: offsets[i] for i, ref in enumerate(all_refs)}

    chrom_to_x = {chrom: i for i, chrom in enumerate(sorted_chromosomes)}

    for ref in all_refs:
        sub = df[df["Reference"] == ref].copy()
        if sub.empty:
            continue
        xs = [chrom_to_x[c] + ref_offsets[ref] for c in sub["Chromosome"]]
        
        # Size based on number of monomers
        sizes = [max(20, min(200, m * 0.5 + 10)) for m in sub["NumMonomers"]]
        
        ax.scatter(
            xs,
            sub["ArraySize"],
            s=sizes,
            alpha=0.7,
            color=color_map.get(ref, "#777777"),
            label=ref,
            edgecolors='white',
            linewidth=0.5
        )

    ax.set_xticks(range(len(sorted_chromosomes)))
    ax.set_xticklabels(sorted_chromosomes, rotation=90, fontsize=9)
    ax.set_xlabel("Sequences", fontsize=12, weight='bold')
    ax.set_ylabel("Array size (bp)", fontsize=12, weight='bold')
    ax.set_title("Array size distribution by chromosome/scaffold", fontsize=14, weight='bold')
    
    # Better grid
    ax.grid(axis="y", which="major", color="#D0D0D0", linewidth=0.5, alpha=0.6)
    ax.grid(axis="y", which="minor", color="#EEEEEE", linewidth=0.2, alpha=0.4)
    
    # Use log scale if data spans multiple orders of magnitude
    sizes = df["ArraySize"].values
    if len(sizes) > 0 and sizes.max() / max(sizes.min(), 1) > 100:
        ax.set_yscale('log')
        ax.set_ylabel("Array size (bp) - Log scale", fontsize=12, weight='bold')

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")

    ax.legend(title="References", bbox_to_anchor=(1.01, 1), loc="upper left", 
              fontsize=9, title_fontsize=10, frameon=True, framealpha=0.9)
    
    fig.tight_layout()
    fig.savefig(output_path("array_chromosome_vs_size_scatter.png"), dpi=300, bbox_inches="tight", facecolor='white')
    fig.savefig(output_path("array_chromosome_vs_size_scatter.pdf"), bbox_inches="tight", facecolor='white')
    plt.close(fig)

def draw_heatmap(arrays_by_reference, sorted_chromosomes, all_refs):
    """Improved heatmap with better color scheme."""
    if arrays_by_reference.empty:
        print("No data available for heatmap. Skipping heatmap.")
        return

    summary_bp = arrays_by_reference.pivot_table(
        index="Chromosome",
        columns="Reference",
        values="ArraySize",
        aggfunc="sum",
        fill_value=0
    )

    summary_count = arrays_by_reference.pivot_table(
        index="Chromosome",
        columns="Reference",
        values="ArraySize",
        aggfunc="count",
        fill_value=0
    )

    summary_bp = summary_bp.reindex(index=sorted_chromosomes, columns=all_refs, fill_value=0)
    summary_count = summary_count.reindex(index=sorted_chromosomes, columns=all_refs, fill_value=0)

    summary_bp.to_csv(output_path("heatmap_total_array_bp_by_chromosome_reference.tsv"), sep="\t")
    summary_count.to_csv(output_path("heatmap_array_count_by_chromosome_reference.tsv"), sep="\t")

    values = np.log10(summary_bp.values.astype(float) + 1)

    fig_width = max(10, len(all_refs) * 0.55)
    fig_height = max(8, len(sorted_chromosomes) * 0.38)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Better color map for heatmap
    cmap = LinearSegmentedColormap.from_list('heat', ['#FEF0D9', '#FDD49E', '#FD8D3C', '#D7301F'])
    im = ax.imshow(values, aspect="auto", interpolation="nearest", cmap=cmap)

    ax.set_xticks(range(len(all_refs)))
    ax.set_xticklabels(all_refs, rotation=90, fontsize=9)
    ax.set_yticks(range(len(sorted_chromosomes)))
    ax.set_yticklabels(sorted_chromosomes, fontsize=9)

    ax.set_xlabel("SatDNA references", fontsize=12, weight='bold')
    ax.set_ylabel("Chromosome/scaffold", fontsize=12, weight='bold')
    ax.set_title("Total array abundance by chromosome and satDNA\n(log10 bp + 1)", fontsize=14, weight='bold')

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(total array bp + 1)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    fig.tight_layout()
    fig.savefig(output_path("heatmap_chromosome_vs_satdna_total_bp.png"), dpi=300, bbox_inches="tight", facecolor='white')
    fig.savefig(output_path("heatmap_chromosome_vs_satdna_total_bp.pdf"), bbox_inches="tight", facecolor='white')
    plt.close(fig)

def array_n50_l50(lengths):
    """Return N50 and L50 for a collection of array lengths."""
    vals = sorted((int(x) for x in lengths if pd.notna(x) and int(x) > 0), reverse=True)
    if not vals:
        return 0, 0
    half = sum(vals) / 2.0
    cumulative = 0
    for i, value in enumerate(vals, start=1):
        cumulative += value
        if cumulative >= half:
            return int(value), int(i)
    return int(vals[-1]), int(len(vals))

def build_nonoverlapping_satdna_composition(arrays_by_reference, pretty_to_length, sorted_chromosomes):
    """Convert possibly overlapping satDNA arrays into disjoint chromosome segments."""
    category_rows = []
    segment_rows = []

    for chrom in sorted_chromosomes:
        chrom_len = int(pretty_to_length[chrom])
        sub = arrays_by_reference[arrays_by_reference["Chromosome"] == chrom].copy()

        events = {}
        for _, row in sub.iterrows():
            start = max(1, int(row["Start"]))
            end = min(chrom_len, int(row["End"]))
            if end < start:
                continue
            ref = str(row["Reference"])
            events.setdefault(start, []).append((ref, 1))
            events.setdefault(end + 1, []).append((ref, -1))

        active_counts = {}
        category_bp = {}
        previous_pos = 1

        for pos in sorted(events):
            pos = min(max(1, int(pos)), chrom_len + 1)
            if pos > previous_pos:
                active_refs = sorted(
                    [ref for ref, count in active_counts.items() if count > 0],
                    key=natural_sort_key,
                )
                segment_len = pos - previous_pos
                if len(active_refs) == 1:
                    category = active_refs[0]
                elif len(active_refs) > 1:
                    category = "Multi_satDNA_overlap"
                else:
                    category = "Non_satDNA"

                category_bp[category] = category_bp.get(category, 0) + segment_len
                if category != "Non_satDNA":
                    segment_rows.append({
                        "Chromosome": chrom,
                        "Start": previous_pos,
                        "End": pos - 1,
                        "LengthBp": segment_len,
                        "Category": category,
                        "ActiveReferences": ",".join(active_refs),
                        "NumActiveReferences": len(active_refs),
                    })

            for ref, delta in events[pos]:
                active_counts[ref] = active_counts.get(ref, 0) + delta
                if active_counts[ref] <= 0:
                    active_counts.pop(ref, None)
            previous_pos = pos

        if previous_pos <= chrom_len:
            category_bp["Non_satDNA"] = category_bp.get("Non_satDNA", 0) + (chrom_len - previous_pos + 1)

        assigned = sum(category_bp.values())
        if assigned < chrom_len:
            category_bp["Non_satDNA"] = category_bp.get("Non_satDNA", 0) + (chrom_len - assigned)
        elif assigned > chrom_len:
            raise ValueError(f"Composition exceeds chromosome length for {chrom}: {assigned} > {chrom_len}")

        for category, bp in category_bp.items():
            category_rows.append({
                "Chromosome": chrom,
                "ChromosomeLengthBp": chrom_len,
                "Category": category,
                "CoveredBp": int(bp),
                "PercentChromosome": 100.0 * float(bp) / chrom_len if chrom_len else 0.0,
            })

    return pd.DataFrame(category_rows), pd.DataFrame(segment_rows)

def summarize_satdna_by_chromosome(arrays_by_reference, composition_long, pretty_to_length, sorted_chromosomes):
    rows = []

    for chrom in sorted_chromosomes:
        chrom_len = int(pretty_to_length[chrom])
        sub = arrays_by_reference[arrays_by_reference["Chromosome"] == chrom].copy()
        comp = composition_long[composition_long["Chromosome"] == chrom].copy()

        sat_comp = comp[comp["Category"] != "Non_satDNA"]
        union_bp = int(sat_comp["CoveredBp"].sum()) if not sat_comp.empty else 0
        overlap_bp = int(sat_comp.loc[sat_comp["Category"] == "Multi_satDNA_overlap", "CoveredBp"].sum()) if not sat_comp.empty else 0

        if sub.empty:
            rows.append({
                "Chromosome": chrom,
                "ChromosomeLengthBp": chrom_len,
                "NumSatDNAFamilies": 0,
                "NumArrays": 0,
                "RawSummedArrayBp": 0,
                "UnionSatDNABp": union_bp,
                "SatDNAPercentChromosome": 100.0 * union_bp / chrom_len if chrom_len else 0.0,
                "MultiSatDNAOverlapBp": overlap_bp,
                "MultiSatDNAOverlapPercent": 100.0 * overlap_bp / chrom_len if chrom_len else 0.0,
                "LargestArrayBp": 0,
                "LargestArraySatDNA": "NA",
                "LargestArrayStart": np.nan,
                "LargestArrayEnd": np.nan,
                "ArrayN50Bp": 0,
                "ArrayL50": 0,
                "MedianArrayBp": 0,
                "MeanArrayBp": 0,
                "ArrayP90Bp": 0,
                "ArraysPerMb": 0,
                "SatDNAFamiliesPerMb": 0,
                "TotalMonomerHits": 0,
                "WeightedMeanPident": np.nan,
                "DominantSatDNA": "NA",
                "DominantSatDNAExclusiveBp": 0,
                "DominantSatDNAPercentChromosome": 0.0,
            })
            continue

        n50, l50 = array_n50_l50(sub["ArraySize"].tolist())
        largest = sub.sort_values(["ArraySize", "Reference"], ascending=[False, True]).iloc[0]

        exclusive = comp[~comp["Category"].isin(["Non_satDNA", "Multi_satDNA_overlap"])].copy()
        if exclusive.empty:
            dominant_name = "NA"
            dominant_bp = 0
            dominant_pct = 0.0
        else:
            dominant = exclusive.sort_values(["CoveredBp", "Category"], ascending=[False, True]).iloc[0]
            dominant_name = str(dominant["Category"])
            dominant_bp = int(dominant["CoveredBp"])
            dominant_pct = float(dominant["PercentChromosome"])

        weights = sub["ArraySize"].astype(float)
        weighted_pident = float(np.average(sub["MeanPident"], weights=weights)) if weights.sum() > 0 else np.nan

        rows.append({
            "Chromosome": chrom,
            "ChromosomeLengthBp": chrom_len,
            "NumSatDNAFamilies": int(sub["Reference"].nunique()),
            "NumArrays": int(len(sub)),
            "RawSummedArrayBp": int(sub["ArraySize"].sum()),
            "UnionSatDNABp": union_bp,
            "SatDNAPercentChromosome": 100.0 * union_bp / chrom_len if chrom_len else 0.0,
            "MultiSatDNAOverlapBp": overlap_bp,
            "MultiSatDNAOverlapPercent": 100.0 * overlap_bp / chrom_len if chrom_len else 0.0,
            "LargestArrayBp": int(largest["ArraySize"]),
            "LargestArraySatDNA": str(largest["Reference"]),
            "LargestArrayStart": int(largest["Start"]),
            "LargestArrayEnd": int(largest["End"]),
            "ArrayN50Bp": n50,
            "ArrayL50": l50,
            "MedianArrayBp": float(sub["ArraySize"].median()),
            "MeanArrayBp": float(sub["ArraySize"].mean()),
            "ArrayP90Bp": float(sub["ArraySize"].quantile(0.90)),
            "ArraysPerMb": float(len(sub) / (chrom_len / 1e6)) if chrom_len else 0.0,
            "SatDNAFamiliesPerMb": float(sub["Reference"].nunique() / (chrom_len / 1e6)) if chrom_len else 0.0,
            "TotalMonomerHits": int(sub["NumMonomers"].sum()),
            "WeightedMeanPident": weighted_pident,
            "DominantSatDNA": dominant_name,
            "DominantSatDNAExclusiveBp": dominant_bp,
            "DominantSatDNAPercentChromosome": dominant_pct,
        })

    return pd.DataFrame(rows)

def summarize_satdna_by_chromosome_reference(arrays_by_reference, pretty_to_length, sorted_chromosomes, all_refs):
    rows = []
    for chrom in sorted_chromosomes:
        chrom_len = int(pretty_to_length[chrom])
        for ref in all_refs:
            sub = arrays_by_reference[
                (arrays_by_reference["Chromosome"] == chrom) &
                (arrays_by_reference["Reference"] == ref)
            ].copy()
            if sub.empty:
                continue
            n50, l50 = array_n50_l50(sub["ArraySize"].tolist())
            largest = sub.sort_values("ArraySize", ascending=False).iloc[0]
            rows.append({
                "Chromosome": chrom,
                "ChromosomeLengthBp": chrom_len,
                "Reference": ref,
                "NumArrays": int(len(sub)),
                "TotalArrayBpRaw": int(sub["ArraySize"].sum()),
                "RawPercentChromosome": 100.0 * float(sub["ArraySize"].sum()) / chrom_len if chrom_len else 0.0,
                "LargestArrayBp": int(largest["ArraySize"]),
                "LargestArrayStart": int(largest["Start"]),
                "LargestArrayEnd": int(largest["End"]),
                "ArrayN50Bp": n50,
                "ArrayL50": l50,
                "MedianArrayBp": float(sub["ArraySize"].median()),
                "MeanArrayBp": float(sub["ArraySize"].mean()),
                "ArrayP90Bp": float(sub["ArraySize"].quantile(0.90)),
                "TotalMonomerHits": int(sub["NumMonomers"].sum()),
                "MeanPidentWeightedByArrayBp": float(np.average(sub["MeanPident"], weights=sub["ArraySize"])) if sub["ArraySize"].sum() > 0 else np.nan,
            })
    return pd.DataFrame(rows)

def format_bp(value):
    """Format a base-pair value for compact plot labels."""
    value = float(value)
    if value >= 1e9:
        return f"{value / 1e9:.2f} Gb"
    if value >= 1e6:
        return f"{value / 1e6:.2f} Mb"
    if value >= 1e3:
        return f"{value / 1e3:.1f} kb"
    return f"{int(round(value))} bp"

def add_complete_genome_to_composition(composition_long, pretty_to_length, sorted_chromosomes):
    """Append a final whole-genome row to the non-overlapping composition table."""
    total_genome_bp = int(sum(int(pretty_to_length[c]) for c in sorted_chromosomes))
    genome_rows = (
        composition_long
        .groupby("Category", as_index=False)["CoveredBp"]
        .sum()
    )
    genome_rows["Chromosome"] = "Complete_genome"
    genome_rows["ChromosomeLengthBp"] = total_genome_bp
    genome_rows["PercentChromosome"] = np.where(
        total_genome_bp > 0,
        100.0 * genome_rows["CoveredBp"].astype(float) / total_genome_bp,
        0.0,
    )
    genome_rows = genome_rows[
        ["Chromosome", "ChromosomeLengthBp", "Category", "CoveredBp", "PercentChromosome"]
    ]
    return pd.concat([composition_long, genome_rows], ignore_index=True)

def build_satellitome_composition(arrays_by_reference, sorted_chromosomes, all_refs):
    """Summarize the raw summed array size of every satDNA family."""
    raw = (
        arrays_by_reference
        .groupby(["Chromosome", "Reference"], as_index=False)["ArraySize"]
        .sum()
        .rename(columns={"ArraySize": "TotalArrayBp"})
    )

    complete = (
        arrays_by_reference
        .groupby("Reference", as_index=False)["ArraySize"]
        .sum()
        .rename(columns={"ArraySize": "TotalArrayBp"})
    )
    complete["Chromosome"] = "Complete_genome"
    raw = pd.concat(
        [raw, complete[["Chromosome", "Reference", "TotalArrayBp"]]],
        ignore_index=True,
    )

    expected_rows = []
    for chrom in list(sorted_chromosomes) + ["Complete_genome"]:
        for ref in all_refs:
            expected_rows.append((chrom, ref))
    full = pd.DataFrame(expected_rows, columns=["Chromosome", "Reference"])
    raw = full.merge(raw, on=["Chromosome", "Reference"], how="left")
    raw["TotalArrayBp"] = raw["TotalArrayBp"].fillna(0).astype(int)

    totals = (
        raw.groupby("Chromosome", as_index=False)["TotalArrayBp"]
        .sum()
        .rename(columns={"TotalArrayBp": "TotalSatellitomeBpRaw"})
    )
    raw = raw.merge(totals, on="Chromosome", how="left")
    raw["PercentSatellitome"] = np.where(
        raw["TotalSatellitomeBpRaw"] > 0,
        100.0 * raw["TotalArrayBp"] / raw["TotalSatellitomeBpRaw"],
        0.0,
    )
    raw["RankWithinChromosome"] = (
        raw.groupby("Chromosome")["TotalArrayBp"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    return raw

def draw_satdna_percentage_composition(
    composition_long,
    satellitome_long,
    sorted_chromosomes,
    all_refs,
    color_map,
):
    """Draw a two-panel figure with improved aesthetics."""
    if composition_long.empty:
        print("No data available for satDNA percentage composition plot. Skipping.")
        return

    composition_with_genome = add_complete_genome_to_composition(
        composition_long, pretty_to_length, sorted_chromosomes
    )

    categories = list(all_refs)
    if (composition_with_genome["Category"] == "Multi_satDNA_overlap").any():
        categories.append("Multi_satDNA_overlap")
    categories.append("Non_satDNA")

    chromosome_order = list(sorted_chromosomes) + ["Complete_genome"]

    pivot = composition_with_genome.pivot_table(
        index="Chromosome",
        columns="Category",
        values="PercentChromosome",
        aggfunc="sum",
        fill_value=0,
    ).reindex(index=chromosome_order, columns=categories, fill_value=0)

    pivot.to_csv(output_path("satdna_percentage_composition_wide.tsv"), sep="\t")
    composition_with_genome.to_csv(
        output_path("satdna_percentage_composition_long_with_complete_genome.tsv"),
        sep="\t",
        index=False,
    )

    sat_wide_bp = satellitome_long.pivot_table(
        index="Chromosome",
        columns="Reference",
        values="TotalArrayBp",
        aggfunc="sum",
        fill_value=0,
    ).reindex(index=chromosome_order, columns=all_refs, fill_value=0)
    sat_wide_pct = satellitome_long.pivot_table(
        index="Chromosome",
        columns="Reference",
        values="PercentSatellitome",
        aggfunc="sum",
        fill_value=0,
    ).reindex(index=chromosome_order, columns=all_refs, fill_value=0)
    sat_wide_bp.to_csv(output_path("satellitome_total_array_bp_wide.tsv"), sep="\t")
    sat_wide_pct.to_csv(output_path("satellitome_percentage_composition_wide.tsv"), sep="\t")
    satellitome_long.to_csv(
        output_path("satellitome_composition_long.tsv"), sep="\t", index=False
    )

    # Improved bar plot with better aesthetics
    plot_chromosomes = ["Complete_genome"] + list(sorted_chromosomes)[::-1]
    fig_height = max(10, len(plot_chromosomes) * 0.45)
    fig, (ax_left, ax_right) = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(32, fig_height),
        sharey=True,
        gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.06},
    )
    y = np.arange(len(plot_chromosomes))

    category_colors = dict(color_map)
    category_colors["Multi_satDNA_overlap"] = (0.2, 0.2, 0.2)
    category_colors["Non_satDNA"] = (0.92, 0.92, 0.92)

    # LEFT PANEL: complete chromosome length = 100%.
    left_values = np.zeros(len(plot_chromosomes), dtype=float)
    for category in categories:
        values = pivot.loc[plot_chromosomes, category].to_numpy(dtype=float)
        ax_left.barh(
            y,
            values,
            left=left_values,
            height=0.7,
            label=category if category in ["Multi_satDNA_overlap", "Non_satDNA"] else "",
            color=category_colors.get(category, (0.5, 0.5, 0.5)),
            edgecolor='white',
            linewidth=0.2,
        )
        left_values += values

    total_sat_chrom_pct = (
        100.0 - pivot.loc[plot_chromosomes, "Non_satDNA"].to_numpy(dtype=float)
    )
    for yi, pct in zip(y, total_sat_chrom_pct):
        if pct >= 0.5:
            ax_left.text(
                min(pct + 0.5, 97), yi, f"{pct:.1f}%",
                va="center", ha="left", fontsize=8, weight='bold'
            )

    ax_left.set_xlim(0, 102)
    ax_left.set_xlabel("Percentage of chromosome length", fontsize=12, weight='bold')
    ax_left.set_ylabel("Chromosome/scaffold", fontsize=12, weight='bold')
    ax_left.set_title("satDNA as a percentage of chromosome length", fontsize=13, weight='bold')
    ax_left.set_yticks(y)
    ax_left.set_yticklabels(plot_chromosomes, fontsize=9)
    ax_left.grid(axis="x", linewidth=0.3, alpha=0.4)
    
    # Add color bar at bottom of left panel
    for spine in ["top", "right"]:
        ax_left.spines[spine].set_visible(False)

    # RIGHT PANEL: summed arrays in each chromosome = 100% of its satellitome.
    total_sat_bp_lookup = (
        satellitome_long[["Chromosome", "TotalSatellitomeBpRaw"]]
        .drop_duplicates("Chromosome")
        .set_index("Chromosome")["TotalSatellitomeBpRaw"]
        .to_dict()
    )

    for yi, chrom in zip(y, plot_chromosomes):
        row = satellitome_long[satellitome_long["Chromosome"] == chrom].copy()
        row = row[row["TotalArrayBp"] > 0].sort_values(
            ["TotalArrayBp", "Reference"], ascending=[False, True]
        )
        current_left = 0.0
        for _, item in row.iterrows():
            ref = str(item["Reference"])
            pct = float(item["PercentSatellitome"])
            ax_right.barh(
                yi,
                pct,
                left=current_left,
                height=0.7,
                color=color_map.get(ref, (0.5, 0.5, 0.5)),
                edgecolor='white',
                linewidth=0.2,
            )
            current_left += pct

        total_bp = int(total_sat_bp_lookup.get(chrom, 0))
        ax_right.text(
            100.5,
            yi,
            format_bp(total_bp),
            va="center",
            ha="left",
            fontsize=8,
            clip_on=False,
            weight='bold'
        )

    ax_right.set_xlim(0, 104)
    ax_right.set_xlabel("Percentage of the satellitome (summed arrays = 100%)", fontsize=12, weight='bold')
    ax_right.set_title("Relative satDNA composition of each satellitome", fontsize=13, weight='bold')
    ax_right.grid(axis="x", linewidth=0.3, alpha=0.4)
    ax_right.tick_params(axis="y", labelleft=False)

    for spine in ["top", "right"]:
        ax_right.spines[spine].set_visible(False)

    # Improved legend
    legend_categories = []
    legend_handles = []
    
    # Add satDNA references with colors
    for ref in all_refs:
        if ref in color_map:
            legend_categories.append(ref)
            legend_handles.append(
                Rectangle((0, 0), 1, 1, facecolor=color_map[ref], edgecolor='white', linewidth=0.2)
            )
    
    # Add overlap and non-satDNA
    legend_categories.append("Multi_satDNA_overlap")
    legend_handles.append(
        Rectangle((0, 0), 1, 1, facecolor=category_colors["Multi_satDNA_overlap"], edgecolor='white', linewidth=0.2)
    )
    legend_categories.append("Non_satDNA")
    legend_handles.append(
        Rectangle((0, 0), 1, 1, facecolor=category_colors["Non_satDNA"], edgecolor='#CCCCCC', linewidth=0.2)
    )

    ax_right.legend(
        legend_handles,
        legend_categories,
        title="SatDNA / category",
        bbox_to_anchor=(1.12, 1),
        loc="upper left",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        framealpha=0.95,
        edgecolor='#888888'
    )

    fig.suptitle(
        "Chromosomal satDNA abundance and satellitome composition",
        fontsize=16,
        y=0.995,
        weight='bold'
    )
    fig.tight_layout(rect=[0, 0, 0.89, 0.985])
    fig.savefig(
        output_path("satdna_percentage_composition_by_chromosome.png"),
        dpi=600,
        bbox_inches="tight",
        facecolor='white'
    )
    fig.savefig(
        output_path("satdna_percentage_composition_by_chromosome.pdf"),
        bbox_inches="tight",
        facecolor='white'
    )
    plt.close(fig)

def save_metrics_workbook(chromosome_metrics, chromosome_reference_metrics, composition_long, composition_segments, satellitome_long):
    out_xlsx = output_path("satdna_metrics_by_chromosome.xlsx")
    try:
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
            chromosome_metrics.to_excel(writer, sheet_name="chromosome_metrics", index=False)
            chromosome_reference_metrics.to_excel(writer, sheet_name="chromosome_x_satDNA", index=False)
            composition_long.to_excel(writer, sheet_name="composition_long", index=False)
            composition_segments.to_excel(writer, sheet_name="nonoverlap_segments", index=False)
            satellitome_long.to_excel(writer, sheet_name="satellitome_composition", index=False)
            arrays_by_reference.to_excel(writer, sheet_name="all_arrays", index=False)
        print("Saved Excel workbook:", out_xlsx)
    except ImportError:
        print("openpyxl is not installed; TSV tables were saved, but the XLSX workbook was skipped.")

def make_top_arrays(arrays_by_reference, top_n):
    if arrays_by_reference.empty:
        return arrays_by_reference.copy()
    return (
        arrays_by_reference
        .sort_values(["Reference","ArraySize"], ascending=[True, False])
        .groupby("Reference", group_keys=False)
        .head(top_n)
        .sort_values(["Chromosome","Reference","Start","End"])
        .reset_index(drop=True)
    )

# Main execution
accession_to_pretty, pretty_to_length, accession_to_length, records_info = build_header_mapping(FASTA_FILE)
save_header_mapping(records_info, output_path("sequence_name_mapping.tsv"))

raw_hits = read_raw_hits(RAW_HITS_FILE)

if raw_hits.empty:
    print("No BLAST hits found. Skipping plots.")
    raise SystemExit(0)

raw_hits["Chromosome_original_accession"] = raw_hits["Chromosome"]
raw_hits["Chromosome"] = raw_hits["Chromosome"].map(lambda x: accession_to_pretty.get(x, x))

valid_names = set(pretty_to_length.keys())
raw_hits_mapped = raw_hits[raw_hits["Chromosome"].isin(valid_names)].copy()

if raw_hits_mapped.empty or len(pretty_to_length) == 0:
    print("No relevant chromosomes/scaffolds after standardization. Skipping plots.")
    raise SystemExit(0)

raw_hits_mapped = raw_hits_mapped.sort_values(["Chromosome","Reference","Start","End"]).reset_index(drop=True)

raw_hits_mapped.to_csv(output_path("raw_blast_hits_mapped.tsv"), sep="\t", index=False)
make_valid_monomers_bed(raw_hits_mapped, output_path("valid_monomers.bed"))

arrays_by_reference, extension_factors = create_arrays_by_reference(raw_hits_mapped)

extension_factors.to_csv(output_path("adaptive_merge_distance_by_reference.tsv"), sep="\t", index=False)

if arrays_by_reference.empty:
    print("No arrays passed the minimum monomer filter. Skipping plots.")
    raise SystemExit(0)

merged_regions = create_merged_regions_multi_satdna(arrays_by_reference, overlap_dist=0)

arrays_by_reference.to_csv(output_path("arrays_by_reference.tsv"), sep="\t", index=False)
arrays_by_reference.to_csv(output_path("arrays_by_reference_for_plot.tsv"), sep="\t", index=False)
merged_regions.to_csv(output_path("merged_regions_multi_satdna.tsv"), sep="\t", index=False)

# Backward-compatible output name
merged_regions.to_csv(output_path("merged_arrays.tsv"), sep="\t", index=False)

all_refs = sorted(arrays_by_reference["Reference"].dropna().unique(), key=natural_sort_key)
sorted_chromosomes = sorted(pretty_to_length.keys(), key=natural_sort_key)
highlight_chromosomes = parse_highlight_chromosomes(HIGHLIGHT_CHROMOSOMES_RAW, sorted_chromosomes)

color_map = get_color_map(all_refs)
with open(os.path.join(GENOME_NAME, "reference_colors.json"), "w") as f:
    json.dump({k: list(v) for k, v in color_map.items()}, f, indent=2)

# Chromosome-level metrics and non-overlapping percentage composition
composition_long, composition_segments = build_nonoverlapping_satdna_composition(
    arrays_by_reference, pretty_to_length, sorted_chromosomes
)
chromosome_metrics = summarize_satdna_by_chromosome(
    arrays_by_reference, composition_long, pretty_to_length, sorted_chromosomes
)
chromosome_reference_metrics = summarize_satdna_by_chromosome_reference(
    arrays_by_reference, pretty_to_length, sorted_chromosomes, all_refs
)

composition_long.to_csv(output_path("satdna_percentage_composition_long.tsv"), sep="\t", index=False)
composition_segments.to_csv(output_path("satdna_nonoverlapping_segments.tsv"), sep="\t", index=False)
chromosome_metrics.to_csv(output_path("satdna_metrics_by_chromosome.tsv"), sep="\t", index=False)
chromosome_reference_metrics.to_csv(output_path("satdna_metrics_by_chromosome_reference.tsv"), sep="\t", index=False)

satellitome_long = build_satellitome_composition(
    arrays_by_reference, sorted_chromosomes, all_refs
)

draw_satdna_percentage_composition(
    composition_long, satellitome_long, sorted_chromosomes, all_refs, color_map
)
save_metrics_workbook(
    chromosome_metrics, chromosome_reference_metrics, composition_long,
    composition_segments, satellitome_long
)

top_arrays = make_top_arrays(arrays_by_reference, TOP_N_ARRAYS)
top_arrays.to_csv(output_path(f"top_{TOP_N_ARRAYS}_arrays_by_reference.tsv"), sep="\t", index=False)

draw_chromosome_array_plot(
    plot_arrays=arrays_by_reference,
    pretty_to_length=pretty_to_length,
    sorted_chromosomes=sorted_chromosomes,
    all_refs=all_refs,
    color_map=color_map,
    out_prefix="chromosomes_with_annotations_exact_scale",
    enhanced_visibility=False,
    min_visible_kb=80,
    title_extra="",
    top_arrays=None,
    combine_with_top=False,
    highlight_chromosomes=highlight_chromosomes
)

draw_chromosome_array_plot(
    plot_arrays=arrays_by_reference,
    pretty_to_length=pretty_to_length,
    sorted_chromosomes=sorted_chromosomes,
    all_refs=all_refs,
    color_map=color_map,
    out_prefix="chromosomes_with_annotations",
    enhanced_visibility=True,
    min_visible_kb=80,
    title_extra="",
    top_arrays=None,
    combine_with_top=False,
    highlight_chromosomes=highlight_chromosomes
)

draw_chromosome_array_plot(
    plot_arrays=arrays_by_reference,
    pretty_to_length=pretty_to_length,
    sorted_chromosomes=sorted_chromosomes,
    all_refs=all_refs,
    color_map=color_map,
    out_prefix=f"chromosomes_with_annotations_plus_top_{TOP_N_ARRAYS}",
    enhanced_visibility=True,
    min_visible_kb=80,
    title_extra="",
    top_arrays=top_arrays,
    combine_with_top=True,
    highlight_chromosomes=highlight_chromosomes
)

draw_scatter(arrays_by_reference, sorted_chromosomes, all_refs, color_map)
draw_heatmap(arrays_by_reference, sorted_chromosomes, all_refs)

summary = (
    arrays_by_reference
    .groupby("Reference")
    .agg(
        NumArrays=("ArraySize", "count"),
        TotalArrayBp=("ArraySize", "sum"),
        MinArrayBp=("ArraySize", "min"),
        MedianArrayBp=("ArraySize", "median"),
        MaxArrayBp=("ArraySize", "max"),
        TotalMonomerHits=("NumMonomers", "sum"),
        MedianMonomersPerArray=("NumMonomers", "median"),
        MeanPident=("MeanPident", "mean"),
        OrientationConsistency=("OrientationConsistency", "mean")
    )
    .reset_index()
)
summary.to_csv(output_path("summary_by_reference.tsv"), sep="\t", index=False)

print(" Analysis completed successfully!")
print(f" Output directory: {GENOME_NAME}")
print(f" Output prefix: {OUTPUT_PREFIX}")
print("")
print(" Generated plots:")
print(f"  - {output_path('chromosomes_with_annotations.png')}")
print(f"  - {output_path('chromosomes_with_annotations_exact_scale.png')}")
print(f"  - {output_path('chromosomes_with_annotations_plus_top_' + str(TOP_N_ARRAYS) + '.png')}")
print(f"  - {output_path('array_chromosome_vs_size_scatter.png')}")
print(f"  - {output_path('heatmap_chromosome_vs_satdna_total_bp.png')}")
print(f"  - {output_path('satdna_percentage_composition_by_chromosome.png')}")
print("")
print(" Generated tables:")
print(f"  - {output_path('raw_blast_hits_mapped.tsv')}")
print(f"  - {output_path('valid_monomers.bed')}")
print(f"  - {output_path('arrays_by_reference.tsv')}")
print(f"  - {output_path('summary_by_reference.tsv')}")
print(f"  - {output_path('satdna_metrics_by_chromosome.xlsx')}")
EOF

done
