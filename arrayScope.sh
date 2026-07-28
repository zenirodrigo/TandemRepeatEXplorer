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
        -qcov_hsp_perc 70 \
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

read -e -p "Enter genome file names (space-separated): " input_biblios
read -p "How many chromosome/scaffold sequences will be used? " num_sequences
read -e -p "Enter reference (satDNA or another tandem repeat MONOMER) files (space-separated): " refs_in
read -p "How many monomers will be used to create an array (minimum monomers in this study)? " multiplier
read -p "Maximum gap allowed for adaptive merging, in bp [default: 10000]: " adaptive_gap_cap
read -p "Fallback gap if adaptive distance cannot be estimated, in bp [default: 2000]: " fallback_gap
read -p "Minimum BLAST percent identity to keep a hit [default: 0 = no extra filter]: " min_pident
read -p "How many largest arrays per satDNA should be highlighted in the TOP plot? [default: 2]: " top_n_arrays
read -p "Comma-separated chromosome/scaffold names to highlight, optional (example: ChrB,B,microB): " highlight_chromosomes
read -p "How many threads do you want to use? (e.g., 4, 8, etc.): " NUM_THREADS

adaptive_gap_cap="${adaptive_gap_cap:-10000}"
fallback_gap="${fallback_gap:-2000}"
min_pident="${min_pident:-0}"
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

    parallel --jobs "$NUM_THREADS" \
        run_blast_for_ref {} "$temp_genome" "$multiplier" "$tmp_parallel_dir/{}.raw_hits.tsv" 1 \
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
    } > "$genome_name/raw_blast_hits.tsv"

    if [ ! -s "$genome_name/raw_blast_hits.tsv" ]; then
        echo "The file raw_blast_hits.tsv is empty for $genome_name."
    else
        echo "The file raw_blast_hits.tsv was created for $genome_name."
    fi

python3 - <<EOF
import os
import json
import re
import math
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap

fasta_file = "$temp_genome"
raw_hits_file = "$genome_name/raw_blast_hits.tsv"
genome_name = "$genome_name"

MIN_MONOMERS_PER_ARRAY = int("$multiplier")
ADAPTIVE_GAP_CAP = int("$adaptive_gap_cap")
FALLBACK_GAP = int("$fallback_gap")
MIN_PIDENT = float("$min_pident")
TOP_N_ARRAYS = int("$top_n_arrays")
HIGHLIGHT_CHROMOSOMES_RAW = "$highlight_chromosomes"

os.makedirs(genome_name, exist_ok=True)

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

    # Preserve simple sequence names already present in the FASTA header,
    # such as B1, B2 or microB.
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
            "StrandMix": strand_mix
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
    cmap_names = ["tab20", "tab20b", "tab20c", "Dark2", "Set1"]
    colors = []
    for cmap_name in cmap_names:
        cmap = plt.get_cmap(cmap_name)
        n = cmap.N if hasattr(cmap, "N") else 20
        for i in range(n):
            rgba = cmap(i)
            colors.append((float(rgba[0]), float(rgba[1]), float(rgba[2])))

    color_map = {}
    for i, ref in enumerate(refs):
        color_map[ref] = colors[i % len(colors)]
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
    highlight_chromosomes = highlight_chromosomes or set()

    sorted_chromosomes_for_axis = list(sorted_chromosomes)[::-1]
    max_len_mb = max(pretty_to_length.values()) / 1e6

    chrom_height = 0.58
    array_height = 0.50
    spacing = 1.08
    fig_width = 38
    fig_height = max(14, len(sorted_chromosomes_for_axis) * 0.50)

    if combine_with_top:
        fig, axes = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(fig_width, fig_height * 1.55),
            sharex=True,
            gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.18}
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
                facecolor = "#fff3cd"
                edgecolor = "#7a5c00"
                linewidth = 1.15
            else:
                facecolor = "#f1f1f1"
                edgecolor = "#bdbdbd"
                linewidth = 0.40

            body = FancyBboxPatch(
                (0, y - chrom_height / 2),
                length_mb,
                chrom_height,
                boxstyle=f"round,pad=0.00,rounding_size={chrom_height * 0.22}",
                linewidth=linewidth,
                edgecolor=edgecolor,
                facecolor=facecolor,
                alpha=1.00,
                zorder=1,
            )
            ax.add_patch(body)

        dataset_sorted = dataset.sort_values("ArraySize", ascending=True).reset_index(drop=True)

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

        ax.set_yticks([y_positions[c] for c in sorted_chromosomes_for_axis])
        ax.set_yticklabels(sorted_chromosomes_for_axis, fontsize=10)
        ax.set_ylabel("Chromosomes/scaffolds", fontsize=12)

        if enhanced_visibility:
            title_suffix = f"visible mode; arrays < {min_visible_kb} kb widened for display"
        else:
            title_suffix = "exact scale"

        ax.set_title(f"{panel_title} ({title_suffix}){title_extra}", fontsize=14, pad=12)
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

    ax_list[-1].set_xlabel("Position on chromosome/scaffold (Mb)", fontsize=12)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=color_map[ref], edgecolor=color_map[ref], alpha=0.98)
        for ref in all_refs
    ]
    labels = [str(ref) for ref in all_refs]

    ax_list[0].legend(
        handles,
        labels,
        title="References",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
        title_fontsize=10,
        frameon=True,
        borderaxespad=0.0,
        ncol=1
    )

    fig.tight_layout()
    fig.savefig(os.path.join(genome_name, f"{out_prefix}.png"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(genome_name, f"{out_prefix}.pdf"), bbox_inches="tight")
    plt.close(fig)

def draw_scatter(arrays_by_reference, sorted_chromosomes, all_refs, color_map):
    df = arrays_by_reference.copy()
    if df.empty:
        print("No data available for the scatter plot. Skipping scatter.")
        return

    df["Chromosome"] = pd.Categorical(df["Chromosome"], categories=sorted_chromosomes, ordered=True)
    df = df.sort_values(["Chromosome","Reference","ArraySize"]).reset_index(drop=True)

    fig_width = max(18, len(sorted_chromosomes) * 0.45)
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    ref_offsets = {}
    if len(all_refs) == 1:
        ref_offsets[all_refs[0]] = 0
    else:
        offsets = np.linspace(-0.34, 0.34, len(all_refs))
        ref_offsets = {ref: offsets[i] for i, ref in enumerate(all_refs)}

    chrom_to_x = {chrom: i for i, chrom in enumerate(sorted_chromosomes)}

    for ref in all_refs:
        sub = df[df["Reference"] == ref].copy()
        if sub.empty:
            continue
        xs = [chrom_to_x[c] + ref_offsets[ref] for c in sub["Chromosome"]]
        ax.scatter(
            xs,
            sub["ArraySize"],
            s=34,
            alpha=0.82,
            color=color_map.get(ref, "#777777"),
            label=ref,
            edgecolors="none"
        )

    ax.set_xticks(range(len(sorted_chromosomes)))
    ax.set_xticklabels(sorted_chromosomes, rotation=90)
    ax.set_xlabel("Sequences")
    ax.set_ylabel("Array size (bp)")
    ax.set_title("Array size distribution by chromosome/scaffold")
    ax.grid(axis="y", which="major", color="#d0d0d0", linewidth=0.45, alpha=0.65)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.legend(title="References", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9, title_fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(genome_name, "array_chromosome_vs_size_scatter.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(genome_name, "array_chromosome_vs_size_scatter.pdf"), bbox_inches="tight")
    plt.close(fig)

def draw_heatmap(arrays_by_reference, sorted_chromosomes, all_refs):
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

    summary_bp.to_csv(os.path.join(genome_name, "heatmap_total_array_bp_by_chromosome_reference.tsv"), sep="\t")
    summary_count.to_csv(os.path.join(genome_name, "heatmap_array_count_by_chromosome_reference.tsv"), sep="\t")

    values = np.log10(summary_bp.values.astype(float) + 1)

    fig_width = max(10, len(all_refs) * 0.55)
    fig_height = max(8, len(sorted_chromosomes) * 0.34)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(values, aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(all_refs)))
    ax.set_xticklabels(all_refs, rotation=90)
    ax.set_yticks(range(len(sorted_chromosomes)))
    ax.set_yticklabels(sorted_chromosomes)

    ax.set_xlabel("satDNA reference")
    ax.set_ylabel("Chromosome/scaffold")
    ax.set_title("Total array abundance by chromosome and satDNA (log10 bp + 1)")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10(total array bp + 1)")

    fig.tight_layout()
    fig.savefig(os.path.join(genome_name, "heatmap_chromosome_vs_satdna_total_bp.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(genome_name, "heatmap_chromosome_vs_satdna_total_bp.pdf"), bbox_inches="tight")
    plt.close(fig)

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

accession_to_pretty, pretty_to_length, accession_to_length, records_info = build_header_mapping(fasta_file)
save_header_mapping(records_info, os.path.join(genome_name, "sequence_name_mapping.tsv"))

raw_hits = read_raw_hits(raw_hits_file)

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

raw_hits_mapped.to_csv(os.path.join(genome_name, "raw_blast_hits_mapped.tsv"), sep="\t", index=False)
make_valid_monomers_bed(raw_hits_mapped, os.path.join(genome_name, "valid_monomers.bed"))

arrays_by_reference, extension_factors = create_arrays_by_reference(raw_hits_mapped)

extension_factors.to_csv(os.path.join(genome_name, "adaptive_merge_distance_by_reference.tsv"), sep="\t", index=False)

if arrays_by_reference.empty:
    print("No arrays passed the minimum monomer filter. Skipping plots.")
    raise SystemExit(0)

merged_regions = create_merged_regions_multi_satdna(arrays_by_reference, overlap_dist=0)

arrays_by_reference.to_csv(os.path.join(genome_name, "arrays_by_reference.tsv"), sep="\t", index=False)
arrays_by_reference.to_csv(os.path.join(genome_name, "arrays_by_reference_for_plot.tsv"), sep="\t", index=False)
merged_regions.to_csv(os.path.join(genome_name, "merged_regions_multi_satdna.tsv"), sep="\t", index=False)

# Backward-compatible output name from the older script.
merged_regions.to_csv(os.path.join(genome_name, "merged_arrays.tsv"), sep="\t", index=False)

all_refs = sorted(arrays_by_reference["Reference"].dropna().unique(), key=natural_sort_key)
sorted_chromosomes = sorted(pretty_to_length.keys(), key=natural_sort_key)
highlight_chromosomes = parse_highlight_chromosomes(HIGHLIGHT_CHROMOSOMES_RAW, sorted_chromosomes)

color_map = get_color_map(all_refs)
with open(os.path.join(genome_name, "reference_colors.json"), "w") as f:
    json.dump({k: list(v) for k, v in color_map.items()}, f, indent=2)

top_arrays = make_top_arrays(arrays_by_reference, TOP_N_ARRAYS)
top_arrays.to_csv(os.path.join(genome_name, f"top_{TOP_N_ARRAYS}_arrays_by_reference.tsv"), sep="\t", index=False)

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
        MeanPident=("MeanPident", "mean")
    )
    .reset_index()
)
summary.to_csv(os.path.join(genome_name, "summary_by_reference.tsv"), sep="\t", index=False)

print("Plots and tables successfully generated in:", genome_name)
print("Saved image:", os.path.join(genome_name, "chromosomes_with_annotations.png"))
print("Saved image:", os.path.join(genome_name, f"chromosomes_with_annotations_plus_top_{TOP_N_ARRAYS}.png"))
print("Saved image:", os.path.join(genome_name, "array_chromosome_vs_size_scatter.png"))
print("Saved image:", os.path.join(genome_name, "heatmap_chromosome_vs_satdna_total_bp.png"))
print("Saved table:", os.path.join(genome_name, "raw_blast_hits.tsv"))
print("Saved table:", os.path.join(genome_name, "raw_blast_hits_mapped.tsv"))
print("Saved table:", os.path.join(genome_name, "valid_monomers.bed"))
print("Saved table:", os.path.join(genome_name, "arrays_by_reference.tsv"))
print("Saved table:", os.path.join(genome_name, "merged_regions_multi_satdna.tsv"))
print("Saved table:", os.path.join(genome_name, "adaptive_merge_distance_by_reference.tsv"))
print("Saved table:", os.path.join(genome_name, "summary_by_reference.tsv"))
print("Saved table:", os.path.join(genome_name, "sequence_name_mapping.tsv"))
EOF

done
