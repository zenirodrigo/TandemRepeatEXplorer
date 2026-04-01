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
    """
    Normaliza IDs de sequência para um formato comparável entre:
    - FASTA headers
    - record.id do Biopython
    - sseqid do BLAST

    Casos tratados:
    - CM038702.1
    - gb|CM038702.1|
    - ref|NC_000001.1|
    - lcl|scaffold_42
    - gi|12345|gb|CM038702.1|
    - strings com descrição após o accession
    """
    if acc is None:
        return ""

    acc = str(acc).strip()
    if not acc:
        return ""

    # pega só o primeiro token antes de descrições
    first = acc.split()[0]

    # remove > caso venha header bruto
    first = first.lstrip(">")

    # caso simples já limpo
    if "|" not in first:
        return first

    parts = [p for p in first.split("|") if p != ""]
    if not parts:
        return first

    db_tags = {
        "gb", "ref", "emb", "dbj", "sp", "tr", "lcl", "gi",
        "gnl", "tpg", "tpe", "tpd", "pdb"
    }

    # prioridade: procurar algo com cara de accession real
    # ex.: CM038702.1, NC_000001.1, NW_..., scaffold_42
    accession_like = []
    for p in parts:
        if p.lower() in db_tags:
            continue
        if re.search(r'[A-Za-z]', p):
            accession_like.append(p)

    if accession_like:
        # em gi|123|gb|CM038702.1| queremos o último accession útil
        return accession_like[-1]

    # fallback
    return parts[-1]

def infer_pretty_name_from_header(header, fallback_index):
    h = str(header).strip()
    h_low = h.lower()

    m_chr = re.search(r'chromosome\s+([0-9]+)', h_low)
    if m_chr:
        return f"Chromosome{m_chr.group(1)}"

    m_chr2 = re.search(r'\bchr(?:omosome)?[_\s-]*([0-9]+)\b', h_low)
    if m_chr2:
        return f"Chromosome{m_chr2.group(1)}"

    m_scaf = re.search(r'scaffold[_\s-]*([0-9]+)', h_low)
    if m_scaf:
        return f"scaffold{m_scaf.group(1)}"

    m_ctg = re.search(r'contig[_\s-]*([0-9]+)', h_low)
    if m_ctg:
        return f"contig{m_ctg.group(1)}"

    # tenta usar accession se tiver algo informativo
    token = sanitize_accession(h)
    if token:
        return token

    return f"Sequence{fallback_index}"

def build_header_mapping(fasta_path):
    records_info = []
    accession_to_pretty = OrderedDict()
    pretty_to_length = OrderedDict()

    idx = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        idx += 1
        raw_header = record.description.strip()

        candidates = []
        for c in [record.id, raw_header]:
            val = sanitize_accession(c)
            if val and val not in candidates:
                candidates.append(val)

        pretty = infer_pretty_name_from_header(raw_header, idx)

        if pretty in pretty_to_length:
            pretty = f"{pretty}_{idx}"

        for acc in candidates:
            accession_to_pretty[acc] = pretty

        pretty_to_length[pretty] = len(record.seq)

        records_info.append({
            "index": idx,
            "record_id": record.id,
            "raw_header": raw_header,
            "normalized_candidates": ";".join(candidates),
            "pretty_name": pretty,
            "length": len(record.seq)
        })

    return accession_to_pretty, pretty_to_length, records_info

def read_bed(filepath):
    if (not os.path.exists(filepath)) or os.path.getsize(filepath) == 0:
        return pd.DataFrame(columns=["Chromosome","Start","End","Reference","References"])

    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "Reference"]
    )

    df["Chromosome_raw"] = df["Chromosome"].astype(str)
    df["Chromosome"] = df["Chromosome_raw"].map(sanitize_accession)
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
                merged.append([chrom, current_start, current_end, sorted(refs, key=natural_sort_key)])
                current_start, current_end = s, e
                refs = set(rlist)

        if current_start is not None:
            merged.append([chrom, current_start, current_end, sorted(refs, key=natural_sort_key)])

    return pd.DataFrame(merged, columns=["Chromosome", "Start", "End", "References"])

def save_header_mapping(records_info, out_tsv):
    pd.DataFrame(records_info).to_csv(out_tsv, sep="\t", index=False)

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
    print("Example normalized BED names:", df["Chromosome_original_accession"].head(10).tolist())
    print("Example FASTA normalized names:", list(accession_to_pretty.keys())[:10])
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

if os.path.exists(color_file):
    with open(color_file, "r") as f:
        color_map = json.load(f)
else:
    color_map = {}

missing_refs = [r for r in all_refs if r not in color_map]
if missing_refs:
    new_palette = sns.color_palette("pastel", len(missing_refs))
    for i, ref in enumerate(missing_refs):
        rgb = new_palette[i]
        color_map[ref] = [float(rgb[0]), float(rgb[1]), float(rgb[2])]

with open(color_file, "w") as f:
    json.dump(color_map, f, indent=2)

color_map = {k: tuple(v) for k, v in color_map.items()}

df_merged.to_csv(os.path.join(genome_name, "merged_arrays.tsv"), sep="\t", index=False)
df_filtered.to_csv(os.path.join(genome_name, "valid_monomers_mapped.tsv"), sep="\t", index=False)

sorted_chromosomes = sorted(pretty_to_length.keys(), key=natural_sort_key)

reference_sizes = {}
for ref in all_refs:
    subset = df_merged[df_merged["References"].apply(lambda x: ref in x)]
    sizes_ = subset["End"] - subset["Start"] + 1
    reference_sizes[ref] = (
        int(sizes_.min()) if not sizes_.empty else 0,
        int(sizes_.max()) if not sizes_.empty else 0
    )

plt.figure(figsize=(36, 24))
spacing = 1.5

for idx, chrom in enumerate(sorted_chromosomes):
    length = pretty_to_length[chrom]
    length_mb = length / 1e6
    y_pos = idx * spacing

    plt.plot([0, length_mb], [y_pos, y_pos], color="lightgray", linewidth=3, zorder=1)
    plt.fill_between(
        [0, length_mb],
        y_pos - 0.4,
        y_pos + 0.4,
        color="lightgray",
        alpha=0.6,
        zorder=1
    )

    chrom_data = df_merged[df_merged["Chromosome"] == chrom]
    for _, row in chrom_data.iterrows():
        start_mb = row["Start"] / 1e6
        end_mb = row["End"] / 1e6
        refs = row["References"]

        if not refs:
            continue

        segment_width = (end_mb - start_mb) / len(refs) if len(refs) > 0 else 0
        for i, ref in enumerate(refs):
            c = color_map.get(ref, "#cccccc")
            plt.fill_between(
                [start_mb + i * segment_width, start_mb + (i + 1) * segment_width],
                y_pos - 0.4,
                y_pos + 0.4,
                color=c,
                alpha=0.85,
                zorder=2
            )

handles = [plt.Line2D([0], [0], color=color_map[ref], lw=5) for ref in all_refs]
labels = [
    f"{ref} (min: {reference_sizes[ref][0]:,} bp, max: {reference_sizes[ref][1]:,} bp)"
    for ref in all_refs
]

plt.legend(handles, labels, title="References", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=12)
plt.yticks([i * spacing for i in range(len(sorted_chromosomes))], sorted_chromosomes, fontsize=12)
plt.xlabel("Base pairs (Mb)", fontsize=14)
plt.ylabel("Sequences", fontsize=14)
plt.title("Distribution of tandem repeat arrays across chromosomes/scaffolds", fontsize=16)
plt.xlim(0, max(pretty_to_length.values()) / 1e6)
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(genome_name, "chromosomes_with_annotations.png"), dpi=300, bbox_inches="tight")
plt.close()

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
