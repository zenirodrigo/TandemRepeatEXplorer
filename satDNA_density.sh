#!/bin/bash
set -euo pipefail

# -------------------------
# Helpers
# -------------------------
remove_extensions() {
    local filename="$1"
    filename="$(basename "$filename")"
    filename="${filename%.fasta}"
    filename="${filename%.fa}"
    filename="${filename%.fna}"
    echo "$filename"
}

safe_id_from_path() {
    local p="$1"
    echo "$p" | sed 's/[^A-Za-z0-9]/_/g'
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found in PATH: $1"
}

# Multiply every FASTA record sequence by "multiplier"
multiply_fasta_python() {
    local in_fa="$1"
    local out_fa="$2"
    local multiplier="$3"

    python3 - "$in_fa" "$out_fa" "$multiplier" <<'PY'
import sys

in_fa  = sys.argv[1]
out_fa = sys.argv[2]
mult   = int(sys.argv[3])

def fasta_iter(path):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)

with open(out_fa, "w", encoding="utf-8") as out:
    for header, seq in fasta_iter(in_fa):
        out.write(header + "\n")
        seq2 = seq * mult
        for i in range(0, len(seq2), 80):
            out.write(seq2[i:i+80] + "\n")
PY
}

run_blast_for_ref_to_bed() {
    # Writes one BED file per reference to avoid parallel race conditions
    local ref_fasta="$1"
    local temp_genome_db="$2"
    local multiplier="$3"
    local out_bed="$4"
    local blast_threads="$5"
    local tmpdir="$6"

    if [ ! -f "$ref_fasta" ]; then
        echo "WARNING: Reference file not found: $ref_fasta" >&2
        : > "$out_bed"
        return 0
    fi

    local ref_label
    ref_label="$(remove_extensions "$ref_fasta")"

    local multiplied_fasta="$tmpdir/multiplied_reference_${ref_label}.fasta"
    multiply_fasta_python "$ref_fasta" "$multiplied_fasta" "$multiplier"

    if [ ! -s "$multiplied_fasta" ]; then
        echo "WARNING: Multiplied FASTA is empty for: $ref_fasta" >&2
        : > "$out_bed"
        return 0
    fi

    local blast_output="$tmpdir/blast_${ref_label}.out"

    blastn \
      -task blastn \
      -outfmt "6" \
      -db "$temp_genome_db" \
      -query "$multiplied_fasta" \
      -out "$blast_output" \
      -evalue 1e-10 \
      -qcov_hsp_perc 70 \
      -num_threads "$blast_threads" || true

    if [ ! -s "$blast_output" ]; then
        : > "$out_bed"
        return 0
    fi

    # BED4: Chrom Start End Reference (Reference = qseqid)
    # Merge hits within 2000 bp ONLY if same chrom AND same reference
    awk '{
        start=($9 < $10) ? $9 : $10;
        end=($9 < $10) ? $10 : $9;
        print $2, start, end, $1
    }' "$blast_output" \
    | sort -k1,1 -k2,2n -k4,4 \
    | awk -v OFS='\t' -v dist=2000 '
        {
            if (NR == 1) {
                chr=$1; start=$2; end=$3; ref=$4
            } else {
                if ($1 == chr && $4 == ref && ($2 <= end + dist)) {
                    end = ($3 > end) ? $3 : end
                } else {
                    print chr, start, end, ref
                    chr=$1; start=$2; end=$3; ref=$4
                }
            }
        }
        END { if (NR>0) print chr, start, end, ref }
    ' > "$out_bed"
}

export -f run_blast_for_ref_to_bed remove_extensions multiply_fasta_python safe_id_from_path

# -------------------------
# Requirements
# -------------------------
need_cmd python3
need_cmd parallel
need_cmd blastn
need_cmd makeblastdb

# -------------------------
# USER INPUTS (English)
# -------------------------
read -e -p "Enter genome FASTA file(s) (space-separated): " input_genomes
read -p "How many chromosome sequences (contigs) should be used from each genome FASTA? " num_sequences
read -e -p "Enter reference FASTA file(s) (satDNA monomers OR a satelitome multi-FASTA) (space-separated): " refs_in
read -p "How many monomers should be used to build an array (repeat multiplier)? " multiplier
read -p "How many parallel jobs should run? (e.g., 4, 8): " NUM_JOBS

echo
echo "Choose how the 10 colored satDNAs are selected:"
echo "  1) First 10 sequences in the reference FASTA (FASTA order)"
echo "  2) Top 10 most abundant (by total bp covered in the genome hits)"
read -p "Type 1 or 2: " TOP_MODE

if [[ "$TOP_MODE" != "1" && "$TOP_MODE" != "2" ]]; then
    die "Invalid option for Top Mode. Please type 1 or 2."
fi

# Expand wildcards for refs
expanded_refs=()
for r in $refs_in; do
    matches=( $(compgen -G "$r" || true) )
    if [ ${#matches[@]} -eq 0 ]; then
        expanded_refs+=( "$r" )
    else
        expanded_refs+=( "${matches[@]}" )
    fi
done

if [ ${#expanded_refs[@]} -eq 0 ]; then
    die "No reference FASTA files were provided / found."
fi

echo
echo "Final reference file list:"
printf ' - %s\n' "${expanded_refs[@]}"
echo

# Avoid oversubscription: keep BLAST threads low when running multiple jobs
BLAST_THREADS=1

# -------------------------
# Main loop per genome
# -------------------------
for genome_fa in $input_genomes; do
    [ -f "$genome_fa" ] || die "Genome FASTA not found: $genome_fa"

    genome_name="$(remove_extensions "$genome_fa")"
    mkdir -p "$genome_name"

    # Temp workspace for this genome
    tmproot="$(mktemp -d -p "$genome_name" tmp_sat_density_XXXXXX)"
    cleanup() { rm -rf "$tmproot"; }
    trap cleanup EXIT

    temp_genome_fa="$tmproot/temp_genome.fasta"
    temp_bed_dir="$tmproot/bed_parts"
    mkdir -p "$temp_bed_dir"

    echo "[INFO] Processing genome: $genome_fa"
    echo "[INFO] Output folder: $genome_name"

    # Keep only first N sequences from genome FASTA (keeps the ORIGINAL FASTA ORDER)
    awk -v num_seq="$num_sequences" '
    BEGIN { count = 0 }
    /^>/ { if (count >= num_seq) exit; count++ }
    { print }
    ' "$genome_fa" > "$temp_genome_fa"

    [ -s "$temp_genome_fa" ] || die "Temporary genome FASTA is empty (check num_sequences)."

    # Build BLAST DB
    makeblastdb -in "$temp_genome_fa" -dbtype nucl -out "$tmproot/temp_genome_db" -parse_seqids >/dev/null 2>&1

    # One BED per reference -> safe parallel
    echo "[INFO] Running BLAST in parallel (jobs=$NUM_JOBS, blast_threads=$BLAST_THREADS per job)..."
    parallel --jobs "$NUM_JOBS" \
        run_blast_for_ref_to_bed \
        {} "$tmproot/temp_genome_db" "$multiplier" \
        "$temp_bed_dir/$(safe_id_from_path {}).bed" \
        "$BLAST_THREADS" "$tmproot" \
        ::: "${expanded_refs[@]}"

    # Merge all BED parts
    out_bed="$genome_name/valid_monomers.bed"
    cat "$temp_bed_dir"/*.bed 2>/dev/null | awk 'NF==4' > "$out_bed" || true

    if [ ! -s "$out_bed" ]; then
        echo "[WARNING] No hits found. valid_monomers.bed is empty for $genome_name."
    else
        echo "[OK] BED created: $out_bed"
    fi

    python3 - "$TOP_MODE" "$temp_genome_fa" "$out_bed" "$genome_name" "$genome_name" "${expanded_refs[@]}" <<'PY'
import sys
import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

TOP_MODE   = sys.argv[1]   # "1" or "2"
FASTA_FILE = sys.argv[2]
BED_FILE   = sys.argv[3]
OUT_DIR    = sys.argv[4]
GENOME_TAG = sys.argv[5]
REF_FASTAS = sys.argv[6:]

BIN_SIZE   = 100_000
TOP_N      = 10

# ---------------------------
# RULER SETTINGS
# ---------------------------
TICK_EVERY_BP  = 10_000_000   # 10 Mb
LABEL_EVERY_BP = 50_000_000   # 50 Mb
TICK_LEN_SMALL = 0.020
TICK_LEN_BIG   = 0.040

def read_bed_or_empty(path):
    if (not os.path.exists(path)) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["Chromosome","Start","End","Reference"])
    df = pd.read_csv(path, sep="\t", header=None, names=["Chromosome","Start","End","Reference"])
    df["Start"] = df["Start"].astype(int)
    df["End"]   = df["End"].astype(int)
    df["Reference"] = df["Reference"].astype(str)
    return df

def get_fasta_lengths_and_order(path):
    """
    IMPORTANT:
    - Preserves the ORIGINAL FASTA order (no sorting)
    - Returns: (chrom_order_list, lengths_dict)
    """
    lengths = {}
    order = []
    cur = None
    cur_len = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = cur_len
                cur = line[1:].split()[0]
                order.append(cur)
                cur_len = 0
            else:
                cur_len += len(line)
        if cur is not None:
            lengths[cur] = cur_len
    return order, lengths

def get_fasta_sequences(path, wanted_set=None):
    seqs = {}
    cur = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    if (wanted_set is None) or (cur in wanted_set):
                        seqs[cur] = "".join(chunks)
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if cur is not None:
            if (wanted_set is None) or (cur in wanted_set):
                seqs[cur] = "".join(chunks)
    return seqs

def compute_gc_track(chrom_order, chrom_lengths, chrom_seqs, bin_size=100_000):
    rows = []
    for chrom in chrom_order:
        clen = int(chrom_lengths[chrom])
        seq = chrom_seqs.get(chrom, "")
        if not seq:
            for bstart in range(0, clen, bin_size):
                bend = min(bstart + bin_size, clen)
                rows.append((chrom, bstart, bend, np.nan))
            continue

        seq_u = seq.upper()
        for bstart in range(0, clen, bin_size):
            bend = min(bstart + bin_size, clen)
            sub = seq_u[bstart:bend]
            if not sub:
                rows.append((chrom, bstart, bend, np.nan))
                continue
            g = sub.count("G"); c = sub.count("C"); a = sub.count("A"); t = sub.count("T")
            denom = a + t + g + c
            gc = (g + c) / denom if denom > 0 else np.nan
            rows.append((chrom, bstart, bend, gc))

    return pd.DataFrame(rows, columns=["Chromosome","Start","End","GC_fraction"])

def build_bins_for_chrom(clen, bin_size):
    return [(bstart, min(bstart + bin_size, clen)) for bstart in range(0, clen, bin_size)]

def read_first_n_refs_from_reference_fastas(ref_fastas, n=10):
    out = []
    seen = set()
    for fp in ref_fastas:
        if not fp or (not os.path.exists(fp)) or os.path.getsize(fp) == 0:
            continue
        with open(fp, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    rid = line[1:].strip().split()[0]
                    if rid and rid not in seen:
                        out.append(rid)
                        seen.add(rid)
                        if len(out) >= n:
                            return out
    return out

def choose_top_refs(bed_df, ref_fastas, top_mode, top_n=10):
    if str(top_mode) == "1":
        top_refs = read_first_n_refs_from_reference_fastas(ref_fastas, n=top_n)
        mode_label = "First 10 in reference FASTA order"
        return top_refs, mode_label

    if bed_df.empty:
        return [], "Top 10 most abundant (no hits found)"
    x = bed_df.copy()
    x["len"] = (x["End"] - x["Start"]).clip(lower=0)
    totals = x.groupby("Reference")["len"].sum().sort_values(ascending=False)
    top_refs = list(totals.head(top_n).index)
    mode_label = "Top 10 most abundant (by bp covered)"
    return top_refs, mode_label

def compute_density_by_ref(bed_df, chrom_order, chrom_lengths, top_refs, bin_size=100_000):
    top_refs = list(top_refs)[:TOP_N]
    ref_groups = top_refs + ["Other"]

    if bed_df.empty:
        rows = []
        for chrom in chrom_order:
            clen = int(chrom_lengths[chrom])
            for bstart, bend in build_bins_for_chrom(clen, bin_size):
                for rg in ref_groups:
                    rows.append((chrom, bstart, bend, rg, 0))
        dens_long = pd.DataFrame(rows, columns=["Chromosome","Start","End","Reference","Coverage_bp"])
        total_by_bin = dens_long.groupby(["Chromosome","Start","End"], as_index=False)["Coverage_bp"].sum()
        return dens_long, total_by_bin

    bed_df = bed_df.copy()
    bed_df["RefGroup"] = bed_df["Reference"].where(bed_df["Reference"].isin(top_refs), other="Other")

    intervals = {}
    for (chrom, rg), part in bed_df.groupby(["Chromosome","RefGroup"], sort=False):
        part = part.sort_values(["Start","End"])
        intervals[(chrom, rg)] = part[["Start","End"]].values.tolist()

    rows = []
    for chrom in chrom_order:
        clen = int(chrom_lengths[chrom])
        bins = build_bins_for_chrom(clen, bin_size)
        for (bstart, bend) in bins:
            bin_len = bend - bstart
            for rg in ref_groups:
                cov = 0
                ints = intervals.get((chrom, rg), [])
                for s, e in ints:
                    if e <= bstart:
                        continue
                    if s >= bend:
                        break
                    ov_s = max(s, bstart)
                    ov_e = min(e, bend)
                    if ov_e > ov_s:
                        cov += (ov_e - ov_s)
                if cov > bin_len:
                    cov = bin_len
                rows.append((chrom, bstart, bend, rg, int(cov)))

    dens_long = pd.DataFrame(rows, columns=["Chromosome","Start","End","Reference","Coverage_bp"])
    total_by_bin = dens_long.groupby(["Chromosome","Start","End"], as_index=False)["Coverage_bp"].sum()
    return dens_long, total_by_bin

def abbreviate_chrom_label(chrom_id):
    """
    ONLY abbreviate ChromosomeN -> Chr N (no reindexing, no reordering).
    Everything else stays the same (B1, B2, scaffold_3...).
    """
    m = re.match(r"^(Chromosome|chromosome)\s*0*([0-9]+)$", str(chrom_id))
    if m:
        return f"Chr {int(m.group(2))}"
    return str(chrom_id)

def plot_circos_like_density(dens_long, top_refs, total_by_bin, gc_df, chrom_order, chrom_lengths, out_png, title_suffix):
    def pastelize(rgb, amount=0.72):
        r, g, b = rgb[:3]
        return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)

    gap = math.radians(2.0)

    r_outer0, r_outer1 = 1.06, 1.20
    r_hist0,  r_hist1  = 0.80, 1.04
    r_gc0, r_gc1 = 0.66, 0.78

    lengths = np.array([chrom_lengths[c] for c in chrom_order], dtype=float)
    total_len = lengths.sum() if lengths.sum() > 0 else 1.0

    total_gap = gap * len(chrom_order)
    avail = 2 * math.pi - total_gap
    chrom_angles = avail * (lengths / total_len)

    cov_max = float(total_by_bin["Coverage_bp"].max()) if len(total_by_bin) else 1.0
    if cov_max <= 0:
        cov_max = 1.0

    gc_vals = gc_df["GC_fraction"].to_numpy(dtype=float) if (gc_df is not None and len(gc_df)) else np.array([])
    gc_vals = gc_vals[np.isfinite(gc_vals)]
    if gc_vals.size == 0:
        gc_min, gc_max = 0.0, 1.0
    else:
        gc_min, gc_max = float(gc_vals.min()), float(gc_vals.max())
        if gc_min == gc_max:
            gc_min, gc_max = 0.0, 1.0

    fig = plt.figure(figsize=(10, 10), dpi=300)
    ax = plt.subplot(111, projection="polar")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2)
    ax.set_axis_off()
    plt.rcParams["font.family"] = "DejaVu Sans"

    cmap_chr = plt.cm.viridis
    chrom_colors = [pastelize(c, amount=0.70) for c in cmap_chr(np.linspace(0.10, 0.90, len(chrom_order)))]

    top_refs = list(top_refs)[:TOP_N]
    ref_groups = top_refs + ["Other"]

    cmap_ref = plt.cm.tab10
    ref_colors = {}
    for i, rg in enumerate(ref_groups):
        if rg == "Other":
            ref_colors[rg] = (0.82, 0.82, 0.82, 0.90)
        else:
            col = cmap_ref(i % 10)
            ref_colors[rg] = (*pastelize(col, amount=0.55), 0.92)

    gc_line_color = (0.55, 0.78, 0.65, 0.95)

    def tangent_rotation_deg(theta, r):
        eps = 1e-4
        p1 = ax.transData.transform((theta, r))
        p2 = ax.transData.transform((theta + eps, r))
        dx, dy = (p2[0] - p1[0]), (p2[1] - p1[1])
        return math.degrees(math.atan2(dy, dx))

    def theta_for_pos(start_theta, ang, clen, pos):
        return start_theta + (pos / clen) * ang

    theta = 0.0
    for i, chrom in enumerate(chrom_order):
        clen = int(chrom_lengths[chrom])
        ang = chrom_angles[i]
        start_theta = theta
        end_theta   = theta + ang

        ax.bar(
            x=(start_theta + end_theta) / 2,
            height=(r_outer1 - r_outer0),
            width=(end_theta - start_theta),
            bottom=r_outer0,
            align="center",
            color=chrom_colors[i],
            edgecolor=(0.35, 0.35, 0.35),
            linewidth=0.6,
        )

        # Ruler ticks + labels
        r_tick_base = r_outer1
        if clen >= TICK_EVERY_BP:
            for pos in range(0, clen + 1, TICK_EVERY_BP):
                t = theta_for_pos(start_theta, ang, clen, pos)
                is_big = (pos % LABEL_EVERY_BP == 0)
                tick_len = TICK_LEN_BIG if is_big else TICK_LEN_SMALL
                lw = 0.75 if is_big else 0.45

                ax.plot([t, t], [r_tick_base, r_tick_base + tick_len],
                        linewidth=lw, color=(0.20, 0.20, 0.20), alpha=0.85)

                if is_big:
                    mb = int(pos / 1_000_000)
                    r_text = r_tick_base + tick_len + 0.014
                    rot = tangent_rotation_deg(t, r_text)
                    if rot < -90: rot += 180
                    elif rot > 90: rot -= 180

                    ax.text(t, r_text, f"{mb} Mb",
                            fontsize=6, color=(0.18, 0.18, 0.18),
                            ha="center", va="center",
                            rotation=rot, rotation_mode="anchor")

        # Sat stacked ring
        chrom_long = dens_long[dens_long["Chromosome"] == chrom].copy()
        if not chrom_long.empty:
            piv = chrom_long.pivot_table(
                index=["Start","End"],
                columns="Reference",
                values="Coverage_bp",
                aggfunc="sum",
                fill_value=0
            ).reset_index()

            starts = piv["Start"].to_numpy(dtype=float)
            ends   = piv["End"].to_numpy(dtype=float)

            theta0 = start_theta + (starts / clen) * ang
            theta1 = start_theta + (ends   / clen) * ang
            widths = np.maximum(theta1 - theta0, 1e-6)

            cols_order = [r for r in top_refs if r in piv.columns] + (["Other"] if "Other" in piv.columns else [])
            base = np.zeros(len(piv), dtype=float)
            ring_h = (r_hist1 - r_hist0)

            for rg in cols_order:
                cov = piv[rg].to_numpy(dtype=float)
                frac = np.clip(cov / cov_max, 0, 1)
                heights = ring_h * frac

                ax.bar(
                    x=(theta0 + theta1) / 2,
                    height=heights,
                    width=widths,
                    bottom=r_hist0 + base,
                    align="center",
                    color=ref_colors.get(rg, (0.7,0.7,0.7,0.8)),
                    edgecolor=None,
                    linewidth=0,
                    alpha=ref_colors.get(rg, (0,0,0,1))[3],
                )
                base += heights

        # GC inner line
        if gc_df is not None and len(gc_df):
            gchrom = gc_df[gc_df["Chromosome"] == chrom].sort_values(["Start","End"])
            if not gchrom.empty:
                gstarts = gchrom["Start"].to_numpy(dtype=float)
                gends   = gchrom["End"].to_numpy(dtype=float)
                gcs     = gchrom["GC_fraction"].to_numpy(dtype=float)

                centers = (gstarts + gends) / 2.0
                thetas = start_theta + (centers / clen) * ang

                gcs2 = np.array(gcs, dtype=float)
                med = float(np.nanmedian(gcs2)) if np.isfinite(gcs2).any() else 0.5
                gcs2[~np.isfinite(gcs2)] = med

                if gc_max > gc_min:
                    gcn = (gcs2 - gc_min) / (gc_max - gc_min)
                else:
                    gcn = gcs2 * 0 + 0.5

                radii = r_gc0 + gcn * (r_gc1 - r_gc0)
                ax.plot(thetas, radii, linewidth=1.0, color=gc_line_color, alpha=0.95)

        # Chrom label (ONLY abbreviation, no reindex)
        mid = (start_theta + end_theta) / 2
        r_label = (r_outer0 + r_outer1) / 2
        rot = tangent_rotation_deg(mid, r_label)
        if rot < -90: rot += 180
        elif rot > 90: rot -= 180

        arc_deg = math.degrees(end_theta - start_theta)
        fs = 8
        if arc_deg < 18: fs = 7
        if arc_deg < 12: fs = 6
        if arc_deg < 9:  fs = 5

        ax.text(mid, r_label, abbreviate_chrom_label(chrom),
                rotation=rot, rotation_mode="anchor",
                ha="center", va="center",
                fontsize=fs, fontweight="bold",
                color=(0.10,0.10,0.10), clip_on=True)

        theta = end_theta + gap

    ax.text(0, 1.38, f"Satelitome density (100 kb bins) — {GENOME_TAG}\n{title_suffix}",
            ha="center", va="center", fontsize=12.5, fontweight="bold")

    legend_handles = []
    for rg in top_refs[:TOP_N]:
        legend_handles.append(Patch(facecolor=ref_colors[rg], edgecolor="none", label=rg))
    legend_handles.append(Patch(facecolor=ref_colors["Other"], edgecolor="none", label="Other"))
    legend_handles.append(Line2D([0],[0], color=gc_line_color, lw=2, label="GC fraction (per 100 kb bin)"))

    fig.text(
        0.02, 0.02,
        "Tracks:\n- Outer ring: contigs/chromosomes (FASTA order; ChromosomeN abbreviated to Chr N)\n- SatDNA ring: stacked coverage per 100 kb bin\n- Inner line: GC fraction per 100 kb bin\n- Ruler: ticks every 10 Mb, labels every 50 Mb",
        ha="left", va="bottom", fontsize=8, color=(0.20, 0.20, 0.20)
    )

    ax.legend(
        handles=legend_handles,
        loc="lower right",
        bbox_to_anchor=(1.20, -0.02),
        frameon=True,
        framealpha=0.93,
        fontsize=7,
        title="Colored satDNAs (Top 10) + tracks",
        title_fontsize=8
    )

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

# -----------------------
# RUN
# -----------------------
bed_df = read_bed_or_empty(BED_FILE)

chrom_order, chrom_lengths = get_fasta_lengths_and_order(FASTA_FILE)
if not chrom_lengths:
    raise SystemExit("ERROR: Could not read any contig lengths from FASTA.")

wanted = set(chrom_order)
chrom_seqs = get_fasta_sequences(FASTA_FILE, wanted_set=wanted)

top_refs, mode_label = choose_top_refs(bed_df, REF_FASTAS, TOP_MODE, top_n=TOP_N)

print("[INFO] Top satDNAs used for coloring:")
if top_refs:
    for i, r in enumerate(top_refs, 1):
        print(f"  {i:02d} - {r}")
else:
    print("  (none; no hits or no references found)")

dens_long, total_by_bin = compute_density_by_ref(
    bed_df, chrom_order, chrom_lengths, top_refs=top_refs, bin_size=BIN_SIZE
)

gc_df = compute_gc_track(chrom_order, chrom_lengths, chrom_seqs, bin_size=BIN_SIZE)

dens_long_out = os.path.join(OUT_DIR, "satelitome_density_100kb_by_ref_long.tsv")
total_out     = os.path.join(OUT_DIR, "satelitome_density_100kb_total.tsv")
gc_out        = os.path.join(OUT_DIR, "gc_track_100kb_full.tsv")

dens_long.to_csv(dens_long_out, sep="\t", index=False)
total_by_bin.to_csv(total_out, sep="\t", index=False)
gc_df.to_csv(gc_out, sep="\t", index=False)

print("[OK] Density tables saved:")
print(f" - {dens_long_out}")
print(f" - {total_out}")
print("[OK] GC track table saved:")
print(f" - {gc_out}")

plot_out = os.path.join(OUT_DIR, "satelitome_density_circos_like.png")
title_suffix = f"Top 10 selection: {mode_label}"
plot_circos_like_density(dens_long, top_refs, total_by_bin, gc_df, chrom_order, chrom_lengths, plot_out, title_suffix)

print(f"[OK] Density plot saved: {plot_out}")
PY

    trap - EXIT
    cleanup
    echo
done

echo "[DONE] All genomes processed."
