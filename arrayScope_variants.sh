#!/bin/bash

shopt -s progcomp

_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}
complete -F _autocomplete read

remove_extensions() {
    local filename="$1"
    filename=$(basename "$filename")
    filename="${filename%.fasta}"
    filename="${filename%.fa}"
    filename="${filename%.fna}"
    echo "$filename"
}

run_blast_for_ref() {
    local ref_no_ext="$1"
    local temp_genome="$2"
    local multiplier="$3"
    local temp_bed="$4"

    if [ -f "${ref_no_ext}.fasta" ]; then
        ref_fasta="${ref_no_ext}.fasta"
    elif [ -f "${ref_no_ext}.fa" ]; then
        ref_fasta="${ref_no_ext}.fa"
    elif [ -f "${ref_no_ext}.fna" ]; then
        ref_fasta="${ref_no_ext}.fna"
    else
        echo "Aviso: Nenhum arquivo fasta encontrado para ${ref_no_ext}"
        return
    fi

    local mult_ref="multiplied_reference_${ref_no_ext}.fasta"
    awk "NR % 2 == 0 { for (i=0;i<${multiplier};i++) printf \$0 } NR % 2 != 0" "$ref_fasta" > "$mult_ref"

    local blast_out="blast_${ref_no_ext}.out"
    blastn -task blastn -outfmt "6" -db "$temp_genome" \
        -query "$mult_ref" -out "$blast_out" \
        -evalue 1e-10 -qcov_hsp_perc 50 -num_threads 1

    if [ ! -s "$blast_out" ]; then
        echo "Sem hits para ${ref_no_ext}"
        return
    fi

    awk '{start=($9<$10)?$9:$10; end=($9<$10)?$10:$9; print $2, start, end}' "$blast_out" \
    | sort -k1,1 -k2,2n \
    | awk -v OFS="\t" -v dist=2000 '
        function print_block() {
            print chr, start, end
        }
        {
            if (NR == 1) { chr=$1; start=$2; end=$3 }
            else {
                if ($1 == chr && ($2 <= end + dist)) {
                    end = ($3 > end) ? $3 : end
                } else {
                    print_block()
                    chr=$1; start=$2; end=$3
                }
            }
        }
        END { print_block() }
    ' | awk -v ref="$ref_no_ext" '{split(ref, a, "_"); print $0"\t"a[1]}' >> "$temp_bed"
}
export -f run_blast_for_ref

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
    if [ -d "$genome_name" ]; then
        echo "Erro: diretório $genome_name já existe."
        exit 1
    fi
    mkdir -p "$genome_name" || { echo "Erro ao criar $genome_name"; exit 1; }

    temp_genome="$genome_name/temp_genome.fasta"
    awk -v num_seq="$num_sequences" '
    BEGIN { count = 0 }
    /^>/ { if (count >= num_seq) exit; count++ }
    { print }
    ' "$input_biblio" > "$temp_genome"

    if [ ! -f "$temp_genome" ]; then
        echo "Erro: $temp_genome não criado."
        exit 1
    fi

    makeblastdb -in "$temp_genome" -dbtype nucl -out "$temp_genome" -parse_seqids

    temp_bed="$genome_name/valid_monomers_temp.bed"
    > "$temp_bed"

    parallel --jobs "$NUM_THREADS" run_blast_for_ref {} "$temp_genome" "$multiplier" "$temp_bed" ::: "${expanded_refs_no_ext[@]}"

    mv "$temp_bed" "$genome_name/valid_monomers.bed"

    if [ ! -f "$genome_name/valid_monomers.bed" ]; then
        echo "Erro: arquivo BED não criado."
        exit 1
    fi

    echo "Arquivo valid_monomers.bed criado para $genome_name."

    # ------------------
    # Script Python embutido com as cores originais
    # ------------------
python3 - <<EOF
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import re
import hashlib

# Paleta fixa definida pelo usuário
distinct_colors = [
    (0.90, 0.12, 0.12), (0.12, 0.55, 0.12), (0.12, 0.12, 0.90), (0.90, 0.75, 0.12),
    (0.54, 0.17, 0.89), (0.89, 0.47, 0.20), (0.20, 0.89, 0.67), (0.20, 0.67, 0.89),
    (0.67, 0.20, 0.89), (0.89, 0.20, 0.67), (0.0, 0.5, 0.0),   (0.5, 0.0, 0.5),
    (0.0, 0.8, 0.8),   (0.8, 0.0, 0.0),   (0.3, 0.7, 0.9),   (0.9, 0.6, 0.0),
    (0.4, 0.2, 0.6),   (0.7, 0.3, 0.3),   (0.1, 0.9, 0.1),   (0.6, 0.4, 0.8),
    (0.95, 0.9, 0.1),  (0.8, 0.5, 0.9),   (0.5, 0.5, 0.5),   (0.0, 0.3, 0.6),
    (0.9, 0.7, 0.4),   (0.2, 0.8, 0.2),   (0.7, 0.0, 0.7),   (0.4, 0.6, 0.2),
    (0.9, 0.2, 0.5),   (0.3, 0.4, 0.9),   (0.6, 0.1, 0.1),   (0.1, 0.6, 0.6),
    (0.9, 0.4, 0.6),   (0.5, 0.3, 0.0),   (0.8, 0.8, 0.0),   (0.0, 0.7, 0.3),
    (0.7, 0.5, 0.0),   (0.29, 0.0, 0.51), (1.0, 0.41, 0.71), (0.0, 0.0, 0.4),
    (0.8, 0.2, 0.8),   (0.4, 0.0, 0.0),   (0.6, 0.8, 0.2),   (0.2, 0.2, 0.2),
    (0.94, 0.94, 0.86),(0.5, 0.0, 0.0),   (0.82, 0.71, 0.55),(0.9, 0.0, 0.9),
    (0.0, 0.5, 0.5)
]
big_palette = distinct_colors
num_colors = len(big_palette)

def get_color_for_ref(ref):
    h = hashlib.md5(ref.encode('utf-8')).hexdigest()
    idx = int(h, 16) % num_colors
    return big_palette[idx]

# Continuação do script (leitura BED, plots etc.) deve ser inserida aqui conforme o original
EOF
done

