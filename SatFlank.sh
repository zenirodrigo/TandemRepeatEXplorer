#!/bin/bash

########################################
# User Input Collection (pre-execution)
########################################

# BLAST Part:
read -e -p "Enter the library file name: " input_biblio
read -e -p "Enter the reference file name: " referencia
read -p "How many times to multiply the reference sequence? " multiplicador
read -e -p "Enter the BLAST output file name: " saida_blast
read -e -p "Enter the final output file name for extraction (seqtk): " output_seqtk

# GFF/Comparison Part:
read -p "Number of CPU cores to use for parallel processing: " nucleos
read -e -p "Enter the .gff file name: " arquivo_gff
if [ ! -f "$arquivo_gff" ]; then
    echo ".gff file not found."
    exit 1
fi
read -e -p "Enter the final output file name for GFF analysis: " arquivo_final

########################################
# Initial Configuration
########################################

# Enable filename autocompletion
shopt -s progcomp
_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}
complete -F _autocomplete read

########################################
# Part 1: BLAST and .bed File Generation
########################################

# Multiply the reference sequence
awk 'NR % 2 == 0 { for (i=0; i<'$multiplicador'; i++) printf $0 } NR % 2 != 0' "$referencia" > referencia_multiplicada.fasta

# Prepare the library for BLAST
makeblastdb -in "$input_biblio" -dbtype nucl

# Run blastn with tabular output format (format 6)
blastn -task blastn -outfmt "6" -db "$input_biblio" -query referencia_multiplicada.fasta -out "$saida_blast" -evalue 1e-10 -qcov_hsp_perc 70

# Check if BLAST output is not empty
if [ ! -s "$saida_blast" ]; then
    echo "BLAST found no matches for the reference."
    exit 1
fi

# Generate .bed file with valid monomers (pre-filtering)
awk '{start=($9 < $10) ? $9 : $10; end=($9 < $10) ? $10 : $9; print $2, start - 5000, end + 5000}' "$saida_blast" > monomeros_validos.bed

########################################
# Filter .bed File to Remove Redundant Regions
########################################

# Sort .bed by chromosome and start coordinate
sort -k1,1 -k2,2n monomeros_validos.bed > monomeros_validos.sorted.bed

# For overlapping intervals on the same chromosome, keep the longest interval.
awk 'BEGIN { OFS="\t" }
{
  chrom = $1; start = $2; end = $3; len = end - start
  if (NR == 1) {
    best_chrom = chrom; best_start = start; best_end = end; best_len = len
  } else {
    if (chrom == best_chrom && start <= best_end) {
      if (len > best_len) {
         best_start = start; best_end = end; best_len = len
      }
      # Discard shorter/equal intervals
    } else {
      print best_chrom, best_start, best_end
      best_chrom = chrom; best_start = start; best_end = end; best_len = len
    }
  }
}
END { print best_chrom, best_start, best_end }' monomeros_validos.sorted.bed > monomeros_validos.filtrado.bed

echo "Generated and filtered .bed file: monomeros_validos.filtrado.bed"

# Use the filtered .bed file for extraction with seqtk
seqtk subseq "$input_biblio" monomeros_validos.filtrado.bed > "$output_seqtk"
echo "BLAST step completed. Extracted sequences saved to $output_seqtk"

########################################
# Part 2: .bed File Comparison with .gff
########################################

# Use the filtered .bed file from Part 1:
arquivo_bed="monomeros_validos.filtrado.bed"

# Convert spaces to tabs in .bed (ensure correct format)
sed -i 's/ \+/\t/g' "$arquivo_bed"

# Preprocess GFF: create a temporary indexed file by chromosome
gff_indexado=$(mktemp)
awk -F'\t' 'BEGIN {OFS=FS} 
    $3 == "gene" {
        gsub(/\r/, "", $0)
        split($9, attr, ";")
        symbol = "No symbol"
        
        # Primeiro verifica por "symbol="
        for (i in attr) {
            if (attr[i] ~ /symbol=/) {
                split(attr[i], parts, /[=]/)
                symbol = parts[2]
                gsub(/"/, "", symbol)
                break
            }
        }
        
        # Se não encontrou "symbol=", procura por "gene="
        if (symbol == "No symbol") {
            for (i in attr) {
                if (attr[i] ~ /gene=/) {
                    split(attr[i], parts, /[=]/)
                    symbol = parts[2]
                    gsub(/"/, "", symbol)
                    break
                }
            }
        }
        
        print $1, $4, $5, symbol
    }' "$arquivo_gff" | sort -k1,1 -k2,2n > "$gff_indexado"
# Parallel processing function for each .bed region
processar_regiao() {
    local cromossomo=$1
    local inicio=$2
    local fim=$3
    local tmpdir=$4
    local region="${inicio}-${fim}"
    local tmp_file
    tmp_file=$(mktemp -p "$tmpdir")
    
    # Search indexed GFF: print overlapping genes in tabular format
    awk -F'\t' -v crom="$cromossomo" -v region="$region" -v ini="$inicio" -v fim="$fim" '
        $1 == crom && !( $3 < ini || $2 > fim ) {
            gene = ($4 == "" ? "No symbol" : $4)
            print crom "\t" region "\t" gene
        }
    ' "$gff_indexado" > "$tmp_file"
    
    # Return non-empty temporary files for concatenation
    if [ -s "$tmp_file" ]; then
        echo "$tmp_file"
    fi
}

export -f processar_regiao
export gff_indexado

# Create temporary directory for partial results
tmpdir=$(mktemp -d)

# Process .bed lines in parallel and generate temporary result files
parallel --jobs "$nucleos" --colsep '\t' \
    processar_regiao {1} {2} {3} "$tmpdir" < <(tr -d '\r' < "$arquivo_bed") > lista_temp_files.txt

# Create final file with header
echo -e "Chromosome\tregion\tgene" > "$arquivo_final"

# Combine results from temporary files (if any)
if [ -s lista_temp_files.txt ]; then
    while IFS= read -r file; do
        cat "$file" >> "$arquivo_final"
    done < lista_temp_files.txt
fi

# Cleanup temporary files
rm -rf "$tmpdir" lista_temp_files.txt
rm "$gff_indexado"

echo "GFF analysis completed. Results saved to: $arquivo_final"
