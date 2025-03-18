#!/bin/bash

########################################
# Inputs
########################################

# BLAST processing:
read -e -p "Enter genome file name: " input_lib
read -e -p "Enter reference file name: " reference
read -p "How many times to multiply the reference monomer to form a array (here is the minimum size for the array? " multiplier
read -e -p "Enter BLAST output file name: " blast_output
read -e -p "Enter final sequence extraction file name (seqtk): " seq_output

# GFF comparison:
read -p "Number of CPU cores for parallel processing: " cores
read -e -p "Enter .gff file name: " gff_file
if [ ! -f "$gff_file" ]; then
    echo ".gff file not found."
    exit 1
fi
read -e -p "Enter final GFF analysis output file name: " final_output

########################################
# Initial Setup
########################################

# Enable filename autocompletion
shopt -s progcomp
_autocomplete() {
    local cur=${COMP_WORDS[COMP_CWORD]}
    COMPREPLY=( $(compgen -f -- "$cur") )
}
complete -F _autocomplete read

########################################
# Part 1: BLAST and .bed Generation
########################################

# Multiply reference sequence
awk 'NR % 2 == 0 { for (i=0; i<'$multiplier'; i++) printf $0 } NR % 2 != 0' "$reference" > multiplied_reference.fasta

# Prepare BLAST database
makeblastdb -in "$input_lib" -dbtype nucl

# Run BLASTn with tabular output (format 6)
blastn -task blastn -outfmt "6" -db "$input_lib" -query multiplied_reference.fasta -out "$blast_output" -evalue 1e-50 -qcov_hsp_perc 70

# Check BLAST results
if [ ! -s "$blast_output" ]; then
    echo "BLAST found no matches for the reference."
    exit 1
fi

# Generate initial .bed file with valid monomers
awk '{start=($9 < $10) ? $9 : $10; end=($9 < $10) ? $10 : $9; print $2, start - 5000, end + 5000}' "$blast_output" > valid_monomers.bed

########################################
# Filter .bed file to remove redundant regions
########################################

# Sort .bed by chromosome and start position
sort -k1,1 -k2,2n valid_monomers.bed > valid_monomers.sorted.bed

# Keep only longest interval in overlapping regions
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
    } else {
      print best_chrom, best_start, best_end
      best_chrom = chrom; best_start = start; best_end = end; best_len = len
    }
  }
}
END { print best_chrom, best_start, best_end }' valid_monomers.sorted.bed > valid_monomers.filtered.bed

echo "Filtered .bed file generated: valid_monomers.filtered.bed"

# Extract sequences using filtered .bed
seqtk subseq "$input_lib" valid_monomers.filtered.bed > "$seq_output"
echo "BLAST processing complete. Extracted sequences saved to $seq_output"

########################################
# Part 2: .bed/.gff Comparison
########################################

bed_file="valid_monomers.filtered.bed"

# Convert spaces to tabs in .bed
sed -i 's/ \+/\t/g' "$bed_file"

# Preprocess GFF: create temp indexed file
gff_indexed=$(mktemp)
awk -F'\t' 'BEGIN {OFS=FS} 
    $3 == "gene" {
        gsub(/\r/, "", $0)
        split($9, attr, ";")
        symbol = ""
        for (i in attr) {
            if (attr[i] ~ /symbol=/) {
                split(attr[i], parts, "=")
                symbol = parts[2]
                gsub(/"/, "", symbol)
                break
            }
        }
        print $1, $4, $5, symbol
    }' "$gff_file" | sort -k1,1 -k2,2n > "$gff_indexed"

# Parallel processing function
process_region() {
    local chrom=$1
    local start=$2
    local end=$3
    local tmpdir=$4
    local region="${start}-${end}"
    local tmp_file
    tmp_file=$(mktemp -p "$tmpdir")
    
    # Find overlapping genes
    awk -F'\t' -v c="$chrom" -v s="$start" -v e="$end" '
        $1 == c && !( $3 < s || $2 > e ) {
            gene = ($4 == "" ? "No_symbol" : $4)
            print c "\t" s"-"e "\t" gene
        }
    ' "$gff_indexed" > "$tmp_file"
    
    [ -s "$tmp_file" ] && echo "$tmp_file"
}

export -f process_region
export gff_indexed

# Create temp directory
tmpdir=$(mktemp -d)

# Process in parallel
parallel --jobs "$cores" --colsep '\t' \
    process_region {1} {2} {3} "$tmpdir" < <(tr -d '\r' < "$bed_file") > temp_files.txt

# Create final output
echo -e "Chromosome\tRegion\tGene" > "$final_output"

# Combine results
[ -s temp_files.txt ] && while IFS= read -r file; do
    cat "$file" >> "$final_output"
done < temp_files.txt

# Cleanup
rm -rf "$tmpdir" temp_files.txt
rm "$gff_indexed"

echo "GFF analysis complete. Results saved to: $final_output"
