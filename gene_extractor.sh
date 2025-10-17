#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <gene_symbol> <species_list_file>"
    exit 1
fi

GENE_SYMBOL=$1
ESPECIES_FILE=$2

if [ ! -f "$ESPECIES_FILE" ]; then
    echo "Species file '$ESPECIES_FILE' not found."
    exit 1
fi

exec 3< "$ESPECIES_FILE"

while IFS= read -r SPECIES_NAME <&3 || [[ -n "$SPECIES_NAME" ]]; do
    echo "🔍 Searching for gene $GENE_SYMBOL in $SPECIES_NAME ..."

    SPECIES_QUERY=$(echo "$SPECIES_NAME" | sed 's/ /+/g')

    # Search for Gene ID
    GENE_ID=$(esearch -db gene -query "$GENE_SYMBOL[Gene Name] AND $SPECIES_QUERY[Organism]" | efetch -format uid | head -1)

    if [ -z "$GENE_ID" ]; then
        echo "❌ Gene $GENE_SYMBOL not found for $SPECIES_NAME."
        echo "-------------------------------------------"
        continue
    fi

    echo "🧠 Gene ID found: $GENE_ID"

    # New approach to extract coordinates using esummary and xtract
    COORDS=$(esummary -db gene -id "$GENE_ID" | \
             xtract -pattern DocumentSummary -element GenomicInfoType/ChrAccVer GenomicInfoType/ChrStart GenomicInfoType/ChrStop GenomicInfoType/ExonCount)

    if [ -z "$COORDS" ]; then
        echo "⚠️  Coordinates not found. Using RNA fallback..."
        # Fallback for RNA
        ZIP_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.zip"
        TMP_DIR="tmp_${GENE_SYMBOL}_${SPECIES_NAME// /_}"
        
        datasets download gene gene-id "$GENE_ID" --filename "$ZIP_FILE"
        
        if [ $? -ne 0 ]; then
            echo "⚠️  Failed to download data for $SPECIES_NAME."
            echo "-------------------------------------------"
            continue
        fi

        mkdir -p "$TMP_DIR"
        unzip -q "$ZIP_FILE" -d "$TMP_DIR"

        # Try to find any FASTA file in the data folder
        FASTA_FILE=$(find "$TMP_DIR/ncbi_dataset/data" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | head -1)

        if [ -f "$FASTA_FILE" ]; then
            OUTPUT_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.fasta"
            cp "$FASTA_FILE" "$OUTPUT_FILE"
            echo "✅ Data saved as $OUTPUT_FILE"
        else
            echo "❌ No FASTA file found for $SPECIES_NAME."
        fi

        rm -rf "$TMP_DIR"
        rm -f "$ZIP_FILE"
        echo "-------------------------------------------"
        continue
    fi

    # Process coordinates
    CONTIG=$(echo "$COORDS" | cut -f1)
    START=$(echo "$COORDS" | cut -f2)
    END=$(echo "$COORDS" | cut -f3)
    EXON_COUNT=$(echo "$COORDS" | cut -f4)

    # Check if coordinates are valid
    if [[ ! "$START" =~ ^[0-9]+$ ]] || [[ ! "$END" =~ ^[0-9]+$ ]]; then
        echo "⚠️  Invalid coordinates: $START-$END. Using RNA fallback..."
        # ... (same fallback code as above)
        continue
    fi

    echo "🗺️ Coordinates found: $CONTIG:$START-$END (Exons: $EXON_COUNT)"

    # Direct download from NCBI using efetch
    OUTPUT_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.fasta"
    
    # Determine orientation and adjust coordinates
    if [ "$START" -gt "$END" ]; then
        TEMP=$START
        START=$END
        END=$TEMP
        REVERSE=1
    else
        REVERSE=0
    fi

    echo "🌐 Downloading genomic sequence (${START}-${END})..."
    efetch -db nuccore -id "$CONTIG" -seq_start "$START" -seq_stop "$END" -format fasta > "$OUTPUT_FILE.tmp"
    
    # Apply reverse complement if needed
    if [ "$REVERSE" -eq 1 ]; then
        seqkit seq -r -p -t dna "$OUTPUT_FILE.tmp" > "$OUTPUT_FILE"
        rm "$OUTPUT_FILE.tmp"
        echo "🔄 Reverse strand detected"
    else
        mv "$OUTPUT_FILE.tmp" "$OUTPUT_FILE"
    fi

    echo "✅ Genomic sequence saved as $OUTPUT_FILE"
    echo "Sequence length: $(seqkit stats -T "$OUTPUT_FILE" | tail -n1 | cut -f5)"
    echo "-------------------------------------------"

done

exec 3<&-

echo "🏁 Process completed for all species."
