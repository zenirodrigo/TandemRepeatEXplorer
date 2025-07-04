#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <gene_symbol> <arquivo_com_lista_de_especies>"
    exit 1
fi

GENE_SYMBOL=$1
ESPECIES_FILE=$2

if [ ! -f "$ESPECIES_FILE" ]; then
    echo "Arquivo de espécies '$ESPECIES_FILE' não encontrado."
    exit 1
fi

exec 3< "$ESPECIES_FILE"

while IFS= read -r SPECIES_NAME <&3 || [[ -n "$SPECIES_NAME" ]]; do
    echo "🔍 Buscando gene $GENE_SYMBOL para $SPECIES_NAME ..."

    SPECIES_QUERY=$(echo "$SPECIES_NAME" | sed 's/ /+/g')

    # Busca o Gene ID
    GENE_ID=$(esearch -db gene -query "$GENE_SYMBOL[Gene Name] AND $SPECIES_QUERY[Organism]" | efetch -format uid | head -1)

    if [ -z "$GENE_ID" ]; then
        echo "❌ Gene $GENE_SYMBOL não encontrado para $SPECIES_NAME."
        echo "-------------------------------------------"
        continue
    fi

    echo "🧠 Gene ID encontrado: $GENE_ID"

    # Nova abordagem para extrair coordenadas usando esummary e xtract
    COORDS=$(esummary -db gene -id "$GENE_ID" | \
             xtract -pattern DocumentSummary -element GenomicInfoType/ChrAccVer GenomicInfoType/ChrStart GenomicInfoType/ChrStop GenomicInfoType/ExonCount)

    if [ -z "$COORDS" ]; then
        echo "⚠️  Coordenadas não encontradas. Usando fallback para RNA..."
        # Fallback para RNA
        ZIP_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.zip"
        TMP_DIR="tmp_${GENE_SYMBOL}_${SPECIES_NAME// /_}"
        
        datasets download gene gene-id "$GENE_ID" --filename "$ZIP_FILE"
        
        if [ $? -ne 0 ]; then
            echo "⚠️  Falha ao baixar os dados para $SPECIES_NAME."
            echo "-------------------------------------------"
            continue
        fi

        mkdir -p "$TMP_DIR"
        unzip -q "$ZIP_FILE" -d "$TMP_DIR"

        # Tentar encontrar qualquer arquivo FASTA na pasta de dados
        FASTA_FILE=$(find "$TMP_DIR/ncbi_dataset/data" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | head -1)

        if [ -f "$FASTA_FILE" ]; then
            OUTPUT_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.fasta"
            cp "$FASTA_FILE" "$OUTPUT_FILE"
            echo "✅ Dados salvos como $OUTPUT_FILE"
        else
            echo "❌ Nenhum arquivo FASTA encontrado para $SPECIES_NAME."
        fi

        rm -rf "$TMP_DIR"
        rm -f "$ZIP_FILE"
        echo "-------------------------------------------"
        continue
    fi

    # Processar coordenadas
    CONTIG=$(echo "$COORDS" | cut -f1)
    START=$(echo "$COORDS" | cut -f2)
    END=$(echo "$COORDS" | cut -f3)
    EXON_COUNT=$(echo "$COORDS" | cut -f4)

    # Verificar se temos coordenadas válidas
    if [[ ! "$START" =~ ^[0-9]+$ ]] || [[ ! "$END" =~ ^[0-9]+$ ]]; then
        echo "⚠️  Coordenadas inválidas: $START-$END. Usando fallback para RNA..."
        # ... (código fallback igual acima)
        continue
    fi

    echo "🗺️ Coordenadas encontradas: $CONTIG:$START-$END (Exons: $EXON_COUNT)"

    # Baixar diretamente do NCBI usando efetch
    OUTPUT_FILE="${GENE_SYMBOL}_${SPECIES_NAME// /_}.fasta"
    
    # Determinar orientação e ajustar coordenadas
    if [ "$START" -gt "$END" ]; then
        TEMP=$START
        START=$END
        END=$TEMP
        REVERSE=1
    else
        REVERSE=0
    fi

    echo "🌐 Baixando sequência genômica (${START}-${END})..."
    efetch -db nuccore -id "$CONTIG" -seq_start "$START" -seq_stop "$END" -format fasta > "$OUTPUT_FILE.tmp"
    
    # Aplicar reverso complementar se necessário
    if [ "$REVERSE" -eq 1 ]; then
        seqkit seq -r -p -t dna "$OUTPUT_FILE.tmp" > "$OUTPUT_FILE"
        rm "$OUTPUT_FILE.tmp"
        echo "🔄 Fita invertida (strand negativo)"
    else
        mv "$OUTPUT_FILE.tmp" "$OUTPUT_FILE"
    fi

    echo "✅ Sequência genômica salva como $OUTPUT_FILE"
    echo "Tamanho da sequência: $(seqkit stats -T "$OUTPUT_FILE" | tail -n1 | cut -f5)"
    echo "-------------------------------------------"

done

exec 3<&-

echo "🏁 Processo concluído para todas as espécies."
