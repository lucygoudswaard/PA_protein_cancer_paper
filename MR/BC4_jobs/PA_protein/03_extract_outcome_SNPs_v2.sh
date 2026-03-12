#!/bin/bash


# Directory containing .gz files
input_dir="/PA_protein/UKB_PPP_sumstats"

# Output file
output_file="${input_dir}/concatenated_output_overall_sedentary.txt"
> "$output_file"  # Clear the output file before starting

# Define SNP positions as a lookup table (CHR -> GENPOS list)
declare -A snp_positions
snp_positions["1"]="2234967 33225097 77984833"
snp_positions["3"]="68477984 18717009"
snp_positions["5"]="88689478 107487207 152566250 152659861 152675265 88653144"
snp_positions["7"]="72258898"
snp_positions["10"]="21531721 21541175"
snp_positions["11"]="104657953"
snp_positions["15"]="74039044"
snp_positions["17"]="45767194 46249498"
snp_positions["18"]="43188344"

# Convert SNP positions into a string format for awk (CHR|GENPOS)
awk_filter=""
for chrom in "${!snp_positions[@]}"; do
    for pos in ${snp_positions[$chrom]}; do
        awk_filter+="$chrom $pos|"
    done
done
awk_filter="${awk_filter%|}"  # Remove the last "|"

header_added=false  # Flag to track header inclusion

# Loop over each .gz file in the directory
for gz_file in "$input_dir"/*.gz; do
    base_name=$(basename "$gz_file" .gz)
    gene_name=$(echo "$base_name" | cut -d'_' -f1)  # Extract gene name

    # Process each file, filtering for matching CHROM & GENPOS
    if zgrep -qE "$awk_filter" "$gz_file"; then
        echo "Processing file: $gz_file"
        
        # Extract the header if not already added
        if [ "$header_added" = false ]; then
            zcat "$gz_file" | head -n 1 | awk '{print $0, "gene"}' >> "$output_file"
            header_added=true
            echo "Header added to output file."
        fi

        # Extract relevant rows and append the gene name column
        zgrep -E "$awk_filter" "$gz_file" | awk -v gene="$gene_name" '{print $0, gene}' >> "$output_file"
        
        echo "Rows with SNPs added from: $gz_file"
    else
        echo "No matching SNPs found in $gz_file"
    fi
done

echo "All matching rows have been concatenated into: $output_file"
