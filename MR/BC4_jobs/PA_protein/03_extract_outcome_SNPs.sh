#!/bin/bash

# Set the directory path where the .gz files are located
input_dir="/PA_protein/UKB_PPP_sumstats"

# SNP positions to search for
positions="45767194|43188344"

# Temporary file to store concatenated results
output_file="${input_dir}/concatenated_output.txt"
> "$output_file"  # Ensure the output file is empty initially

# Flag to track if the header has been added
header_added=false

# Loop over each .gz file in the directory
for gz_file in "$input_dir"/*.gz; do
    # Extract the base filename (without the .gz extension)
    base_name=$(basename "$gz_file" .gz)
    
    # Get the gene name from the file (first part before "_")
    gene_name=$(echo "$base_name" | cut -d'_' -f1)

    # Check if positions exist in the file and extract the rows
    if zgrep -qE "$positions" "$gz_file"; then
        echo "Processing file: $gz_file"
        
        # Extract the header if not already added
        if [ "$header_added" = false ]; then
            # Add the header from the first file to the output file, and append "gene" as the last column
            header=$(zcat "$gz_file" | head -n 1)
            echo -e "${header}\tgene" >> "$output_file"  # Append 'gene' column to header
            header_added=true
            echo "Header added to output file."
        fi
        
        # Extract rows containing the SNP positions and add the 'outcome' column (gene)
        zgrep -E "$positions" "$gz_file" | awk -v gene="$gene_name" 'BEGIN {FS=OFS="\t"} {print $0, gene}' >> "$output_file"
        
        echo "Rows with SNPs added to: $output_file"
    else
        echo "No matching SNPs found in $gz_file"
    fi
done

# Inform the user about the concatenated output
echo "All matching rows have been concatenated into: $output_file"
