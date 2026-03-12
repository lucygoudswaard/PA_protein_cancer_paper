#!/bin/bash

# --- Set directories ---
input_dir="/PA_protein/UKB_PPP_sumstats"
snp_file="/PA_protein/cancer_sumstats/cancer_clump.txt"
output_file="${input_dir}/concatenated_output_cancer.txt"

# --- Clear previous output ---
> "$output_file"

# --- Create chr:pos column for SNP lookup ---
awk 'NR > 1 { print $7":"$8 }' "$snp_file" > snp_chrpos.txt

# --- Set flag for header ---
header_added=false

# --- Loop over each .gz file ---
for gz_file in "$input_dir"/*.gz; do
    base_name=$(basename "$gz_file" .gz)
    gene_name=$(echo "$base_name" | cut -d'_' -f1)

    # Check if file contains any matching SNPs by chr:pos in ID column
    if zcat "$gz_file" | awk '
        BEGIN {
            while ((getline < "snp_chrpos.txt") > 0) snp[$1] = 1
        }
        NR > 1 {
            split($3, parts, ":")
            key = parts[1] ":" parts[2]
            if (key in snp) { found=1; exit }
        }
        END { exit !found }
    '; then
        echo "Processing file: $gz_file"

        # Add header once
        if [ "$header_added" = false ]; then
            zcat "$gz_file" | head -n 1 | awk '{print $0, "gene"}' >> "$output_file"
            header_added=true
            echo "Header added to output file."
        fi

        # Append matching rows with gene name
        zcat "$gz_file" | awk -v gene="$gene_name" '
            BEGIN {
                while ((getline < "snp_chrpos.txt") > 0) snp[$1] = 1
            }
            NR == 1 { next }  # skip header
            {
                split($3, parts, ":")
                key = parts[1] ":" parts[2]
                if (key in snp) {
                    print $0, gene
                }
            }
        ' >> "$output_file"

        echo "Rows with matching SNPs added from: $gz_file"
    else
        echo "No matching SNPs found in $gz_file"
    fi
done

# --- Remove temporary file ---
rm snp_chrpos.txt

echo "All matching rows have been concatenated into: $output_file"
