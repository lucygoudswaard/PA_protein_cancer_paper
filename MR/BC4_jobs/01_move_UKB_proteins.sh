#!/bin/bash

# Define the path to the protein names list file
protein_names_file="/PA_protein/protein_names_MR.txt"

# Define the source and destination directories
source_directory="/UKB-PPP/UKB-PPP pGWAS summary statistics/European (discovery)/"
destination_directory="/PA_protein/UKB_PPP_sumstats/"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_directory"

# Initialize counters
copied_files_count=0
error_count=0

# Loop through each protein name in the protein names file
while IFS= read -r protein_name; do
  # Trim any leading/trailing whitespace from protein_name
  protein_name=$(echo "$protein_name" | xargs)

  # Skip empty protein names (if any)
  if [[ -z "$protein_name" ]]; then
    continue
  fi

  # Use find to locate files that match the protein name exactly before the first underscore
  matching_files=$(find "$source_directory" -maxdepth 1 -type f -iname "$protein_name"_*)

  # Check if matching files were found
  if [[ -z "$matching_files" ]]; then
    echo "No files found for protein: $protein_name" >> "$destination_directory/copy_errors.log"
    continue
  fi

  # Loop through the found files and copy them
  while IFS= read -r file; do
    # Check if the part before the first underscore in the filename matches the protein name
    base_filename=$(basename "$file")
    gene_part="${base_filename%%_*}"

    if [[ "$gene_part" == "$protein_name" ]]; then
      if cp "$file" "$destination_directory"; then
        copied_files_count=$((copied_files_count + 1))  # Increment the counter for successful copies
      else
        echo "Failed to copy: $file" >> "$destination_directory/copy_errors.log"
        error_count=$((error_count + 1))  # Increment the error counter
      fi
    else
      # If the part before the underscore doesn't match, skip this file
      echo "File does not match the protein name: $file" >> "$destination_directory/copy_errors.log"
    fi
  done <<< "$matching_files"
done < "$protein_names_file"

# Final log summary
echo "Copy operation completed."
echo "Total files copied: $copied_files_count"
echo "Total errors: $error_count"

if [ "$copied_files_count" -gt 0 ]; then
  echo "Successfully copied $copied_files_count files." >> "$destination_directory/copy_summary.log"
else
  echo "No files were copied. Please check the error logs for details." >> "$destination_directory/copy_summary.log"
fi
