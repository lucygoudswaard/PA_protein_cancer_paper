#!/bin/bash

DIRECTORY_PROCESSED=/PA_protein/
DIRECTORY_ANCESTRY="UKB_PPP_sumstats/"
cd "${DIRECTORY_PROCESSED}${DIRECTORY_ANCESTRY}"

tmp=$(mktemp) || { ret="$?"; printf 'Failed to create temp file\n'; exit "$ret"; }
# Loop through the files and unzip
for FILE_NAME_EXTENSION in *.tar; do
  DIRECTORY="${FILE_NAME_EXTENSION%.tar}"
tar -xvf ${FILE_NAME_EXTENSION}
# extract GWAS name
## Set the pattern to match and extract the desired string
pattern="(.*)_v[0-9]+_"
## Loop through the folders and extract the string
  if [[ ${DIRECTORY} =~ $pattern ]]; then
    GWAS_NAME="${BASH_REMATCH[1]}"
    result="$GWAS_NAME"
  fi
# unzip GWAS
echo "unzip"
gzip -d ${DIRECTORY}/*
# cat GWAS
echo "cat"
cat ${DIRECTORY}/* > ${DIRECTORY}/${GWAS_NAME}
# move to processed
echo "move"
mv ${DIRECTORY}/${GWAS_NAME} ${DIRECTORY_PROCESSED}${DIRECTORY_ANCESTRY}
# zip
echo "zip"
gzip ${DIRECTORY_PROCESSED}${DIRECTORY_ANCESTRY}${GWAS_NAME}
# clean up
rm -rf ${DIRECTORY}/
rm -rf ${DIRECTORY}.tar
done
