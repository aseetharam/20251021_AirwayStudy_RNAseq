#!/bin/bash
# files can be edited directly from the vscode!
# 00_get-raw-files.sh
#
# Downloads, compresses, and renames a list of SRA runs for the HASM project.
# This script requires sra-tools (for fasterq-dump) and pigz to be installed.
# It should be run from the project's base directory.

# --- Configuration ---
# Stop on error, unset variable, or pipe failure
set -euo pipefail

# Define project base directory and output directory
BASE_DIR=$(pwd) # Assumes script is run from the project root
RAW_DATA_DIR="${BASE_DIR}/01_data/raw"

# Number of threads for downloading and compression.
THREADS=8

# Use an associative array to map SRR accessions to their descriptive names.
declare -A SAMPLES
SAMPLES=(
    ["SRR1039508"]="N61311_untreated"
    ["SRR1039509"]="N61311_Dex"
    ["SRR1039510"]="N61311_Alb"
    ["SRR1039511"]="N61311_Alb_Dex"
    ["SRR1039512"]="N052611_untreated"
    ["SRR1039513"]="N052611_Dex"
    ["SRR1039514"]="N052611_Alb"
    ["SRR1039515"]="N052611_Alb_Dex"
    ["SRR1039516"]="N080611_untreated"
    ["SRR1039517"]="N080611_Dex"
    ["SRR1039518"]="N080611_Alb"
    ["SRR1039519"]="N080611_Alb_Dex"
    ["SRR1039520"]="N061011_untreated"
    ["SRR1039521"]="N061011_Dex"
    ["SRR1039522"]="N061011_Alb"
    ["SRR1039523"]="N061011_Alb_Dex"
)

# --- Setup ---
# Create the raw data directory if it doesn't exist
mkdir -p "$RAW_DATA_DIR"
echo "INFO: Raw data will be placed in ${RAW_DATA_DIR}"
echo "--------------------------------------------------"

# --- Main ---
# Process each sample defined in the associative array.
for srr in "${!SAMPLES[@]}"; do
    echo "INFO: Processing accession: ${srr}"
    desc="${SAMPLES[$srr]}"
    
    final_r1_path="${RAW_DATA_DIR}/${srr}_${desc}_1.fastq.gz"
    
    # Check if the final file already exists to make the script resumable
    if [ -f "$final_r1_path" ]; then
        echo "  -> Final file already exists. Skipping."
        continue
    fi

    # Download the paired-end FASTQ files into the current directory.
    echo "  -> Downloading with fasterq-dump..."
    fasterq-dump --split-files --threads "${THREADS}" --progress "${srr}"

    # Compress the output files using parallel gzip.
    echo "  -> Compressing with pigz..."
    pigz --force --processes "${THREADS}" "${srr}_1.fastq"
    pigz --force --processes "${THREADS}" "${srr}_2.fastq"

    # Rename and move the compressed files to their final destination.
    echo "  -> Renaming and moving files for: ${desc}"
    mv "${srr}_1.fastq.gz" "${RAW_DATA_DIR}/${srr}_${desc}_1.fastq.gz"
    mv "${srr}_2.fastq.gz" "${RAW_DATA_DIR}/${srr}_${desc}_2.fastq.gz"
    
    echo "INFO: Completed processing for ${srr}."
    echo "--------------------------------------------------"
done

echo "INFO: All download tasks completed successfully."
