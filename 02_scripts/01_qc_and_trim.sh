#!/bin/bash
# files can be edited directly!
# 01_qc_and_trim.sh
#
# This script performs two main functions:
# 1. Runs FastQC on the raw FASTQ files to assess initial quality.
# 2. Uses Trimmomatic to remove adapters and low-quality reads.
#
# It expects the project directory structure to be in place.

# --- Configuration ---
# Stop script on any error
set -e

# Define project base directory relative to the script location
BASE_DIR="/scratch/scholar/aseethar/20250915_AirwayStudy_RNAseq"

# Input/Output Directories
RAW_DATA_DIR="${BASE_DIR}/01_data/raw"
PROCESSED_DATA_DIR="${BASE_DIR}/01_data/processed"
LOG_DIR="${BASE_DIR}/99_logs"
FASTQC_RAW_DIR="${PROCESSED_DATA_DIR}/fastqc_raw"
TRIMMED_READS_DIR="${PROCESSED_DATA_DIR}/trimmed_reads"

# Trimmomatic specific settings
# NOTE: You must provide the path to your Trimmomatic JAR file and adapter sequences.
TRIMMOMATIC_JAR="/home/aseethar/trimmomatic-0.40/trimmomatic-0.40.jar"
ADAPTERS="//home/aseethar/trimmomatic-0.40/adapters/TruSeq3-PE-2.fa"

# Number of threads to use
THREADS=32

# --- Setup ---
# Create output directories if they don't exist
mkdir -p "$FASTQC_RAW_DIR"
mkdir -p "$TRIMMED_READS_DIR"
mkdir -p "$LOG_DIR"

echo "Starting QC and Trimming Pipeline..."
echo "Project Base Directory: ${BASE_DIR}"
echo "------------------------------------"

# --- Step 1: Run FastQC on Raw Data ---
echo "Step 1: Running FastQC on raw data..."
ml --force purge
ml biocontainers
ml fastqc
fastqc --threads "$THREADS" --outdir "$FASTQC_RAW_DIR" "$RAW_DATA_DIR"/*.fastq.gz
echo "FastQC on raw data complete. Reports are in ${FASTQC_RAW_DIR}"
echo "------------------------------------"

# --- Step 2: Run Trimmomatic ---
echo "Step 2: Running Trimmomatic for adapter and quality trimming..."

for r1_file in "$RAW_DATA_DIR"/*_1.fastq.gz; do
    # Derive the R2 filename and the base sample name from the R1 file
    r2_file="${r1_file/_1.fastq.gz/_2.fastq.gz}"
    sample_base=$(basename "$r1_file" _1.fastq.gz)
    
    echo "Processing sample: ${sample_base}"

    # Define output files for Trimmomatic
    r1_paired_out="${TRIMMED_READS_DIR}/${sample_base}_1_paired.fastq.gz"
    r1_unpaired_out="${TRIMMED_READS_DIR}/${sample_base}_1_unpaired.fastq.gz"
    r2_paired_out="${TRIMMED_READS_DIR}/${sample_base}_2_paired.fastq.gz"
    r2_unpaired_out="${TRIMMED_READS_DIR}/${sample_base}_2_unpaired.fastq.gz"

    # Trimmomatic command
    java -jar "$TRIMMOMATIC_JAR" PE \
        -threads "$THREADS" \
        -phred33 \
        "$r1_file" \
        "$r2_file" \
        "$r1_paired_out" \
        "$r1_unpaired_out" \
        "$r2_paired_out" \
        "$r2_unpaired_out" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2> "${LOG_DIR}/${sample_base}_trimmomatic.log"
done

echo "Trimmomatic processing complete. Trimmed files are in ${TRIMMED_READS_DIR}"
echo "------------------------------------"

echo "Pipeline finished successfully."


