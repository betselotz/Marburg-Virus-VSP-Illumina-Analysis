#!/bin/bash
set -euo pipefail


# --- Configuration & Utility ---
# Calculate available CPU cores, reserving 2 (min 1 core).
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "⚙️ Running with ${THREADS} threads (Max available: $(nproc))"

# Display progress separator (User Requested)
echo "Processing... \------------------------------------------------------"
echo " STEP: Host Removal and Viral Mapping (3 Steps)"
echo "------------------------------------------------------"
echo "--- Host Removal and Viral Mapping: Started ---"


# Base directories
# PROJECT_DIR assumes the script is run from a subdirectory like 'scripts' (or similar structure).
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CLEAN_READS="$PROJECT_DIR/results/02_clean_reads"
# NON_HOST_READS_DIR is the final destination for the non-host reads
NON_HOST_READS_DIR="$PROJECT_DIR/results/03_nonhuman_reads"

# --- Reference and Index Configuration (ASSUMED TO EXIST) ---
# Assuming reference files and index are directly under reference_genomes/
HUMAN_REF_DIR="$PROJECT_DIR/reference_genomes" 

# FASTA filename (for reference, index must be built from this)
reference_genome="$HUMAN_REF_DIR/GCA_000001405.28_GRCh38.p13_genomic.fna"

# Define the index prefix variable (must match the existing index files, e.g., human_index.1.bt2)
index_prefix="$HUMAN_REF_DIR/human_index" 

# Activate Conda Environment
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate host_removal_env
set -u

# --- Prerequisites ---
# Requires: Bowtie2, samtools, gzip. The index files must be present at "$index_prefix.*.bt2".

### 3. Host Genome Removal with Bowtie2

# Check for required Bowtie2 index files
if [[ ! -f "${index_prefix}.1.bt2" ]]; then
  echo "Error: Bowtie2 index files not found at ${index_prefix}.*.bt2."
  echo "Please ensure the index is built and named 'human_index' in the reference_genomes directory."
  exit 1
fi
echo "✅ Bowtie2 index found. Skipping build/download."

# --- Output Directory Setup ---
# Create the main output directory (final destination for non-host reads)
mkdir -p "$NON_HOST_READS_DIR"


# --- Sample Initialization and Pre-Check ---

# Initialize counters and list for samples that need processing
PROCESSED_SAMPLES=0
SKIPPED_SAMPLES=0
FAILED_SAMPLES=0
SAMPLES_TO_PROCESS=()

# Check if CLEAN_READS directory contains any gzipped R1 files to process
shopt -s nullglob
r1_files=("$CLEAN_READS"/*_R1.trimmed.fastq.gz)
shopt -u nullglob

if [ ${#r1_files[@]} -eq 0 ]; then
    echo "No matching R1 files found in $CLEAN_READS. Expected pattern: *_R1.trimmed.fastq.gz"
else
    # 1. Pre-Check Loop: Determine which samples need processing and print skip messages
    for forward_read in "${r1_files[@]}"; do
        sample=$(basename "$forward_read" _R1.trimmed.fastq.gz)

        # Define final expected output paths using the new naming convention
        final_unmapped_R1="$NON_HOST_READS_DIR/${sample}_1.nonhost.fastq.gz"
        final_unmapped_R2="$NON_HOST_READS_DIR/${sample}_2.nonhost.fastq.gz"
        
        # Check if both final output files already exist
        if [[ -f "$final_unmapped_R1" ]] && [[ -f "$final_unmapped_R2" ]]; then
            echo "✅ nonhuman_reads output for ${sample} exists, skipping..."
            SKIPPED_SAMPLES=$((SKIPPED_SAMPLES + 1))
        else
            SAMPLES_TO_PROCESS+=("$forward_read")
        fi
    done
    
    # 2. Processing Loop: Only process samples identified in the list
    if [ ${#SAMPLES_TO_PROCESS[@]} -gt 0 ]; then
        echo "Starting read mapping to human host genome..."

        for forward_read in "${SAMPLES_TO_PROCESS[@]}"; do
          # Derive sample name and define paths again (safer scope)
          sample=$(basename "$forward_read" _R1.trimmed.fastq.gz)
          reverse_read="$CLEAN_READS/${sample}_R2.trimmed.fastq.gz"

          # Check if the corresponding R2 read exists (already checked R1 existence implicitly)
          if [[ ! -f "$reverse_read" ]]; then
            echo "Error: Reverse read $reverse_read not found for sample $sample. Skipping this sample."
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue
          fi

          # --- Output Definitions ---
          # Directory for alignment stats and BAM files
          SAMPLE_ALIGNMENT_DIR="$NON_HOST_READS_DIR/alignment_stats/${sample}"
          mkdir -p "$SAMPLE_ALIGNMENT_DIR"
          
          # Final compressed output files (NEW NAMING CONVENTION)
          bam_output="$SAMPLE_ALIGNMENT_DIR/${sample}.bam"
          flagstat_output="$SAMPLE_ALIGNMENT_DIR/${sample}.flagstat"
          final_unmapped_R1="$NON_HOST_READS_DIR/${sample}_1.nonhost.fastq.gz"
          final_unmapped_R2="$NON_HOST_READS_DIR/${sample}_2.nonhost.fastq.gz"
          
          # Temporary uncompressed files for unmapped reads
          temp_unmapped_prefix="$SAMPLE_ALIGNMENT_DIR/${sample}_temp_unmapped"

          # --- Sample Processing Header (User Requested) ---
          echo "🚀 Processing sample: ${sample}"
          
          # 1. Run Bowtie2 (Alignment to BAM, Unmapped to Temp Files)
          # Use process substitution (2> >(tee ...)) to separate stderr (logs) from stdout (SAM data).
          bowtie2 \
            --very-sensitive-local \
            --score-min L,0,-0.6 \
            --ma 0 \
            -x "$index_prefix" \
            -1 "$forward_read" \
            -2 "$reverse_read" \
            -p "$THREADS" \
            --un-conc "$temp_unmapped_prefix" \
            -S - \
            2> >(tee "$SAMPLE_ALIGNMENT_DIR/${sample}_bowtie2_log.txt" >&2) \
            | samtools view -@ "$THREADS" -bS - > "$bam_output"

          # Check for successful Bowtie2/Samtools run
          if [[ $? -ne 0 ]]; then
            echo "Error in Bowtie2/Samtools pipe for sample $sample. Review log file."
            rm -f "${temp_unmapped_prefix}.1" "${temp_unmapped_prefix}.2" || true # '|| true' prevents failure if files don't exist
            FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
            continue
          fi

          echo "Mapping completed. BAM file saved to ${bam_output}"
          
          # 2. Compress and move unmapped reads
          echo "Compressing non-host reads..."
          # Note: Bowtie2 names the un-conc files based on the prefix: .1 and .2
          gzip -c "${temp_unmapped_prefix}.1" > "$final_unmapped_R1"
          gzip -c "${temp_unmapped_prefix}.2" > "$final_unmapped_R2"
          
          # 3. Clean up temporary uncompressed files
          rm "${temp_unmapped_prefix}.1" "${temp_unmapped_prefix}.2"

          # 4. Print mapping statistics
          samtools flagstat -@ "$THREADS" "$bam_output" > "${flagstat_output}"
          echo "Mapping statistics saved to ${flagstat_output}"
          echo "Non-host reads saved to ${final_unmapped_R1} and ${final_unmapped_R2}"
          
          PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1))
        done
    fi
fi

# --- Final Summary (User Requested) ---
echo "--- Host removal: Completed ---"
echo "Summary:"
echo "Processed Successfully: $((PROCESSED_SAMPLES + SKIPPED_SAMPLES)) samples."
echo "Failed/Errored: ${FAILED_SAMPLES} samples."
echo "--- Host Removal and Viral Mapping: Finished ---"
echo "✅ Host removal process complete."
