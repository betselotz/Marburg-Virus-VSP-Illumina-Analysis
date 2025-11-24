#!/bin/bash
set -euo pipefail

# ========================================
# STEP 4: Variant calling (iVar Pipeline)
# ========================================

# Base directories
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
VAR_DIR="$PROJECT_DIR/results/05_variants"
MARBURG_REFERENCE="$PROJECT_DIR/reference_genomes/Marburg_reference.fasta"

# Number of threads
THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "Using $THREADS threads"
echo "----------------------------------------------------"

# Activate Conda environment
echo "📦 Activating Conda environment 'ivar_env'..."
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u
echo "Active conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Create variants output directory
mkdir -p "$VAR_DIR"

# Check required tools
SAMTOOLS_EXEC=$(which samtools)
IVAR_EXEC=$(which ivar)

if [[ -z "$SAMTOOLS_EXEC" || -z "$IVAR_EXEC" ]]; then
    echo "❌ FATAL: samtools or ivar not found in the environment. Exiting."
    exit 1
fi

# ========================================
# VARIANT CALLING: Skip markdup, use .sorted.bam
# ========================================
for sorted_bam in "$BAM_DIR"/*.sorted.bam; do

    if [[ ! -f "$sorted_bam" ]]; then
        echo "⚠️ WARNING: No .sorted.bam files found in $BAM_DIR. Skipping."
        continue
    fi

    sample=$(basename "$sorted_bam" .sorted.bam)
    output_prefix="$VAR_DIR/${sample}_variants"
    output_file="${output_prefix}.tsv"

    # --- Skip if output already exists ---
    if [[ -s "$output_file" ]]; then
        echo "✅ Variant file already exists for $sample: $output_file. Skipping."
        continue
    fi

    echo "----------------------------------------------------"
    echo "🔹 Calling variants for $sample"
    echo "Input BAM: $sorted_bam"
    echo "Output prefix: $output_prefix"

    # Run iVar variant calling
    if ! "$SAMTOOLS_EXEC" mpileup -A -d 1000000 -B -Q 0 -f "$MARBURG_REFERENCE" "$sorted_bam" | \
        "$IVAR_EXEC" variants -r "$MARBURG_REFERENCE" -p "$output_prefix"; then
        echo "❌ ERROR: Variant calling failed for $sample. Skipping."
        continue
    fi

    echo "✅ Variant calling completed for $sample. Output: $output_file"

done

# Deactivate environment
conda deactivate
echo "Conda environment deactivated."
echo "🎉 Variant calling complete."

