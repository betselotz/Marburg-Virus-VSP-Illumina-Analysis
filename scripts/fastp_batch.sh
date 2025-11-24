#!/bin/bash
set -euo pipefail

# --- Configuration ---
MIN_QUAL=20       # Minimum quality score (Q) for trimming
MIN_LEN=50        # Minimum length (L) for retained reads
MAX_THREADS=8     # Absolute maximum threads to use for fastp (prevents resource hogging)
SIZE_THRESHOLD=0.90 # R1/R2 size must be >= 90% of the other file's size
PROCESSED_SAMPLES=0
FAILED_SAMPLES=0    

# --- Conda Environment Activation ---
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate fastp_env || { echo "❌ ERROR: Conda environment 'fastp_env' not found."; exit 1; }

# Safety check: make sure fastp is available
if ! command -v fastp &> /dev/null; then
    echo "❌ ERROR: fastp not found. Make sure it's installed in the fastp_env environment."
    exit 1
fi

# Determine threads dynamically, capping at MAX_THREADS
if command -v nproc &> /dev/null; then
    AVAILABLE_THREADS=$(nproc)
else
    AVAILABLE_THREADS=2 # Safe fallback
fi
THREADS=$(( AVAILABLE_THREADS > MAX_THREADS ? MAX_THREADS : AVAILABLE_THREADS ))
THREADS=$(( THREADS < 1 ? 1 : THREADS )) # Ensure at least 1 thread
echo "⚙️ Using $THREADS threads for fastp (Max set to $MAX_THREADS)."
echo "--- Marburg fastp QC: Started ---"

# --- Directory Setup (Sticking to your naming convention) ---
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

RAW_DIR="$PROJECT_DIR/raw_reads"
RESULTS_DIR="$PROJECT_DIR/results"
FASTP_QC_DIR="$RESULTS_DIR/01_fastp"
CLEAN_FASTQ_DIR="$RESULTS_DIR/02_clean_reads"

mkdir -p "$FASTP_QC_DIR" "$CLEAN_FASTQ_DIR"

# --- Main Processing Loop ---
while read -r fq1; do
    
    # 1. Robust Sample Name Extraction
    sample=$(basename "$fq1" | sed -E 's/(_R1|_R2).*\.f(ast)?q(\.gz)?$//')

    # 2. Safer R2 Filename Construction using Parameter Expansion
    fq2="${fq1/_R1/_R2}"

    out1="$CLEAN_FASTQ_DIR/${sample}_R1.trimmed.fastq.gz"
    out2="$CLEAN_FASTQ_DIR/${sample}_R2.trimmed.fastq.gz"
    html="$FASTP_QC_DIR/${sample}_fastp.html"
    json="$FASTP_QC_DIR/${sample}_fastp.json"
    log="$FASTP_QC_DIR/${sample}_fastp.log"
    html_title="fastp QC $sample"
    
    # Check if R2 file exists
    if [[ ! -f "$fq2" ]]; then
        echo "⚠️ WARNING: Paired-end file $fq2 not found for sample $sample. Skipping."
        continue
    fi
    
    # 3. Skip Empty Input Files (using -s for non-empty check)
    if [[ ! -s "$fq1" || ! -s "$fq2" ]]; then
        echo "⚠️ WARNING: Input file(s) for $sample is/are empty. Skipping."
        continue
    fi
    
    # 4. Critical: Check for significant file size discrepancy (Truncation/Corruption Check)
    fq1_size=$(stat -c%s "$fq1")
    fq2_size=$(stat -c%s "$fq2")

    if ! awk "BEGIN { exit ( $fq1_size < $fq2_size * $SIZE_THRESHOLD || $fq2_size < $fq1_size * $SIZE_THRESHOLD ) }" ; then
        echo "⚠️ WARNING: R1 and R2 files for $sample have a major size discrepancy."
        echo "    R1 size: $fq1_size bytes, R2 size: $fq2_size bytes. Continuing, but investigate this sample!"
    fi

    # Skip if final outputs exist
    if [[ -s "$out1" && -s "$out2" ]]; then
        echo "✅ fastp output for $sample exists, skipping..."
        PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1))
        continue
    fi

    echo "⚙️ Processing $sample (Q>$MIN_QUAL, L>$MIN_LEN)..."
    
    # Execute fastp with Graceful Error Handling and Optimized Params
    (
        fastp -i "$fq1" -I "$fq2" \
            -o "$out1" -O "$out2" \
            -h "$html" -j "$json" \
            --report_title "$html_title" \
            --thread "$THREADS" \
            --detect_adapter_for_pe \
            --length_required "$MIN_LEN" \
            --qualified_quality_phred "$MIN_QUAL" \
            --trim_poly_g \
            --trim_poly_x \
            --low_complexity_filter \
            --cut_window_size 4 --cut_mean_quality 20 \
            --trim_front1 0 --trim_front2 0
    ) 2> "$log" && PROCESSED_SAMPLES=$((PROCESSED_SAMPLES + 1)) || { 
        echo "❌ ERROR: fastp failed for $sample. See log file $log for details."
        FAILED_SAMPLES=$((FAILED_SAMPLES + 1))
    }

    echo "Finished $sample."
done < <(find "$RAW_DIR" -maxdepth 1 -name "*_R1*.fastq*")

# --- Summary Report ---
echo "--- Marburg fastp QC: Completed ---"
echo "Summary:"
echo "Processed Successfully: $PROCESSED_SAMPLES samples."
echo "Failed/Errored: $FAILED_SAMPLES samples."
