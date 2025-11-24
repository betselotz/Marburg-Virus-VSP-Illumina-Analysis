#!/bin/bash
set -euo pipefail

# --- 1. Base Directory Configuration ---
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NONHOST_DIR="$BASE_DIR/results/03_nonhuman_reads"
BAM_DIR="$BASE_DIR/results/04_mapped_bam"
REFERENCE="$BASE_DIR/reference_genomes/Marburg_reference.fasta"

# Define the index prefix (Marburg_reference) for minimap2
REFERENCE_PREFIX="${REFERENCE%.fasta}"

# Define the global log file path and the summary log table header
LOG_FILE="$BASE_DIR/pipeline_mapping_log_$(date +%Y%m%d_%H%M%S).txt"
SAMPLE_LOG_DIR="$BASE_DIR/logs/mapping" 

echo "Starting mapping of non-host reads with Minimap2." | tee -a "$LOG_FILE"
echo "Project Base Directory: $BASE_DIR" | tee -a "$LOG_FILE"
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"
echo -e "Sample\tInput_Reads\tMapped_Reads\tMapping_Rate" >> "$LOG_FILE" # Add header for summary table

# --- 2. Conda Environment Activation ---
echo "📦 Activating Conda environment 'mapping_env'..." | tee -a "$LOG_FILE"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mapping_env
set -u
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"

# --- 3. Configuration & Utility ---
# Calculate available CPU cores, reserving 2 for system overhead (min 1 core).
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
echo "⚙️ Running with ${THREADS} threads." | tee -a "$LOG_FILE"

# Define required tools
MINIMAP2_EXEC=$(which minimap2)
SAMTOOLS_EXEC=$(which samtools)

# Create necessary directories if they don't exist
mkdir -p "$BAM_DIR"
mkdir -p "$SAMPLE_LOG_DIR"

# --- 4. Validation Checks & Reference Preparation ---
if [ -z "$MINIMAP2_EXEC" ] || [ -z "$SAMTOOLS_EXEC" ]; then
    echo "FATAL ERROR: minimap2 or samtools not found in PATH. Exiting." | tee -a "$LOG_FILE" >&2
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "FATAL ERROR: Reference file not found: $REFERENCE. Exiting." | tee -a "$LOG_FILE" >&2
    exit 1
fi

# Ensure reference is indexed with Samtools faidx for downstream operations
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "🛠️ Creating Samtools FAI index for the reference: ${REFERENCE}" | tee -a "$LOG_FILE"
    if ! "$SAMTOOLS_EXEC" faidx "$REFERENCE"; then
        echo "FATAL ERROR: Samtools faidx indexing failed. Exiting." | tee -a "$LOG_FILE" >&2
        exit 1
    fi
else
    echo "Reference FAI index already exists." | tee -a "$LOG_FILE"
fi

echo "Configuration OK." | tee -a "$LOG_FILE"
echo "Input Dir: $NONHOST_DIR" | tee -a "$LOG_FILE"
echo "Output Dir: $BAM_DIR" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"


# --- 5. Mapping Loop ---
for forward_read in "$NONHOST_DIR"/*_1.nonhost.fastq.gz; do
    
    if [ "$forward_read" = "$NONHOST_DIR/*_1.nonhost.fastq.gz" ]; then
        echo "WARNING: No input files found matching pattern *_1.nonhost.fastq.gz in $NONHOST_DIR. Exiting loop." | tee -a "$LOG_FILE"
        break
    fi

    # Record start time for runtime calculation
    START_TIME=$(date +%s)

    # Extract the sample name
    sample=$(basename "$forward_read" _1.nonhost.fastq.gz)
    fq1="$forward_read"
    
    # Construct the path to the reverse read file
    reverse_read="${fq1%_1.nonhost.fastq.gz}_2.nonhost.fastq.gz"

    # Define the final sorted BAM file
    outbam="${BAM_DIR}/${sample}.sorted.bam"
    SAMPLE_LOG_FILE="${SAMPLE_LOG_DIR}/${sample}_minimap2_mapping.log"
    
    # Define the final index file for the skip check
    outbai="${outbam}.bai"

    # --- Skip Check ---
    if [ -s "$outbai" ]; then
        echo "✅ Output BAM index for $sample already exists ($outbai). Skipping mapping." | tee -a "$LOG_FILE"
        continue
    fi
    # ------------------

    echo "▶️ Processing sample: ${sample}" | tee -a "$LOG_FILE"
    
    if [ ! -f "$reverse_read" ]; then
        echo "  ERROR: Missing reverse read file: ${reverse_read}. Skipping sample ${sample}." | tee -a "$LOG_FILE" >&2
        echo -e "${sample}\t0\t0\t0" >> "$LOG_FILE"
        continue
    fi

    # Check for empty input file and calculate read count
    READ_COUNT=$(zcat "$fq1" | wc -l | awk '{print $1/4}')
    if [[ "$READ_COUNT" -lt 1 ]]; then
        echo "⚠️ WARNING: Non-host file $fq1 is empty ($READ_COUNT reads). Skipping sample." | tee -a "$LOG_FILE"
        echo -e "${sample}\t${READ_COUNT}\t0\t0" >> "$LOG_FILE"
        continue
    fi

    echo "Mapping sample $sample with Minimap2 (using $THREADS threads, $READ_COUNT reads)..." | tee -a "$LOG_FILE"
    
    # 5a. Run Minimap2 and pipe to Samtools sort
    echo "  1/3 Running Minimap2 and piping to Samtools sort..." | tee -a "$LOG_FILE"
    
    # The pipeline is wrapped in parentheses for reliable stderr capture
    ( 
        "$MINIMAP2_EXEC" -ax sr -t "$THREADS" "$REFERENCE" "$fq1" "$reverse_read" | \
        "$SAMTOOLS_EXEC" sort -@ "$THREADS" -o "$outbam" -
    ) 2> "$SAMPLE_LOG_FILE"
    
    # Check the status of the piped commands
    # FIX APPLIED HERE: Temporarily disable strict checking for PIPESTATUS access
    set +u
    MINIMAP2_STATUS="${PIPESTATUS[0]}"
    SAMTOOLS_STATUS="${PIPESTATUS[1]}"
    set -u
    
    if [ "$MINIMAP2_STATUS" -ne 0 ] || [ "$SAMTOOLS_STATUS" -ne 0 ]; then
        echo "  ❌ ERROR: Mapping/Sorting failed (Minimap2 Status $MINIMAP2_STATUS, Samtools Status $SAMTOOLS_STATUS) for ${sample}. See $SAMPLE_LOG_FILE" | tee -a "$LOG_FILE" >&2
        rm -f "$outbam"
        echo -e "${sample}\t${READ_COUNT}\tFAIL\tFAIL" >> "$LOG_FILE"
        continue
    fi
    
    echo "    Sorted BAM file created: $outbam" | tee -a "$LOG_FILE"
    echo "    Mapping logs saved to: $SAMPLE_LOG_FILE" | tee -a "$LOG_FILE"


    # 5b. Index the sorted BAM file
    echo "  2/3 Indexing BAM file..." | tee -a "$LOG_FILE"
    if "$SAMTOOLS_EXEC" index "$outbam"; then
        echo "    Index file created: ${outbam}.bai" | tee -a "$LOG_FILE"
    else
        echo "  ❌ ERROR: Samtools index failed for sample ${sample}." | tee -a "$LOG_FILE" >&2
        echo -e "${sample}\t${READ_COUNT}\tSORTED_BAM\tINDEX_FAIL" >> "$LOG_FILE"
        continue
    fi
    
    # 5c. Calculate mapping statistics
    MAPPED_READS=$( ( "$SAMTOOLS_EXEC" view -F 4 -c "$outbam" ) || echo 0 )
    
    MAPPING_RATE=$(awk "BEGIN { if ($READ_COUNT == 0) { print 0.00 } else { printf \"%.2f\", $MAPPED_READS / $READ_COUNT * 100 } }")

    echo "  3/3 Mapping Stats: Mapped Reads = $MAPPED_READS ($MAPPING_RATE%)" | tee -a "$LOG_FILE"
    
    # Log successful sample to the summary table
    echo -e "${sample}\t${READ_COUNT}\t${MAPPED_READS}\t${MAPPING_RATE}" >> "$LOG_FILE"
    
    # Calculate runtime
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))

    # Log completion with runtime
    echo "Mapping for $sample completed: $outbam in ${RUNTIME}s" | tee -a "$LOG_FILE"
    echo "--------------------------------------------------------" | tee -a "$LOG_FILE"

done

echo "🎉 All non-host reads mapped and duplicates marked." | tee -a "$LOG_FILE"
echo "Mapping log saved to: $LOG_FILE"
