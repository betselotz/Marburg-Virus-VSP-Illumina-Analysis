#!/bin/bash
set -euo pipefail

# ---------------------------------------
# CONFIGURATION & DIRECTORIES
# ---------------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
CONS_DIR="$PROJECT_DIR/results/06_consensus"
LOG_FILE="$PROJECT_DIR/pipeline_consensus_log_$(date +%Y%m%d_%H%M%S).txt"

# Create consensus directory if it doesn't exist
mkdir -p "$CONS_DIR"

# Determine number of threads (reserve 2 cores minimum 1)
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

echo "Starting consensus generation workflow." | tee -a "$LOG_FILE"
echo "⚙️ Running with $THREADS threads using samtools + ivar." | tee -a "$LOG_FILE"
echo "--------------------------------------------------------" | tee -a "$LOG_FILE"

# ---------------------------------------
# ACTIVATE CONDA ENVIRONMENT
# ---------------------------------------
echo "📦 Activating Conda environment 'ivar_env'..." | tee -a "$LOG_FILE"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ivar_env
set -u

SAMTOOLS_EXEC=$(which samtools)
IVAR_EXEC=$(which ivar)

# ---------------------------------------
# FIND BAM FILES
# ---------------------------------------
BAM_FILES=("$BAM_DIR"/*.sorted.bam)
if [[ ! -e "${BAM_FILES[0]}" ]]; then
    echo "❌ FATAL: No BAM files (*.sorted.bam) found in $BAM_DIR. Exiting." | tee -a "$LOG_FILE" >&2
    conda deactivate
    exit 1
fi

# ---------------------------------------
# LOOP THROUGH BAM FILES
# ---------------------------------------
for bam in "${BAM_FILES[@]}"; do
    START_TIME=$(date +%s)
    sample=$(basename "$bam" .sorted.bam)
    consensus_file="$CONS_DIR/${sample}.fa"

    # Skip if consensus already exists
    if [[ -s "$consensus_file" ]]; then
        echo "✅ Consensus for $sample exists, skipping..." | tee -a "$LOG_FILE"
        continue
    fi

    echo "----------------------------------------------------" | tee -a "$LOG_FILE"
    echo "🔹 Generating consensus for sample: $sample" | tee -a "$LOG_FILE"

    # Generate consensus (metagenomic-friendly thresholds)
    if ! (
        "$SAMTOOLS_EXEC" mpileup -A -Q 0 "$bam" |
        "$IVAR_EXEC" consensus -p "$CONS_DIR/$sample" -q 20 -t 0.7 -m 1
    ) 2>> "$LOG_FILE"; then
        echo "❌ ERROR: Consensus generation failed for $sample. Skipping to next sample." | tee -a "$LOG_FILE" >&2
        rm -f "$consensus_file"
        continue
    fi

    # Rename FASTA header
    if [[ -f "$CONS_DIR/${sample}.fa" ]]; then
        sed -i "1s/.*/>${sample} | Homo sapiens | Ethiopia | 2025/" "$CONS_DIR/${sample}.fa"
        echo "✅ Consensus completed for $sample: $CONS_DIR/${sample}.fa" | tee -a "$LOG_FILE"
    else
        echo "⚠️ Consensus file not found after processing $sample!" | tee -a "$LOG_FILE"
    fi

    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    echo "⏱ Runtime for $sample: ${RUNTIME}s" | tee -a "$LOG_FILE"
done

# ---------------------------------------
# CLEANUP
# ---------------------------------------
conda deactivate
echo "Conda environment deactivated." | tee -a "$LOG_FILE"

TOTAL_CONSENSUS=$(find "$CONS_DIR" -maxdepth 1 -name "*.fa" -type f | wc -l)
echo -e "\nConsensus generation completed." | tee -a "$LOG_FILE"
echo "Summary: Successfully generated ${TOTAL_CONSENSUS} consensus FASTAs in $CONS_DIR." | tee -a "$LOG_FILE"
echo "Process log saved to: $LOG_FILE" | tee -a "$LOG_FILE"

