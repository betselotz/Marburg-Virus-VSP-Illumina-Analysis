#!/bin/bash
set -euo pipefail

# ========================================
# MSA Pipeline for Marburg Virus
# Corrected for variable passing using shell arrays.
# ========================================

# --- Configuration & Utility ---
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
LOG_FILE_BASE="$(date +%Y%m%d_%H%M%S)_mafft"
echo "⚙️ Running with ${THREADS} threads."
# -------------------------------

# Base project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Directories
CONS_DIR="$PROJECT_DIR/results/06_consensus"
REF_DIR="$PROJECT_DIR/reference_genomes/MARV_compare"
MSA_DIR="$PROJECT_DIR/results/09_msa"
mkdir -p "$MSA_DIR"

# Outgroup (fixed absolute path)
OUTGROUP="$PROJECT_DIR/reference_genomes/EF446131.1.fasta"

# Log file
LOG_FILE="$MSA_DIR/${LOG_FILE_BASE}.log"

# Output files
ALL_SEQ_INPUT="$MSA_DIR/all_sequences_input.fasta"
ALN_RAW="$MSA_DIR/all_sequences_raw_aligned.fasta"
ALN_TRIM="$MSA_DIR/all_sequences_aligned_trimmed.fasta" # Final Output

# Toggle for thorough alignment: "true" = SLOW (globalpair), "false" = FAST (recommended)
THOROUGH_ALIGNMENT=false

# --- Global Cleanup TRAP ---
cleanup() {
    echo "🧹 Running emergency cleanup..." | tee -a "$LOG_FILE"
    rm -f "$ALL_SEQ_INPUT" "$ALN_RAW"
}
trap cleanup EXIT

# Skip if final trimmed MSA already exists
if [[ -s "$ALN_TRIM" ]]; then
    echo "Trimmed MSA already exists, skipping..." | tee -a "$LOG_FILE"
    exit 0
fi

echo "--- MSA RUN STARTED ---" | tee "$LOG_FILE"
echo "THREADS: $THREADS" | tee -a "$LOG_FILE"
echo "Consensus: $CONS_DIR, Reference: $REF_DIR, Outgroup: $OUTGROUP" | tee -a "$LOG_FILE"

# --- Validate Input Files ---
CONS_INPUT_FILES=("$CONS_DIR/M261M.fa" "$CONS_DIR/M262M.fa")
for f in "${CONS_INPUT_FILES[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ FATAL: Consensus file not found: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
done

REF_FILES=("$REF_DIR"/*.fasta)
if [[ ${#REF_FILES[@]} -eq 0 ]]; then
    echo "❌ FATAL: No reference FASTAs found in $REF_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

if [[ ! -f "$OUTGROUP" ]]; then
    echo "❌ FATAL: Outgroup file not found at $OUTGROUP" | tee -a "$LOG_FILE"
    exit 1
fi

echo "✅ All input files validated." | tee -a "$LOG_FILE"

# --- Combine sequences ---
echo "Combining consensus, reference, and outgroup sequences into single input file..." | tee -a "$LOG_FILE"
cat "${CONS_INPUT_FILES[@]}" "${REF_FILES[@]}" "$OUTGROUP" > "$ALL_SEQ_INPUT"

# --- REQUIRED CHANGE: Sanitize FASTA Headers ---
echo "Sanitizing FASTA headers: Replacing spaces with underscores for single-token IDs..." | tee -a "$LOG_FILE"
# This ensures the entire descriptive line is treated as the sequence ID by MAFFT and downstream tools.
# It prevents truncation where tools only read the text before the first space.
# We use a temporary file extension (.bak) for sed's in-place edit, then remove it.
sed -i.bak '/^>/s/ /_/g' "$ALL_SEQ_INPUT"
rm -f "$ALL_SEQ_INPUT.bak"

# --- Activate MAFFT environment ---
echo "Activating MAFFT environment..." | tee -a "$LOG_FILE"
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mafft_env
set -u

# --- Run MAFFT Alignment (Optimized for speed) ---
echo "Running MAFFT alignment..." | tee -a "$LOG_FILE"
if [[ "$THOROUGH_ALIGNMENT" = true ]]; then
    # SLOWEST, highest-accuracy setting (For short sequences)
    echo "-> Using thorough alignment (--globalpair --maxiterate 1000)..." | tee -a "$LOG_FILE"
    MAFFT_OPTS=(--globalpair --maxiterate 1000) # Defined as array
else
    # RECOMMENDED FAST SETTING for long genomes (>10kb)
    echo "-> Using fast, standard genomic alignment (--retree 2 --maxiterate 2)..." | tee -a "$LOG_FILE"
    MAFFT_OPTS=(--retree 2 --maxiterate 2) # Defined as array
fi

# Execute MAFFT with corrected argument passing
# Note: ${MAFFT_OPTS[@]} must be double-quoted to handle spaces/special characters correctly
if ! mafft "${MAFFT_OPTS[@]}" --adjustdirection --thread "$THREADS" \
      "$ALL_SEQ_INPUT" > "$ALN_RAW" 2> >(tee -a "$LOG_FILE" >&2); then

    echo "❌ FATAL: MAFFT alignment failed! Check the error messages above/in the log file." | tee -a "$LOG_FILE" >&2
    exit 1
fi

# --- Trim alignment ---
echo "Trimming terminal gaps/poorly aligned regions with trimAl (-automated1)..." | tee -a "$LOG_FILE"
if ! trimal -in "$ALN_RAW" -out "$ALN_TRIM" -automated1 2>> "$LOG_FILE"; then
    echo "❌ FATAL: trimAl failed! Check $LOG_FILE for details." | tee -a "$LOG_FILE" >&2
    exit 1
fi

# --- Deactivate environment ---
conda deactivate

# --- Cleanup ---
echo "Cleaning up intermediate alignment files..." | tee -a "$LOG_FILE"
rm -f "$ALL_SEQ_INPUT" "$ALN_RAW"

echo "✅ MSA completed, outgroup included, trimmed. Final output: $ALN_TRIM" | tee -a "$LOG_FILE"
