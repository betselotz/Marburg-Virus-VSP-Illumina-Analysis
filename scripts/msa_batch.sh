#!/bin/bash
set -euo pipefail

# ========================================
# High-Confidence, Memory-Safe MSA Pipeline for Marburg Virus
# ========================================

# -------- CPU AUTO-DETECTION (Safe Mode) --------
TOTAL_THREADS=$(nproc)
# Use all-minus-two for normal tasks, but constrain MAFFT for memory safety
SAFE_THREADS=$(( TOTAL_THREADS > 4 ? TOTAL_THREADS - 2 : 2 ))
echo "🧠 System Threads: $TOTAL_THREADS  |  Safe Threads: $SAFE_THREADS"

# --- Project Directories ---
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONS_DIR="$PROJECT_DIR/results/06_consensus_selected"     # <=== UPDATED
REF_DIR="$PROJECT_DIR/reference_genomes/MARV_compare"
MSA_DIR="$PROJECT_DIR/results/09_msa"
OUTGROUP="$PROJECT_DIR/reference_genomes/EF446131.1.fasta"

mkdir -p "$MSA_DIR"

LOG_FILE="$MSA_DIR/$(date +%Y%m%d_%H%M%S)_mafft.log"

ALL_SEQ_INPUT="$MSA_DIR/all_sequences_input.fasta"
ALN_RAW="$MSA_DIR/all_sequences_raw_aligned.fasta"
ALN_TRIM="$MSA_DIR/all_sequences_aligned_trimmed.fasta"

echo "⚙️ Starting MSA Pipeline" | tee "$LOG_FILE"

# -------- Emergency Cleanup --------
cleanup() {
    echo "🧹 Cleanup triggered..." | tee -a "$LOG_FILE"
    rm -f "$ALL_SEQ_INPUT" "$ALN_RAW"
}
trap cleanup EXIT

# -------- Skip if final MSA exists --------
if [[ -s "$ALN_TRIM" ]]; then
    echo "MSA already exists. Exiting." | tee -a "$LOG_FILE"
    exit 0
fi

# -------- Validate Inputs --------
# UPDATED: Load ALL .fa files from 06_consensus_selected
CONS_INPUT_FILES=("$CONS_DIR"/*.fa)

for f in "${CONS_INPUT_FILES[@]}" "$REF_DIR"/*.fasta "$OUTGROUP"; do
    if [[ ! -s "$f" ]]; then
        echo "❌ FATAL: Missing or empty FASTA: $f" | tee -a "$LOG_FILE"
        exit 1
    fi

    SEQ_LEN=$(grep -v ">" "$f" | tr -d '\n' | wc -c)
    if [[ $SEQ_LEN -lt 10000 ]]; then
        echo "❌ FATAL: Sequence too short: $f (${SEQ_LEN} bp)" | tee -a "$LOG_FILE"
        exit 1
    fi
done

echo "✅ All sequences validated." | tee -a "$LOG_FILE"

# -------- Combine Sequences --------
cat "${CONS_INPUT_FILES[@]}" "$REF_DIR"/*.fasta "$OUTGROUP" > "$ALL_SEQ_INPUT"

# -------- Ensure Header Uniqueness --------
awk '/^>/{if(++c[$0]>1) $0=$0"_"c[$0]}1' "$ALL_SEQ_INPUT" > "$ALL_SEQ_INPUT.tmp"
mv "$ALL_SEQ_INPUT.tmp" "$ALL_SEQ_INPUT"
echo "✅ Headers unique." | tee -a "$LOG_FILE"

# -------- Activate MAFFT env --------
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mafft_env
set -u

# ========================================
#        MAFFT HIGH-CONFIDENCE MODE
# ========================================

# High accuracy, but memory safe
MAFFT_OPTS=(--globalpair --maxiterate 1000)
MAFFT_THREADS=2     # <== CRITICAL: prevents RAM crash on 32 GB

echo "🔥 Running MAFFT High-Accuracy Mode" | tee -a "$LOG_FILE"
echo "   Settings: --globalpair --maxiterate 1000" | tee -a "$LOG_FILE"
echo "   Threads limited to: $MAFFT_THREADS (safe mode)" | tee -a "$LOG_FILE"

if ! mafft "${MAFFT_OPTS[@]}" --adjustdirection --thread "$MAFFT_THREADS" \
      "$ALL_SEQ_INPUT" > "$ALN_RAW" 2>>"$LOG_FILE"; then
    echo "❌ MAFFT FAILED in high-accuracy mode" | tee -a "$LOG_FILE"
    echo "⏳ Retrying with SAFER mode (--localpair --maxiterate 1000)" | tee -a "$LOG_FILE"

    MAFFT_OPTS=(--localpair --maxiterate 1000)

    if ! mafft "${MAFFT_OPTS[@]}" --adjustdirection --thread "$SAFE_THREADS" \
          "$ALL_SEQ_INPUT" > "$ALN_RAW" 2>>"$LOG_FILE"; then
        echo "❌ FATAL: MAFFT failed in fallback mode." | tee -a "$LOG_FILE"
        exit 1
    fi
fi

echo "✅ MAFFT alignment complete." | tee -a "$LOG_FILE"

# -------- TRIMAL --------
echo "✂️ Trimming with trimAl..." | tee -a "$LOG_FILE"
if ! trimal -in "$ALN_RAW" -out "$ALN_TRIM" -automated1 2>>"$LOG_FILE"; then
    echo "❌ FATAL: trimAl failed" | tee -a "$LOG_FILE"
    exit 1
fi

# -------- Cleanup --------
conda deactivate
rm -f "$ALL_SEQ_INPUT" "$ALN_RAW"

echo "🎉 COMPLETE: High-confidence MSA generated."
echo "👉 Output: $ALN_TRIM" | tee -a "$LOG_FILE"

