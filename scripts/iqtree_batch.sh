#!/bin/bash
set -euo pipefail

# ========================================
# STEP: Phylogenetic Tree Generation (IQ-TREE)
# CORRECTED: Uses hard copy instead of symbolic link for better file system compatibility.
# ========================================

# --- Configuration & Utility ---
export THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))
IQ_PREFIX="marv_phylogeny"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Directories
MSA_DIR="$PROJECT_DIR/results/09_msa"              
PHYLO_DIR="$PROJECT_DIR/results/10_phylogeny"      
INTERMEDIATES_DIR="$PHYLO_DIR/intermediates"     

# Input/Output Files
ALN_TRIM="$MSA_DIR/all_sequences_aligned_trimmed.fasta"
# Define the temporary input link path (used for IQ-TREE to name files)
TEMP_INPUT_COPY="$PHYLO_DIR/$IQ_PREFIX.fasta" # Renamed variable for clarity
TREE_FILE="$PHYLO_DIR/${IQ_PREFIX}.treefile"       
SUMMARY_FILE="$PHYLO_DIR/iqtree_summary.txt"      
LOG_FILE="$PHYLO_DIR/${IQ_PREFIX}_$(date +%Y%m%d_%H%M%S).log" 

mkdir -p "$PHYLO_DIR" "$INTERMEDIATES_DIR"

echo "⚙️ Running IQ-TREE (v3.0.1) with ${THREADS} threads." | tee "$LOG_FILE"
echo -e "\n========================================" | tee -a "$LOG_FILE"
echo "      STEP: Phylogenetic Tree Generation (IQ-TREE)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

# --- Global Cleanup TRAP ---
cleanup_on_error() {
    echo "❌ IQ-TREE failed. Intermediate files kept for debugging. Review $LOG_FILE for details." | tee -a "$LOG_FILE" >&2
    rm -f "$TEMP_INPUT_COPY" 2>/dev/null || true # Ensure temp copy is cleaned
}
trap cleanup_on_error ERR
# ---------------------------

# --- 1. Check Input File ---
if [[ ! -s "$ALN_TRIM" ]]; then
    echo "❌ FATAL: Trimmed MSA file not found or is empty at $ALN_TRIM. Exiting." | tee -a "$LOG_FILE"
    exit 1
fi

# --- 2. Check for Final Output ---
if [[ -s "$TREE_FILE" ]]; then
    echo "✅ Final tree file already exists at $TREE_FILE, skipping IQ-TREE." | tee -a "$LOG_FILE"
    exit 0
fi

# --- 3. Run IQ-TREE (Tool Check) ---
echo "--------------------------------------------" | tee -a "$LOG_FILE"
IQTREE_EXEC=$(which iqtree3)
if [[ -z "$IQTREE_EXEC" ]]; then
    echo "❌ FATAL: iqtree3 not found in PATH. Exiting." | tee -a "$LOG_FILE"
    exit 1
fi
echo "Tool check: iqtree3 found at $IQTREE_EXEC" | tee -a "$LOG_FILE"
echo "--------------------------------------------" | tee -a "$LOG_FILE"

# --- PREP: Copy input file into the target directory to control output naming ---
echo "Preparing input file copy: $TEMP_INPUT_COPY" | tee -a "$LOG_FILE"

# Use 'cp -f' to make a hard copy, avoiding the symbolic link error.
if ! cp -f "$ALN_TRIM" "$TEMP_INPUT_COPY"; then
    echo "❌ FATAL: Failed to copy input file. Check file system write permissions or disk space." | tee -a "$LOG_FILE"
    exit 1
fi

echo "🔹 Running IQ-TREE (Best-fit model selection + 1000 ultrafast bootstraps)..." | tee -a "$LOG_FILE"

START=$(date +%s)

# Execute IQ-TREE. Pass the copied file path as the input file (-s).
set -o pipefail
"$IQTREE_EXEC" \
    -s "$TEMP_INPUT_COPY" \
    -m MFP \
    -bb 1000 \
    -nt "$THREADS" 2>&1 | tee -a "$LOG_FILE"
set +o pipefail

END=$(date +%s)
RUNTIME=$((END-START))

# --- POST-RUN CLEANUP AND SUMMARY ---
rm -f "$TEMP_INPUT_COPY" # Remove the temporary hard copy

if [[ -s "$TREE_FILE" ]]; then
    
    echo -e "\n--------------------------------------------" | tee -a "$LOG_FILE"
    echo "✅ IQ-TREE run completed successfully." | tee -a "$LOG_FILE"
    echo "Final tree saved: $TREE_FILE" | tee -a "$LOG_FILE"
    echo "IQ-TREE runtime: ${RUNTIME} seconds" | tee -a "$LOG_FILE"

    # --- REFINED CLEANUP: Move useful files, delete true intermediates ---
    echo "Processing and moving intermediate files to $INTERMEDIATES_DIR..." | tee -a "$LOG_FILE"
    
    # 4a. Move ALL IQ-TREE outputs (*.prefix*) to the intermediate folder for review
    find "$PHYLO_DIR" -maxdepth 1 -type f -name "$IQ_PREFIX.*" \
        ! -name "$IQ_PREFIX.treefile" \
        ! -name "$IQ_PREFIX.log" \
        -exec mv {} "$INTERMEDIATES_DIR/" \; 2>/dev/null || true
    
    # 4b. SUMMARY TABLE GENERATION
    if [[ ! -f "$SUMMARY_FILE" ]]; then
        echo -e "Treefile\tRuntime_seconds\tBest_Model" > "$SUMMARY_FILE"
    fi
    
    # Extract best-fit model from the IQ-TREE log file
    BEST_MODEL=$(grep "ModelFinder best model:" "$INTERMEDIATES_DIR/$IQ_PREFIX.log" 2>/dev/null | awk '{print $NF}' || echo "N/A")
    
    echo -e "$TREE_FILE\t$RUNTIME\t$BEST_MODEL" >> "$SUMMARY_FILE"
    echo "Summary appended to $SUMMARY_FILE" | tee -a "$LOG_FILE"

else
    echo "❌ ERROR: IQ-TREE failed to generate the final treefile ($TREE_FILE)." | tee -a "$LOG_FILE"
    echo "Intermediate files are retained for manual inspection in $PHYLO_DIR." | tee -a "$LOG_FILE"
    exit 1
fi
