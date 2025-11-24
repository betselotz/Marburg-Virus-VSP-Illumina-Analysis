#!/bin/bash
set -euo pipefail

##############################################
# Full MARG2 Pipeline
# Betselot Zerihun Ayano
# GitHub-ready version: runs scripts in scripts/
# Interactive terminal with spinners for long steps
# Preserves the results directory structure exactly
##############################################

# ===== Terminal colors =====
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

# ===== Base directories =====
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="$PROJECT_DIR/scripts"

FASTP_DIR="$PROJECT_DIR/results/01_fastp"
HOSTREM_DIR="$PROJECT_DIR/results/02_clean_reads"
NONHOST_DIR="$PROJECT_DIR/results/03_nonhuman_reads"
MAP_DIR="$PROJECT_DIR/results/04_mapped_bam"
VAR_DIR="$PROJECT_DIR/results/05_variants"
CON_DIR="$PROJECT_DIR/results/06_consensus"
COV_DIR="$PROJECT_DIR/results/07_coverage"
MSA_DIR="$PROJECT_DIR/results/08_msa"
IQTREE_DIR="$PROJECT_DIR/results/09_iqtree"
MULTIQC_DIR="$PROJECT_DIR/results/10_multiqc"

mkdir -p "$FASTP_DIR" "$HOSTREM_DIR" "$NONHOST_DIR" "$MAP_DIR" "$VAR_DIR" "$CON_DIR" "$COV_DIR" "$MSA_DIR" "$IQTREE_DIR" "$MULTIQC_DIR"

# ===== Spinner function =====
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while ps -p $pid > /dev/null 2>&1; do
        for i in $(seq 0 3); do
            printf "\r${CYAN}Processing... %s${NC}" "${spinstr:$i:1}"
            sleep $delay
        done
    done
    printf "\r"
}

# ===== Print section header =====
print_step() {
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW} $1 ${NC}"
    echo -e "${YELLOW}========================================${NC}"
}

echo -e "🚀 ${CYAN}Starting full MARG2 pipeline at $(date)${NC}\n"

########################################
# STEP 1: Fastp QC and trimming
########################################
print_step "STEP 1: Fastp QC and trimming"
bash "$SCRIPTS_DIR/fastp_batch.sh"
echo ""

########################################
# STEP 2: Host removal
########################################
print_step "STEP 2: Host removal"
bash "$SCRIPTS_DIR/host_removal_batch.sh" &
HOST_PID=$!
spinner $HOST_PID
wait $HOST_PID
echo -e "\r${GREEN}✅ Host removal completed.${NC}\n"

########################################
# STEP 3: Mapping non-host reads
########################################
print_step "STEP 3: Mapping non-host reads"
bash "$SCRIPTS_DIR/mapping_batch.sh"
echo -e "${GREEN}✅ All non-host reads mapped.${NC}\n"

########################################
# STEP 4: Variant calling (iVar)
########################################
print_step "STEP 4: Variant calling (iVar)"
bash "$SCRIPTS_DIR/variant_calling_batch.sh"
echo ""

########################################
# STEP 5: Consensus generation (iVar)
########################################
print_step "STEP 5: Consensus generation (iVar)"
bash "$SCRIPTS_DIR/consensus_batch.sh"
echo -e "${GREEN}✅ Consensus generation completed.${NC}\n"

########################################
# STEP 6: Coverage calculation
########################################
print_step "STEP 6: Coverage calculation"
bash "$SCRIPTS_DIR/coverage_batch.sh" &
COV_PID=$!
spinner $COV_PID
wait $COV_PID
echo -e "\r${GREEN}✅ Coverage calculation completed.${NC}\n"

########################################
# STEP 7: MultiQC report
########################################
print_step "STEP 9: MultiQC report"
bash "$SCRIPTS_DIR/multiqc_batch.sh"
echo -e "${GREEN}✅ MultiQC report completed.${NC}\n"

echo -e "🎉 ${MAGENTA}All steps completed successfully!${NC}"

