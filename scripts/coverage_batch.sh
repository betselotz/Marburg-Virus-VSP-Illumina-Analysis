#!/bin/bash
set -euo pipefail

#######################################
# MARBURG DEPTH & COVERAGE ANALYSIS
# Updated version: robust, multi-threaded, improved logging
#######################################

# ---------------------------
# CONFIGURATION & DIRECTORIES
# ---------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_DIR="$PROJECT_DIR/results/04_mapped_bam"
CONS_DIR="$PROJECT_DIR/results/06_consensus"
COVERAGE_DIR="$PROJECT_DIR/results/07_coverage"
mkdir -p "$COVERAGE_DIR"

# Threads
THREADS=$(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))

echo "Starting coverage calculation..."
echo "Using $THREADS threads for samtools depth"

# ---------------------------
# Logging helpers
# ---------------------------
log() {
    echo -e "$(date '+%Y-%m-%d %H:%M:%S') | $*"
}
log_sample() {
    echo -e "  → $*"
}

# ---------------------------
# Count Ns in FASTA
# ---------------------------
count_N_bases() {
    local fasta_file="$1"
    if [[ ! -f "$fasta_file" ]]; then
        echo "0 0 0.00"
        return
    fi
    local total_len n_count percent_n
    total_len=$(grep -v "^>" "$fasta_file" | tr -d '\n' | wc -c)
    n_count=$(grep -v "^>" "$fasta_file" | tr -d '\n' | tr 'a-z' 'A-Z' | tr -cd 'N' | wc -c)
    percent_n=$(awk -v n="$n_count" -v len="$total_len" 'BEGIN{printf "%.2f", (n/len)*100}')
    echo "$n_count $total_len $percent_n"
}

# ---------------------------
# Coverage statistics
# ---------------------------
calculate_coverage_stats() {
    local depth_file="$1"
    local glen="$2"
    awk -v glen="$glen" '
    BEGIN {
        min=999999; max=0; sum=0; sum_sq=0; count=0;
        covered1=0; covered10=0; covered20=0; covered30=0; covered50=0;
    }
    {
        d=$3;
        sum+=d; sum_sq+=d*d; count++;
        if(d<min) min=d;
        if(d>max) max=d;
        if(d>=1) covered1++;
        if(d>=10) covered10++;
        if(d>=20) covered20++;
        if(d>=30) covered30++;
        if(d>=50) covered50++;
    }
    END {
        mean=(count>0)?sum/count:0;
        stddev=(count>1)?sqrt((sum_sq-(sum*sum)/count)/(count-1)):0;
        median=mean; # approximate
        cov1=(covered1/glen)*100;
        cov10=(covered10/glen)*100;
        cov20=(covered20/glen)*100;
        cov30=(covered30/glen)*100;
        cov50=(covered50/glen)*100;
        printf("%.2f\t%.0f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
               mean, median, min, max, stddev, cov1, cov10, cov20, cov30, cov50);
    }' "$depth_file"
}

# ---------------------------
# MAIN
# ---------------------------
main() {
    summary_file="$COVERAGE_DIR/depth_summary.tsv"
    echo -e "Sample\tMeanDepth\tMedianDepth\tMinDepth\tMaxDepth\tStdDevDepth\tCoverage1x\tCoverage10x\tCoverage20x\tCoverage30x\tCoverage50x\tN_count\tGenome_length\tPercent_N" > "$summary_file"

    for bam in "$BAM_DIR"/*.sorted.bam; do
        sample=$(basename "$bam" .sorted.bam)
        depth_file="$COVERAGE_DIR/${sample}.depth"
        fasta_file="$CONS_DIR/${sample}.fa"

        # Skip if already done
        if [[ -s "$depth_file" && -s "$summary_file" && $(grep -c "^${sample}\t" "$summary_file") -gt 0 ]]; then
            log_sample "Coverage already computed for $sample. Skipping."
            continue
        fi

        if [[ ! -f "$bam" ]]; then
            log_sample "BAM not found: $bam, skipping."
            continue
        fi
        if [[ ! -f "$fasta_file" ]]; then
            log_sample "Consensus FASTA not found for $sample, skipping."
            continue
        fi

        log_sample "Processing $sample..."

        # Generate depth file
        samtools depth -a -@ "$THREADS" "$bam" > "$depth_file"

        # Compute stats
        genome_len=$(grep -v "^>" "$fasta_file" | tr -d '\n' | wc -c)
        stats=$(calculate_coverage_stats "$depth_file" "$genome_len")
        read n_count genome_len percent_n <<< $(count_N_bases "$fasta_file")

        # Append to summary
        echo -e "${sample}\t${stats}\t${n_count}\t${genome_len}\t${percent_n}" >> "$summary_file"

        # Quick summary log
        mean_depth=$(echo "$stats" | cut -f1)
        cov1=$(echo "$stats" | cut -f6)
        cov30=$(echo "$stats" | cut -f9)
        log_sample "Done: Mean=${mean_depth}x, >=1x=${cov1}%, >=30x=${cov30}%, N_count=${n_count}, Percent_N=${percent_n}%"
    done

    echo -e "\n✅ Coverage summary generated at $summary_file"
}

# ---------------------------
# Check dependencies
# ---------------------------
for cmd in samtools awk grep tr wc; do
    command -v $cmd >/dev/null 2>&1 || { echo "$cmd not found, please install"; exit 1; }
done

# Run main
main

