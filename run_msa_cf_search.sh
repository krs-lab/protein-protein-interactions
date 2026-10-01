#!/bin/bash
set -euo pipefail # pipeline fails fi any cmd fails 

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUT_DIR="$SCRIPT_DIR/pair_fastas"
OUTPUT_DIR="$SCRIPT_DIR/msa_output"
LOG_FILE="$SCRIPT_DIR/colabfold_search_biomix.log"

export PATH="/home/pritom/schmitzLab/main/pipBinders/colabfold/localcolabfold/.pixi/envs/default/bin:$PATH"
export MMSEQS_IGNORE_INDEX=1

mkdir -p "$OUTPUT_DIR"

echo "==========================================" | tee -a "$LOG_FILE"
echo "Starting ColabFold MSA Search on BioMix at $(date)" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

# Count input pair FASTA files without splitting them
count=$(ls -1 "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fas 2>/dev/null | wc -l || true)

if [ "$count" -eq 0 ]; then
    echo "ERROR: No interaction FASTA files found in $INPUT_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Found $count interaction FASTA file(s) in $INPUT_DIR. Running colabfold_search..." | tee -a "$LOG_FILE"

# Run colabfold_search on intact pair FASTA files against the local BioMix DB
colabfold_search \
	--threads "${SLURM_CPUS_PER_TASK:-32}" \
    "$INPUT_DIR" \
    /mnt/dbases/colabfold_db \
    "$OUTPUT_DIR" >> "$LOG_FILE" 2>&1

echo "MSA generation complete. Outputs saved in $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Finished processing at $(date)" | tee -a "$LOG_FILE"
