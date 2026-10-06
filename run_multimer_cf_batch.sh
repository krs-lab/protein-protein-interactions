#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MSA_DIR="$SCRIPT_DIR/msa_output"
OUTPUT_DIR="$SCRIPT_DIR/output_predictions"
LOG_FILE="$SCRIPT_DIR/colabfold_batch_multimer_biomix.log"

export PATH="/home/pritom/schmitzLab/main/pipBinders/colabfold/localcolabfold/.pixi/envs/default/bin:$PATH"

mkdir -p "$OUTPUT_DIR"

echo "==========================================" | tee -a "$LOG_FILE"
echo "Starting ColabFold AF2-Multimer Batch on BioMix: $(date)" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

# Run colabfold_batch directly on pre-computed MSAs
colabfold_batch \
    --model-type alphafold2_multimer_v3 \
    --num-models 2 \
    --rank iptm \
    "$MSA_DIR" \
    "$OUTPUT_DIR" >> "$LOG_FILE" 2>&1

echo "Processing outputs in $OUTPUT_DIR..." | tee -a "$LOG_FILE"

# Safely rename top-ranked PDB files to <Pair_ID>.pdb inside output/
for unrelaxed_pdb in "$OUTPUT_DIR"/*_unrelaxed_rank_001_*.pdb; do
    [ -f "$unrelaxed_pdb" ] || continue

    filename=$(basename "$unrelaxed_pdb")
    pair_id="${filename%%_unrelaxed*}"

    relaxed_pdb=$(ls "$OUTPUT_DIR"/${pair_id}_relaxed_rank_001_*.pdb 2>/dev/null | head -n 1 || true)

    if [ -n "$relaxed_pdb" ] && [ -f "$relaxed_pdb" ]; then
        target_pdb="$relaxed_pdb"
    else
        target_pdb="$unrelaxed_pdb"
    fi

    mv "$target_pdb" "$OUTPUT_DIR/${pair_id}.pdb"
    echo "Renamed $(basename "$target_pdb") -> ${pair_id}.pdb" >> "$LOG_FILE"
done

# Clean up rank 2-5 models to conserve storage
rm -f "$OUTPUT_DIR"/*_rank_00[2-5]_*.pdb

echo "Finished processing at $(date)" | tee -a "$LOG_FILE"
