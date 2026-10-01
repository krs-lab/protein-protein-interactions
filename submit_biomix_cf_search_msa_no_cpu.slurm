#!/bin/bash
#SBATCH --job-name=cf_search_msa
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

set -e # exit immediately if any cmd fails
# set -x # for line-by-line execution tracing 

chmod +x ./run_msa_cf_search.sh
./run_msa_cf_search.sh
