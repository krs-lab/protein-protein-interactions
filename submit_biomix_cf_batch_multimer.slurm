#!/bin/bash
#SBATCH --job-name=cf_batch_multimer
#SBATCH --partition=gpu 
#SBATCH --qos=gpu   
#SBATCH --nodes=1
#SBATCH --gres=gpu:1  
#SBATCH --nodelist=biomix62,biomix63,biomix64,biomix66 
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=slurm_batch_%j.log
#SBATCH --error=slurm_batch_%j.err

set -e

chmod +x ./run_multimer_cf_batch.sh
./run_multimer_cf_batch.sh
