#!/bin/bash
# Runs after the array finishes. Every checkpoint already exists, so
# model_selection.R only reads them, writes the summaries and plot, and
# refits the selected models on all ten years.

#SBATCH --job-name=weather_aggregate
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=greatlakes/logs/aggregate_%j.log

set -euo pipefail

module load R

export OMP_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"

Rscript model_selection.R
