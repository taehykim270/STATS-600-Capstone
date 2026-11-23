#!/bin/bash
# One array task per (outcome, held-out year, seasonal structure, h) checkpoint.
# Submitted by greatlakes/submit.sh, which supplies --account and --array.

#SBATCH --job-name=weather_cv
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=greatlakes/logs/cv_%A_%a.log

set -euo pipefail

module load R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"

Rscript greatlakes/cv_task.R "$SLURM_ARRAY_TASK_ID"
