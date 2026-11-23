#!/bin/bash
# Run ON GREAT LAKES, from the project folder:
#   bash greatlakes/submit.sh [slurm_account]
#
# Installs missing R packages, then submits the cross-validation array
# and a follow-up job that writes the summaries once the array is done.

set -euo pipefail

cd "$(dirname "$0")/.."

account="${1:-}"

if [ -z "$account" ]; then

  accounts=$(sacctmgr -nP show associations user="$USER" format=account | sort -u)

  if [ "$(echo "$accounts" | grep -c .)" -ne 1 ]; then
    echo "Pick a Slurm account and rerun: bash greatlakes/submit.sh <account>"
    echo "Your accounts:"
    echo "$accounts"
    exit 1
  fi

  account="$accounts"
fi

echo "Slurm account: $account"

module load R

mkdir -p greatlakes/logs

Rscript greatlakes/setup_packages.R

# Write run_settings.rds once, so array tasks never race to create it.
Rscript greatlakes/cv_task.R --init

n_tasks=$(Rscript greatlakes/cv_task.R --count | tr -d ' ')

echo "Submitting $n_tasks array tasks."

array_id=$(sbatch --parsable --account="$account" \
  --array="1-${n_tasks}" greatlakes/cv_array.sbatch)

aggregate_id=$(sbatch --parsable --account="$account" \
  --dependency="afterok:${array_id}" --kill-on-invalid-dep=yes \
  greatlakes/aggregate.sbatch)

echo "Array job:     $array_id"
echo "Aggregate job: $aggregate_id (starts when the array succeeds)"
echo "Check with:    squeue -u $USER"
