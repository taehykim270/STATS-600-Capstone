#!/bin/bash
# Run ON YOUR MAC, from anywhere:
#   bash greatlakes/launch.sh <uniqname> [slurm_account]
#
# Copies the project (including finished checkpoints, which the cluster
# then skips) to ~/capstone on Great Lakes and submits the jobs.
# Needs the U-M VPN when off campus. Reuses an open login if there is one
# (see ControlPath below); otherwise you are asked for password + Duo.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: bash greatlakes/launch.sh <uniqname> [slurm_account]"
  exit 1
fi

user="$1"
account="${2:-}"
host="greatlakes.arc-ts.umich.edu"
remote_dir="capstone"
ssh_options="-o ControlPath=$HOME/.ssh/cm-greatlakes"

cd "$(dirname "$0")/.."

rsync -av -e "ssh $ssh_options" \
  --include='model_selection.R' \
  --include='weather.rda' \
  --include='greatlakes/' \
  --include='greatlakes/*.R' \
  --include='greatlakes/*.sh' \
  --include='greatlakes/*.sbatch' \
  --include='weather_selection_*/' \
  --include='weather_selection_*/*.rds' \
  --exclude='*' \
  ./ "$user@$host:$remote_dir/"

ssh $ssh_options "$user@$host" \
  "cd $remote_dir && bash greatlakes/submit.sh $account"
