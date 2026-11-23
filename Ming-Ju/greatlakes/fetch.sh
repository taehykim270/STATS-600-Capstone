#!/bin/bash
# Run ON YOUR MAC:
#   bash greatlakes/fetch.sh <uniqname>
#
# Shows the queue, then copies results and logs back into this project.
# Safe to run repeatedly while jobs are still going.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: bash greatlakes/fetch.sh <uniqname>"
  exit 1
fi

user="$1"
host="greatlakes.arc-ts.umich.edu"
remote_dir="capstone"
ssh_options="-o ControlPath=$HOME/.ssh/cm-greatlakes"

cd "$(dirname "$0")/.."

ssh $ssh_options "$user@$host" \
  "squeue -u $user -h -o '%T' | sort | uniq -c; \
   echo \"Checkpoints: \$(ls $remote_dir/weather_selection_*/*_year_*_h_*.rds 2>/dev/null | wc -l)\""

rsync -av -e "ssh $ssh_options" \
  --include='weather_selection_*/' \
  --include='weather_selection_*/*' \
  --include='greatlakes/' \
  --include='greatlakes/logs/' \
  --include='greatlakes/logs/*' \
  --exclude='*' \
  "$user@$host:$remote_dir/" ./
