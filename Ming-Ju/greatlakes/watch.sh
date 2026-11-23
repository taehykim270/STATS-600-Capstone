#!/bin/bash
# Run ON YOUR MAC, after launch.sh:
#   bash greatlakes/watch.sh <uniqname> [seconds_between_checks]
#
# Shows a live progress bar, queue counts and estimated time remaining.
# Stops by itself when the final summary exists. Ctrl-C to quit early;
# that only stops the display, not the cluster jobs.

set -uo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: bash greatlakes/watch.sh <uniqname> [seconds_between_checks]"
  exit 1
fi

user="$1"
interval="${2:-30}"
host="greatlakes.arc-ts.umich.edu"
remote_dir="capstone"
ssh_options="-o ControlPath=$HOME/.ssh/cm-greatlakes"

cd "$(dirname "$0")/.."

total=$(Rscript greatlakes/cv_task.R --count 2>/dev/null | tr -d ' ')

if ! [[ "$total" =~ ^[0-9]+$ ]]; then
  echo "Could not work out the number of tasks."
  exit 1
fi

# Only count the run that model_selection.R currently points at.
output_dir=$(sed -n 's/^output_directory <- "\(.*\)"$/\1/p' model_selection.R)

if [ -z "$output_dir" ]; then
  echo "Could not read output_directory from model_selection.R."
  exit 1
fi

echo "Run: $output_dir ($total checkpoints)"

# One line back from the cluster: done running pending failed finished
remote_command="
  cd $remote_dir 2>/dev/null || { echo 'ERR no project folder'; exit; }
  done_count=\$(ls $output_dir/*_year_*_h_*.rds 2>/dev/null | wc -l)
  states=\$(squeue -u $user -h -r -o '%T' 2>/dev/null)
  running=\$(echo \"\$states\" | grep -c RUNNING)
  pending=\$(echo \"\$states\" | grep -c PENDING)
  failed=\$(grep -l -i 'error\|killed\|cancelled' greatlakes/logs/cv_*.log 2>/dev/null | wc -l)
  finished=0
  ls $output_dir/best_settings.csv >/dev/null 2>&1 && finished=1
  echo \$done_count \$running \$pending \$failed \$finished
"

format_duration() {
  local seconds=$1
  if [ "$seconds" -ge 3600 ]; then
    printf '%dh %02dm' $((seconds / 3600)) $((seconds % 3600 / 60))
  elif [ "$seconds" -ge 60 ]; then
    printf '%dm %02ds' $((seconds / 60)) $((seconds % 60))
  else
    printf '%ds' "$seconds"
  fi
}

start_time=$(date +%s)
start_done=""
bar_width=40

echo "Watching Great Lakes every ${interval}s. Ctrl-C to stop watching."

while true; do

  reply=$(ssh $ssh_options -o BatchMode=yes -o ConnectTimeout=15 \
    "$user@$host" "$remote_command" 2>/dev/null | tail -1)

  read -r done_count running pending failed finished <<< "$reply"

  if ! [[ "${done_count:-}" =~ ^[0-9]+$ ]]; then
    printf '\r\033[K%s  no answer from the cluster (%s); retrying' \
      "$(date +%H:%M:%S)" "${reply:-login expired?}"
    sleep "$interval"
    continue
  fi

  now=$(date +%s)

  if [ -z "$start_done" ]; then
    start_done=$done_count
  fi

  elapsed=$((now - start_time))
  gained=$((done_count - start_done))
  remaining=$((total - done_count))

  if [ "$remaining" -le 0 ]; then
    eta="fits done, writing summaries"
  elif [ "$gained" -gt 0 ]; then
    eta="about $(format_duration $((remaining * elapsed / gained))) left"
  elif [ "$running" -eq 0 ] && [ "$pending" -gt 0 ]; then
    eta="waiting in the queue"
  else
    eta="estimating"
  fi

  filled=$((done_count * bar_width / total))
  [ "$filled" -gt "$bar_width" ] && filled=$bar_width
  bar=$(printf '%*s' "$filled" '' | tr ' ' '#')
  bar="$bar$(printf '%*s' $((bar_width - filled)) '' | tr ' ' '.')"

  printf '\r\033[K%s [%s] %d/%d (%d%%) | run %d, queue %d, failed %d | %s' \
    "$(date +%H:%M:%S)" "$bar" "$done_count" "$total" \
    $((done_count * 100 / total)) "$running" "$pending" "$failed" "$eta"

  if [ "$finished" -eq 1 ]; then
    printf '\nFinished. Get the results with: bash greatlakes/fetch.sh %s\n' "$user"
    exit 0
  fi

  if [ "$running" -eq 0 ] && [ "$pending" -eq 0 ] && [ "$remaining" -gt 0 ]; then
    printf '\nNo jobs are queued but %d checkpoints are missing.\n' "$remaining"
    printf 'Look at greatlakes/logs on the cluster, or rerun launch.sh to resubmit.\n'
    exit 1
  fi

  sleep "$interval"
done
