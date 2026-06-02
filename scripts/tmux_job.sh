#!/usr/bin/env bash
# Launch a long job in its own detached tmux session, logging to logs/campaign/<name>.log.
# Survives client/network disconnects. tmux-1.8 compatible (no -c flag).
# Usage: scripts/tmux_job.sh <name> '<command>'
set -euo pipefail
REPO="/data/djxg096/SquarePXPDynamics.jl"
name="$1"; shift
cmd="$*"
log="$REPO/logs/campaign/${name}.log"
mkdir -p "$REPO/logs/campaign"
tmux kill-session -t "job-$name" 2>/dev/null || true
tmux new-session -d -s "job-$name" \
  "cd '$REPO' && { echo '=== START $name @ '\$(date -Iseconds)' ==='; $cmd; echo \"=== EXIT:\$? @ \$(date -Iseconds) ===\"; } 2>&1 | tee '$log'"
echo "launched job-$name -> tail -f $log"
