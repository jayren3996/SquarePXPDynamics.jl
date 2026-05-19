#!/usr/bin/env bash
# Autonomous Claude improvement loop, detached in tmux so it survives logout.
#
# Usage:
#   scripts/run_autonomous_loop.sh start    # start the loop
#   scripts/run_autonomous_loop.sh attach   # view the live loop (Ctrl-b d to detach)
#   scripts/run_autonomous_loop.sh status   # is it running?
#   scripts/run_autonomous_loop.sh stop     # kill the session
#   scripts/run_autonomous_loop.sh tail     # tail the log without attaching
#
# Env knobs:
#   SLEEP_SECONDS   seconds to sleep between passes (default 1800 = 30 min)
#   MAX_PASSES      stop after N passes (default 0 = forever)

set -euo pipefail

SESSION="claude-improve"
REPO="/data/djxg096/SquarePXPDynamics.jl"
PROMPT_FILE="$REPO/prompts/autonomous-improve-prompt.md"
LOG_FILE="$REPO/logs/autonomous-loop.log"
SLEEP_SECONDS="${SLEEP_SECONDS:-1800}"
MAX_PASSES="${MAX_PASSES:-0}"

cmd="${1:-status}"

case "$cmd" in
  start)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "Session '$SESSION' is already running. Use 'attach' or 'stop'."
      exit 1
    fi
    if [[ ! -f "$PROMPT_FILE" ]]; then
      echo "Prompt file missing: $PROMPT_FILE"
      exit 1
    fi
    if ! command -v claude >/dev/null 2>&1; then
      echo "claude CLI not found on PATH"
      exit 1
    fi
    if ! command -v tmux >/dev/null 2>&1; then
      echo "tmux not found on PATH"
      exit 1
    fi
    mkdir -p "$(dirname "$LOG_FILE")"

    # The loop body. Each pass is a fresh `claude -p` invocation so context
    # cannot grow unbounded. Memory files + git provide continuity.
    #
    # Auth survival across user logout: we strip the VS Code git askpass
    # env vars (their IPC socket dies when the user disconnects) and force
    # git to fail fast instead of hanging on a TTY prompt. Push credentials
    # come from ~/.git-credentials populated via `credential.helper=store`.
    read -r -d '' LOOP <<EOF || true
cd "$REPO"
unset GIT_ASKPASS SSH_ASKPASS
unset VSCODE_GIT_ASKPASS_NODE VSCODE_GIT_ASKPASS_MAIN VSCODE_GIT_ASKPASS_EXTRA_ARGS VSCODE_GIT_IPC_HANDLE
export GIT_TERMINAL_PROMPT=0
pass=0
while true; do
  pass=\$((pass + 1))
  echo "=== pass \$pass @ \$(date -Iseconds) ===" | tee -a "$LOG_FILE"
  claude -p "\$(cat "$PROMPT_FILE")" --dangerously-skip-permissions 2>&1 | tee -a "$LOG_FILE"
  echo "--- end pass \$pass; sleeping ${SLEEP_SECONDS}s ---" | tee -a "$LOG_FILE"
  if [[ "$MAX_PASSES" -gt 0 && "\$pass" -ge "$MAX_PASSES" ]]; then
    echo "=== reached MAX_PASSES=$MAX_PASSES, exiting ===" | tee -a "$LOG_FILE"
    break
  fi
  sleep "$SLEEP_SECONDS"
done
EOF

    tmux new-session -d -s "$SESSION" "bash -c '$LOOP'"
    echo "Started tmux session '$SESSION'."
    echo "  Log:    $LOG_FILE"
    echo "  Attach: tmux attach -t $SESSION   (Ctrl-b d to detach)"
    echo "  Tail:   $0 tail"
    echo "  Stop:   $0 stop"
    ;;

  stop)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      tmux kill-session -t "$SESSION"
      echo "Killed '$SESSION'."
    else
      echo "No session '$SESSION'."
    fi
    ;;

  attach)
    tmux attach -t "$SESSION"
    ;;

  status)
    if tmux has-session -t "$SESSION" 2>/dev/null; then
      echo "Running. Attach: tmux attach -t $SESSION"
    else
      echo "Not running."
    fi
    ;;

  tail)
    tail -f "$LOG_FILE"
    ;;

  *)
    echo "Usage: $0 {start|attach|status|stop|tail}"
    exit 1
    ;;
esac
