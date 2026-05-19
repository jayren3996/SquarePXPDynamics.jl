# Autonomous Improvement Pass — SquarePXPDynamics.jl

You are working on `/data/djxg096/SquarePXPDynamics.jl`, a Julia package for
square-lattice PXP dynamics. The long-term objective is a working **ScarFinder
for the PXP model on a 2D infinite square lattice** built on iPEPS + CTM
(PEPSKit). This invocation is one pass in an autonomous loop — each pass runs
in fresh context, so the `memory/` directory and git history are your only
continuity.

## Read first (do not skip)

1. `AGENTS.md`
2. `memory/README.md` and the files it lists, in order:
   - `memory/mid_term/project_goals.md`
   - `memory/mid_term/architecture.md`
   - `memory/mid_term/decision_log.md`
   - `memory/short_term/current_state.md`
   - `memory/short_term/handoff.md`
3. Open the **Known Problems** and **Next Recommended Actions** sections of
   `memory/short_term/handoff.md` — your slice for this pass should normally
   come from there.

Then run `git status -s` and `git pull --ff-only origin main` so you start
from a clean, up-to-date `main`. If any **tracked** files are modified or
staged (lines starting with ` M`, `M `, `A `, `D `, `R `, `??` with a path
under `src/`, `test/`, `memory/`, `prompts/`, or `scripts/`), or if the pull
is not fast-forward, write `NOOP: dirty tree / non-FF pull, skipping` to
stdout and exit without changes. Stray untracked artifacts/logs are fine.

## Pick one slice

Choose **one** concrete, well-scoped slice — small enough to finish, test, and
commit in this pass. Good candidate shapes:

- A small refactor with no behavior change (clarify naming, extract a helper).
- A new unit test that exercises a currently uncovered branch.
- A docstring tightening or a corrected error message.
- A documented invariant added as a comment where the WHY is non-obvious.
- A short performance probe whose run completes in under ~5 minutes.
- A memory file update reconciled with what is actually in the code (drift
  cleanup) — but only if you also fix the underlying drift.

Prefer slices that move us toward the ScarFinder-on-iPEPS goal: CTM
observable correctness/throughput, gauge fixing, PXP energy/density
operators, finite-vs-infinite reconciliation, and the test harness around
those.

If no well-scoped slice fits, write a single line starting with `NOOP:` and a
short reason, then exit. Doing nothing is the correct outcome for some
passes.

## Hard limits for one pass

- Up to ~20 minutes of wall clock total.
- At most **one** git commit, touching at most ~5 files.
- No multi-hour benchmarks — if a slice needs one, instead add a TODO entry to
  `memory/short_term/handoff.md` describing the run, and exit.
- No tool call you expect to exceed 10 minutes; split it or defer it.
- Do not run shell commands with `run_in_background: true` unless you also
  set an explicit, finite `timeout` on the command itself. The launcher
  wraps each pass in a 25-minute hard timeout, but a backgrounded poll for
  a never-arriving exit-file is the most common way passes get killed
  unproductively. If you must run something long, prefer foreground with
  `timeout 540 <cmd>`.

## Things you must not do

- Do not delete tracked files, force-push, rebase published commits, or amend
  pushed commits.
- Do not skip git hooks (`--no-verify`) or bypass signing.
- Do not disable, skip, or `@test_skip` existing tests to get a change to
  pass. If a test legitimately needs to change, change it deliberately and
  explain in the commit message.
- Do not restart 7x7 ED dynamics or treat simple/local observables as
  CTMRG-quality physics measurements (see `handoff.md` → Things Not To Do).
- Do not run `project-memory-curator`.
- Do not start a second long-running background process; this loop is the
  only one.
- Do not edit, rename, or delete the files that drive this loop:
  `prompts/autonomous-improve-prompt.md` and
  `scripts/run_autonomous_loop.sh`. If they look wrong, write a NOOP line
  describing the concern and exit.

## After making changes

1. Run the targeted test that covers what you changed, e.g.
   `julia --project=test test/runtests.jl test_pepskit_measurements.jl`.
2. If tests fail, attempt one fix. If still failing, `git restore` your
   changes and write a single-line `NOOP: <reason>` to stdout and exit.
3. If tests pass: `git add` only the files you intentionally modified (no
   wildcards), commit with a concise message in the existing style
   (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`), then
   `git push origin main`.
4. If the slice represented a non-trivial decision (algorithmic choice,
   convention, scope cut), append a short dated entry to
   `memory/mid_term/decision_log.md` with `Confirmed:` / `Inferred:` /
   `Open question:` labels per `memory/README.md`.
5. If you finished or materially changed a "Next Recommended Action" or
   "Known Problem", update `memory/short_term/handoff.md` and
   `memory/short_term/current_state.md` accordingly. Stage and commit these
   together with the code change when they belong to the same slice.

## Output discipline

Stdout from this pass becomes the loop log the user reads when they return.
Be terse. Useful lines to emit at minimum:

- One line stating the slice you chose, with the file path it touches.
- One line for the test command and its pass/fail result.
- One line for the commit hash and message, or the `NOOP:` line if you
  exited early.

Do not narrate internal deliberation. Do not summarize the repo or restate
the prompt back. Skip end-of-turn pleasantries.
