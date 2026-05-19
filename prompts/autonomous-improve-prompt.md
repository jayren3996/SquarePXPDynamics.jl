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

### Stranded commits from a prior timed-out pass

Right after the pull, run `git log --oneline origin/main..main`. If it
returns anything, a previous pass committed but didn't push (it was killed
by the 25-min wall before its push). Before picking a new slice:

- Inspect the stranded commit's diff with `git show HEAD`. If it looks
  reasonable (small diff, plausibly correct, tests not obviously broken),
  push it (`git push origin main`) as the entire slice for this pass and
  exit cleanly with a one-line summary. Do not rebase, reorder, or amend.
- If the diff looks wrong (truncated, mid-edit, accidentally deletes
  important code), `git reset --hard origin/main` to drop it, write
  `NOOP: dropped stranded commit <hash> — <reason>` and exit.

Ship-stranded or new-slice, **one per pass** — do not stack them.

## Active priority: code-quality simplification

The user has flagged that the codebase is **unnecessarily long** and asked
for a sustained code-quality / simplification campaign. Until this section
is removed from the prompt, **prefer slices that shorten or simplify code
without changing behavior**. Net-line-removal slices are the most valuable
shape right now. Look for:

- **Dead code**: unreferenced functions, unused exports, branches that
  cannot trigger, error paths that no caller can ever hit.
- **Over-engineered abstractions**: single-use helpers that could be
  inlined; layered indirection that hides nothing meaningful; "config
  objects" that wrap a single field.
- **Duplicate logic**: similar code in 2+ places that should be a shared
  helper, or vice-versa — a helper used in one place that should be
  inlined.
- **Redundant validation**: the same input check repeated at multiple
  layers when one is enough.
- **Comments that restate WHAT** the code does (the identifiers already
  say what). Keep WHY comments — hidden constraints, workarounds,
  surprising invariants.
- **Verbose error handling** for cases that can't happen given internal
  callers and framework guarantees (validate only at system boundaries).
- **Over-broad try/catch** that swallows or rewraps errors with no added
  information.

Within this campaign, prefer simplifying code that is on the path to
**ScarFinder on the 2D infinite square PXP** — i.e. iPEPS/CTM observables,
gauge fixing, PXP energy/density operators, finite-vs-infinite
reconciliation. Cleaning those modules pays off twice.

## Pick one slice

Choose **one** concrete, well-scoped slice — small enough to finish, test,
and commit in this pass. In addition to the simplification shapes above,
these are still valid:

- A small refactor with no behavior change (clarify naming, extract a
  helper).
- A new unit test that exercises a currently uncovered branch.
- A docstring tightening or a corrected error message.
- A documented invariant added as a comment where the WHY is non-obvious.
- A short performance probe whose run completes in under ~5 minutes.
- A memory file update reconciled with what is actually in the code (drift
  cleanup) — but only if you also fix the underlying drift.

If no well-scoped slice fits, write a single line starting with `NOOP:` and
a short reason, then exit. Doing nothing is the correct outcome for some
passes.

## Simplification guardrails

When the slice is a code-removal:

- Before deleting a function or method, grep `src/`, `test/`, `scripts/`,
  and `Notes/` for any remaining callers. If you cannot prove zero
  external usage, do not delete it — propose the deletion in a memory
  TODO instead.
- Before removing an export from `src/SquarePXPDynamics.jl`, grep for the
  symbol across the whole repo. Tests and downstream scripts count.
- Behavior-preserving means: same return values, same observable side
  effects, same error semantics (modulo error-message wording, which may
  improve). If you are unsure, the slice is not behavior-preserving and
  belongs in a different pass.
- Net deletions over ~50 lines should be split. Pick a tight cluster
  (one helper + its inlined call sites, or one module's dead branches)
  per pass.

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

## After making changes — commit first, then test

The 25-minute pass wall lands hardest when an agent is mid-Julia-test:
the commit never happens and the slice is lost. To make slices durable,
**commit before running tests**. If the test then kills the pass, the
commit survives locally and pass N+1 ships it via the stranded-commit
path above.

1. `git add` **only** the files you intentionally modified (no `git add
   .` or wildcards), and commit immediately with a concise message in the
   existing style (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`). Do
   not push yet.

2. Run the targeted test that covers what you changed, e.g.
   `julia --project=test test/runtests.jl test_pepskit_measurements.jl`.
   Always wrap it: `timeout 540 julia ...`. Prefer a focused test file
   over the full `runtests.jl`.

3. If tests **pass**: `git push origin main`.

4. If tests **fail**:
   - Attempt one fix. Because the original commit is still local-only,
     amending is safe: stage the fix and run `git commit --amend
     --no-edit`. Re-run the test.
   - If the test still fails: `git reset --hard HEAD~1` to drop the
     commit entirely. Write `NOOP: <slice> failed test <name>` and exit.

5. If the pass is killed mid-test (rc 124/137), the local commit survives
   and pass N+1 picks it up via the stranded-commit path.

If the slice represented a non-trivial decision (algorithmic choice,
convention, scope cut), append a short dated entry to
`memory/mid_term/decision_log.md` and stage+commit it as part of the same
slice — keep "one commit per pass".

If you finished or materially changed a "Next Recommended Action" or
"Known Problem", update `memory/short_term/handoff.md` and
`memory/short_term/current_state.md` accordingly, in the same commit.

## Output discipline

Stdout from this pass becomes the loop log the user reads when they return.
Be terse. Useful lines to emit at minimum:

- One line stating the slice you chose, with the file path it touches.
- One line for the test command and its pass/fail result.
- One line for the commit hash and message, or the `NOOP:` line if you
  exited early.

Do not narrate internal deliberation. Do not summarize the repo or restate
the prompt back. Skip end-of-turn pleasantries.
