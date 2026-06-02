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

## Active priority: execute the 3-stage roadmap

The consolidated review and backlog live in
`docs/superpowers/notes/2026-06-02-package-review-and-roadmap.md`. **Read it
first** — it is the source of slices for this campaign and lists exact
file:line targets. The three user stages are: (1) reliable iPEPS dynamics that
match ED, (2) bond truncation + proper regauging, (3) ScarFinder that finds a
better-than-Néel initial state.

The owner made three decisions on 2026-06-02:

1. **Stage-3 "better than Néel" metric = staggered-magnetization revival.**
   Score the longer-lived / larger *return* of the `(n_even − n_odd)` order
   parameter over the trajectory time series — NOT instantaneous imbalance
   (the current RevivalObjective is wrong: it rewards t=0 Néel itself).
2. **Slimming = aggressive delete.** These are SANCTIONED for removal once you
   `grep -rn` confirm no remaining caller in `src/`, `test/`, `scripts/`
   (delete tests that only cover the removed code in the same slice):
   - the TFIM benchmark axis (`Benchmarks.jl`, `FiniteTFIMReference.jl`,
     `FiniteMPSTFIMReference.jl`, `TFIMStarModel` + `tfim_*` observables);
   - the orphaned CTM gauge module (`CTMGaugeReadiness.jl` / `fix_bond_gauge!`
     and its 9 exports) — never called by evolution/compression;
   - placeholder ScarFinder objectives (`TargetEnergyObjective`,
     `LowVarianceObjective`) + `energy_variance_proxy` plumbing; `JSONCandidateStore`;
   - the duplicate PXP campaign driver (`run_pxp_audit_campaign` + its struct
     family) — keep `run_pxp_larger_d_benchmark`;
   - `SquarePEPS.jl`, `GaugeDiagnostics.jl` (unused), dead exports (census in
     the note), and stale/kagome docs (archive to `docs/superpowers/archive/`).
3. **Stage-3 geometry = even-Lx/Ly cells with `schedule = :serial`** (Néel
   works there today; `:five_color` forces 5-divisible odd-tiling cells that
   cannot host a perfect checkerboard).

### Slice priority order — work the highest item with a well-scoped slice

1. **Harden the test gate FIRST** (so later physics changes are catchable):
   convert every `density_error_simple > 1e-4` lower-bound assertion to
   `@test_broken` and add `density_error_exact_finite < 1e-6` upper-bound gates;
   add a `@test_throws` pinning the `:z_up` blockade-forbidden error; add a
   nightly (`SQUAREPXP_EXTENDED_TESTS`) real-`measure_ctm`-vs-ED D=2 test.
2. **Aggressive slimming** — one sanctioned target cluster per pass.
3. **Stage-1 simple-observable D>1 bug**: the simple-gauge density contraction
   in `Observables.jl` reports ~1.9e-4 error at D=2 while `exact_density_finite`
   is ~1e-7. Fix it, validated against `exact_density_finite`. Small steps.
4. **Stage-3 staggered-mag revival**: a per-iteration `(n_even−n_odd)` time
   series + redefine `RevivalObjective` to score the post-collapse return; then
   an `scarfinder_search` outer layer ranking candidate states vs a Néel
   baseline. Large feature — build in small, separately-committed slices.
5. **Stage-2 Vidal √λ regauge** (principled successor to `rel_floor`). Largest;
   approach last, with ED/exact-finite verification at tight cutoff for D≤4.

**Physics-changing slices (items 3-5) MUST verify against ED or
`exact_density_finite`** — never claim a physics improvement from
`measure_simple` alone (it is wrong for D>1). If a physics slice cannot be
verified within the pass budget, leave a TODO in `handoff.md` and pick a
slimming or test-hardening slice instead.

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

2. Run the fast physics-correctness gate plus the test file(s) covering your
   change. The gate runs the Stage-1/2 regression tests in one warm process:
   `timeout 540 julia --project=test scripts/fast_gate.jl`. For a change under
   `src/` that touches evolution, truncation, observables, or scarfinder, ALSO
   run the relevant focused test file, e.g.
   `timeout 540 julia --project=test test/runtests.jl test_ipeps_evolution.jl`.
   Always wrap with `timeout`. Reserve the full `Pkg.test()` (~10 min) for a
   slice you cannot otherwise convince yourself is safe; it will not fit the
   pass budget alongside much else.

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
