# Agent Handoff

## Current Objective

No active implementation branch. The repo is ready for the next
research/engineering milestone selection after the GPT PXP roadmap completion
was merged, pushed, verified, and branch-cleaned.

## What Was Just Done

Merged the GPT PXP roadmap completion branch into `main`, pushed `main` to
GitHub, verified the merged tree, removed the feature worktree, and deleted
merged cleanup branches. Then refreshed short-term and mid-term memory plus a
new completion note so future sessions start from the current pushed state.

## Important Files Touched

- `memory/mid_term/project_goals.md`
- `memory/mid_term/architecture.md`
- `memory/mid_term/milestones.md`
- `memory/mid_term/open_questions.md`
- `memory/mid_term/experiments_and_results.md`
- `memory/short_term/*.md`
- `docs/superpowers/notes/2026-05-17-gpt-roadmap-completion.md`

## Commands Run

- `git status --short --branch`
- `git worktree list`
- `git merge --no-ff codex/complete-gpt-pxp-roadmap`
- `julia --project=. -e 'using Pkg; Pkg.test()'`
- `git diff --check HEAD~1 HEAD`
- `git push origin main`
- `git worktree remove ...`
- `git branch -d ...`

## Tests/Results

- Feature-branch full-suite verification:
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65356/65356` in
  `7m28.3s`.
- Post-merge full-suite verification on `main`:
  `julia --project=. -e 'using Pkg; Pkg.test()'` passed `65356/65356` in
  `4m57.8s`.
- `git diff --check HEAD~1 HEAD` passed.
- `origin/main` points to `98e1ad7`.

## Known Problems

- CTM-aware/full-update evolution is not implemented.
- Full tensor snapshot persistence for ScarFinder candidates is not implemented.
- Publication-grade physics audit sweeps across `dt`, `D`, `chi`, cutoff, unit
  cell, and update scheme remain future work.
- Return/fidelity proxy, richer CTM observables, transfer-matrix correlation
  length, and energy-variance-quality observables remain future work.

## Next Recommended Actions

Choose the next milestone: CTM-trusted ScarFinder audit campaign, tensor
snapshot persistence, expanded CTM observables, or CTM-aware/full-update design.

## Things Not To Do

- Do not treat simple/local observables as CTMRG-quality physics measurements.
- Do not claim publication-grade ScarFinder validation from trusted CTM plumbing
  alone; convergence/audit campaigns are still required.
- Do not rerun `project-memory-curator` unless explicitly requested.
- Do not delete older notes just because they are summarized here.
