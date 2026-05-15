# Agent Handoff

## Current objective

Create a curated repository-local project memory system for
`SquarePXPDynamics`.

## What was just done

Inspected README, project files, notes, docs, source/test layout, git status,
recent commits, and the TFIM benchmark handoff note. Created `memory/` with
long-term, mid-term, and short-term memory files plus a top-level `AGENTS.md`
pointer.

## Important files touched

- `AGENTS.md`
- `memory/README.md`
- `memory/long_term/*.md`
- `memory/mid_term/*.md`
- `memory/short_term/*.md`

## Commands run

- `git status --short --branch`
- `git log --oneline -10`
- `rg --files`
- `rg -n '^#|^##|^###' notes docs/superpowers`
- `sed -n ...` on README, project files, notes, docs, source entrypoint, and
  test runner

## Tests/results

No Julia tests were run for this documentation-only memory update. Verification
was limited to reading the generated Markdown, checking repository status, and
running `git diff --check`.

## Known problems

- The memory system is newly created and uncommitted.
- The completed TFIM benchmark work is in a separate worktree/branch and is not
  merged into the current `main` checkout.
- Existing notes contain some superseded descriptions of earlier repo state;
  current README/source should be preferred for architecture status.

## Next recommended actions

Review the memory files, then stage and commit them if they should be kept. Next
decide whether to merge, push, or further extend `codex/infinite-tfim-benchmark`.

## Things not to do

- Do not treat simple/local observables as CTMRG-quality physics measurements.
- Do not delete older notes just because they are summarized here.
- Do not rerun `project-memory-curator` unless explicitly requested.
- Do not silently overwrite contradictions; record them in
  `memory/mid_term/open_questions.md`.
