# Agent Handoff

## Current objective

Create repository-local curated project memory for the triangular PXP ScarFinder PEPS project.

## What was just done

Added `memory/` with long-term, mid-term, and short-term files summarizing physics context, conventions, literature direction, architecture, decisions, milestones, open questions, current state, active tasks, next steps, and this handoff.

## Important files touched

- `memory/README.md`
- `memory/long_term/*.md`
- `memory/mid_term/*.md`
- `memory/short_term/*.md`
- `AGENTS.md`

## Commands run

- `rg --files ...` to inspect repo docs/source/tests.
- `git status --short && git log --oneline -8`
- Several `sed -n` reads of README, AGENTS, Notes, docs, source, and tests.

## Tests/results

No full Julia test suite was run during memory creation. A test/source contradiction was found around `D>1` non-product Simple Update behavior.

## Known problems

`test/test_simple_update.jl` expects a direct `D>1` non-product star update to throw, while source/docs/evolution tests describe or use a `D>1` local product-projection path.

## Next recommended actions

Resolve the `D>1` Simple Update test/source drift, then run `julia --project=. -e 'using Pkg; Pkg.test()'`.

## Things not to do

Do not create a nested Julia package or separate PEPS environment. Do not treat this PEPS layer as a general tensor-network library. Do not run project-memory-curator again unless explicitly requested.
