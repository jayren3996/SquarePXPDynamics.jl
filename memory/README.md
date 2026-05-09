# Project Memory

This directory is curated project memory for future Codex and AI-agent sessions working on the triangular-lattice PXP ScarFinder PEPS code. It is not a raw notes archive and should stay smaller than the source notes.

Before substantial work, normally read these files in order:

1. `memory/README.md`
2. `memory/mid_term/project_goals.md`
3. `memory/mid_term/architecture.md`
4. `memory/mid_term/decision_log.md`
5. `memory/short_term/current_state.md`
6. `memory/short_term/handoff.md`

Then read relevant long-term files for the task, and only then inspect the original notes, source, and tests.

## Memory Boundaries

- `long_term/`: stable scientific context, physics definitions, basis conventions, theoretical constraints, and literature framing expected to remain valid across many sessions.
- `mid_term/`: project goals, architecture, decisions, milestones, known results, open questions, and implementation direction.
- `short_term/`: immediate state, active tasks, recent commands/results, blockers, and a practical handoff. These files are rewriteable and should stay brief.

Do not put whole papers, raw extracted text, large implementation plans, speculative design branches, or transient debugging output here. Keep detailed derivations and long notes in `Notes/` or `docs/`, and reference them from memory.

## Updating Rules

- Run the project-memory-curator skill only when explicitly requested.
- Preserve source references for nontrivial claims, using paths such as `Notes/literature_review.md`, `src/Evolution.jl`, or `test/test_evolution.jl`.
- Mark uncertainty with labels such as `Confirmed:`, `Inferred:`, `Speculative:`, `Open question:`, `Deprecated / stale:`, or `Superseded:`.
- Do not silently resolve contradictions. Record them in `memory/mid_term/open_questions.md`.
- Put important scientific, architectural, implementation, or workflow choices in `memory/mid_term/decision_log.md`.
- After substantial work, update short-term handoff files only when asked or when explicitly curating memory again.

## Decision Log

Use `memory/mid_term/decision_log.md` for decisions future agents might otherwise reverse, including package structure, basis conventions, update-backend boundaries, and ScarFinder workflow choices. Mark superseded decisions rather than deleting them.
