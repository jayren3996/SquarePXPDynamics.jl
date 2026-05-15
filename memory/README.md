# Project Memory

This directory is curated project memory for future Codex and AI-agent sessions
working on `SquarePXPDynamics`. It is a recovery layer for scientific context,
architecture, decisions, current state, and handoff notes. It is not a dump of
all repository notes.

Before substantial work, agents should normally read these files in order:

1. `memory/README.md`
2. `memory/mid_term/project_goals.md`
3. `memory/mid_term/architecture.md`
4. `memory/mid_term/decision_log.md`
5. `memory/short_term/current_state.md`
6. `memory/short_term/handoff.md`

Then read relevant long-term files for the task:

- `memory/long_term/physics_context.md`
- `memory/long_term/literature_context.md`
- `memory/long_term/definitions_and_conventions.md`
- `memory/long_term/theoretical_constraints.md`

## Memory Tiers

- Long-term memory contains stable physics context, definitions, conventions,
  literature anchors, and theoretical constraints. It should not contain
  transient TODOs.
- Mid-term memory contains project goals, architecture, milestones, decisions,
  experiment summaries, and unresolved questions.
- Short-term memory contains the current working state, active tasks, immediate
  next steps, blockers, commands, and handoff notes. Keep it concise and easy
  to rewrite.

## What Not To Store

Do not copy whole design documents, large logs, generated artifacts, raw test
output, or speculative ideas without labeling them. Preserve detailed source
material in `notes/`, `docs/`, source files, tests, and git history.

## Update Rules

- Add source references for nontrivial entries, such as `Source:
  notes/foo.md`, `Source: src/bar.jl`, or `Source: git log`.
- Use confidence labels where useful: `Confirmed:`, `Inferred:`,
  `Speculative:`, `Open question:`, `Deprecated / stale:`, `Superseded:`.
- Record important scientific, architectural, implementation, and workflow
  choices in `memory/mid_term/decision_log.md`.
- Mark superseded decisions clearly instead of deleting them.
- If sources contradict each other, record the contradiction in
  `memory/mid_term/open_questions.md`; do not silently overwrite it.
- Update short-term handoff files when explicitly asked after substantial work.

## Skill Invocation Rule

Run `project-memory-curator` only when explicitly requested by the user, for
example to curate project memory, update repo memory, or refresh memory files.
Do not run it for ordinary coding, debugging, documentation, planning, or code
review unless the user asks for it.
