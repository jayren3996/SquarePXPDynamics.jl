# notes/

Single source of truth for durable project knowledge. Replaces the previous
split across a Claude auto-memory, `memory/{long,mid,short}_term/`, and `Notes/`
(retired 2026-06-02). Design specs and plans remain under `docs/superpowers/`.

## Layout

- `methodology/` — how to evaluate the dynamics (validation rules, lessons). Read
  these before trusting any benchmark number.
- `conventions/` — the model, states, physics context, theoretical constraints.
- `stage1-dynamics/` — Stage 1: reliable iPEPS dynamics matching ED.
- `stage2-truncation/` — Stage 2: bond truncation, conditioning, environment.
- `stage3-scarfinder/` — Stage 3: better-than-Néel initial-state search.
- `workflow/` — how to work on this repo (conventions for the agent + author).
- `archive/` — superseded/historical notes kept for provenance.
- `project_goals.md`, `open_questions.md` — the three-stage goal and live questions.

Each topical `*.md` is a concise current summary; the dated `YYYY-MM-DD-*.md`
files alongside are the detailed experiment records behind them.

## Start here

1. `methodology/revival-validation.md` — **the load-bearing lesson**: judge
   revival dynamics by the n(t) TRAJECTORY, never a single time.
2. `project_goals.md` — the three-stage objective.
3. `stage2-truncation/improvement-roadmap.md` — the current top of the backlog.
