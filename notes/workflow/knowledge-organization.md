# Knowledge organization: notes/, not memory

Decision (2026-06-02, author): this repo does NOT use a separate "memory" layer.
A single systematic `notes/` tree (this directory, with topic subfolders) is the
one source of truth for durable knowledge. The previous split — a Claude Code
auto-memory at `~/.claude/...`, the repo's `memory/{long,mid,short}_term/`, and a
`Notes/` dir — was retired and consolidated here; it caused duplicated, drifting
records (the same finding rewritten in several places as conclusions changed).

**How to apply (for the agent):**
- Do NOT maintain the Claude auto-memory for this project. Put durable findings,
  decisions, conventions, and lessons in the relevant `notes/<subfolder>/` file.
- Each topical `*.md` is a concise current summary; keep the dated
  `YYYY-MM-DD-*.md` experiment records alongside as the detailed provenance.
- When a conclusion is overturned, EDIT the topical note in place (don't leave the
  stale claim in a second location). Move superseded notes to `notes/archive/`.
- Design specs/plans stay under `docs/superpowers/{specs,plans}/`.
