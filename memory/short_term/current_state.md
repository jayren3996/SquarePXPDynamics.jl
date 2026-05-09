# Current State

- Confirmed: Repository has uncommitted user work in notes, README, source, tests, deleted old literature files, new PEPS literature PDFs, and a new `.claude/` directory.
  - Source: `git status --short` on 2026-05-09

- Confirmed: No `memory/` directory existed before this curation; this session created the memory structure.
  - Source: `find memory -maxdepth 3 -type f`

- Confirmed: Current code exports geometry, gates, schedules, state constructors, observables, Simple Update, evolution helpers, and ScarFinder APIs from `src/TriangularPEPSDynamics.jl`.
  - Source: `src/TriangularPEPSDynamics.jl`

- Known issue: Test/source drift exists around direct `D>1` non-product star update behavior. See `memory/mid_term/open_questions.md`.
  - Source: `test/test_simple_update.jl`
  - Source: `src/SimpleUpdate.jl`
