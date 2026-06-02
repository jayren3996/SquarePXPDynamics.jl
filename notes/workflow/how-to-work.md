# How to work on this repo (agent + author conventions)

## Engineering principles (always)

1. **Think before coding** — state assumptions, surface tradeoffs, ask when
   genuinely unclear; don't hide confusion.
2. **Simplicity first** — minimum code that solves the actual problem; nothing
   speculative.
3. **Surgical changes** — touch only what's needed, match surrounding style; only
   remove orphans your change created.
4. **Goal-driven** — define verifiable success criteria up front, and be willing
   to report that a change does NOT help rather than ship dead complexity.

## Autonomy

When the author explicitly opts into "push straight through", proceed without
per-commit pauses. Otherwise commit only when asked, and surface genuine forks
(production defaults, valued tests) for a decision rather than deciding unilaterally.

## Compute

64-core large-memory shared cluster, full usage rights. Under high load (`uptime`
≳ 100), run the Julia test suite with `JULIA_TEST_NWORKERS=8` (not 16) to avoid
worker starvation, and always capture Julia's real exit code (`PIPESTATUS`), since
the pmap harness can under-report failures. The exact 16-site contraction is
memory-heavy at D≥5 — `GC.gc()` before each call and expect OOM at D≥6.

## Validation discipline

Never make a physics claim from a single time point or from simple/local
observables. See `../methodology/revival-validation.md`.
