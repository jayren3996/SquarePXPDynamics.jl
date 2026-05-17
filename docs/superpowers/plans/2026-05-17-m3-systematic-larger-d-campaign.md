# M3 Systematic Larger-D PXP Campaign Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:systematic-debugging` for any benchmark failure, schema issue, runtime failure, or unexpected trend. Use `superpowers:verification-before-completion` before reporting conclusions. If source changes are needed, use `superpowers:test-driven-development` before implementation.

**Goal:** Run a systematic all-zero/all-down PXP dynamics campaign comparing larger-D iPEPS trajectories with translationally invariant finite PBC ED references over sizes `3x3` through `7x7`.

**Architecture:** The initial state is translationally symmetric, so there is no central-region comparison. Compare translation-invariant one-qubit/local data through the ED global excitation density and the iPEPS density fields. Use exact finite iPEPS contraction for `3x3`; for larger PBC sizes, use symmetric PBC global observables and clearly label `density_simple` as diagnostic for `D > 1`. Energy spectra/levels are optional and require adding or using an explicit ED spectral path because the current M3 report does not expose levels.

**Tech Stack:** Julia 1.12, existing `scripts/pxp_larger_d_ed_benchmark.jl`, existing M3 JSON/CSV report schema, shell/Julia/Python summaries as convenient.

---

## Context And Constraints

- Use the newer Git if `/usr/bin/git` is too old:

```bash
export PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH
```

- Do not blindly set `JULIA_NUM_THREADS=42`. Linear algebra may already use BLAS/OpenBLAS threading.
- Before heavy runs, record the runtime threading setup:

```bash
julia --project=. -e 'using LinearAlgebra; println("BLAS threads = ", BLAS.get_num_threads()); println("Julia threads = ", Threads.nthreads())'
```

- Do not claim publication-grade physics.
- Do not implement CTM-aware/full-update evolution as part of this campaign.
- Do not use `density_simple` as exact finite truth for `D > 1`.
- Do not describe any PBC result as a central-region observable.
- All dynamics start from the all-zero/all-down state only.
- The main direct local observable is the translation-invariant one-site excitation density.
- For `3x3`, compare ED against exact finite iPEPS observables.
- For `4x4`, `5x5`, `6x6`, and `7x7`, compare valid symmetric PBC global one-site density data.
- Keep `log_norm_delta_abs` as a diagnostic only.

## Preflight

- [ ] **Step 1: Read project context**

Read:

```text
AGENTS.md
memory/README.md
memory/mid_term/project_goals.md
memory/mid_term/architecture.md
memory/mid_term/decision_log.md
memory/short_term/current_state.md
memory/short_term/handoff.md
docs/superpowers/notes/2026-05-17-m3-larger-d-pxp-ed-benchmark.md
```

- [ ] **Step 2: Check branch state**

Run:

```bash
git status --short --branch
git log --oneline --decorate -8
```

Expected: branch includes the M3 implementation commits ending near:

```text
e6eebff docs: document M3 PXP ED benchmark semantics
```

- [ ] **Step 3: Run focused preflight tests**

Run:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
julia --project=. test/runtests.jl test_pxp_ed_benchmark.jl test_pxp_validation.jl
```

Expected: both commands pass. If either fails, use `superpowers:systematic-debugging` before running campaigns.

## Campaign A: 3x3 Exact-Finite Longer Dynamics

- [ ] **Step 1: Run 3x3 exact-finite campaign**

Run:

```bash
mkdir -p artifacts/m3-systematic

SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02,0.01,0.005 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.20 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/3x3-exact-long.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/3x3-exact-long.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected: exact finite fields are populated.

- [ ] **Step 2: Optionally extend 3x3**

Run only if the `total_time = 0.20` campaign is practical:

```bash
SQUAREPXP_LARGERD_N=3 \
SQUAREPXP_LARGERD_DT=0.02,0.01,0.005 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.50 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/3x3-exact-longer.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/3x3-exact-longer.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

- [ ] **Step 3: Analyze 3x3 trend**

For each `dt` and `D`, summarize:

```text
abs(density_error_exact_finite)
abs(return_probability_error)
density_error_simple
max_truncerr
log_norm_delta_abs
reversibility_density_drift
```

Required conclusion:

- Does `D = 2,3,4` improve over `D = 1`?
- Is improvement monotonic, saturated, or degraded at longer time?
- Does smaller `dt` change the conclusion?
- Are exact finite density and return probability telling the same story?

## Campaign B: Size Sweep For Translationally Invariant One-Site Dynamics

- [ ] **Step 1: Run `3x3` through `6x6` size sweep**

Run:

```bash
SQUAREPXP_LARGERD_N=3,4,5,6 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2,3 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.10 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/size-sweep-3to6.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/size-sweep-3to6.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected:

- `3x3` rows use `observable_mode = exact_finite`.
- `4x4`, `5x5`, and `6x6` rows use `observable_mode = symmetric_pbc_ed_global`.
- The ED one-site observable is `ed_excitation_density`.

- [ ] **Step 2: Optionally extend size sweep to `D = 4`**

Run only if the first sweep is practical:

```bash
SQUAREPXP_LARGERD_N=3,4,5,6 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2,3,4 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.20 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/size-sweep-3to6-longer.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/size-sweep-3to6-longer.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

- [ ] **Step 3: Analyze finite-size trend**

For each `n`, `dt`, and `D`, summarize:

```text
ed_excitation_density
ipeps_simple_density
ipeps_exact_finite_density for 3x3 only
density_error_simple diagnostic only for D > 1
density_error_exact_finite for 3x3 only
return_probability_error when populated
ed_runtime_seconds
ipeps_runtime_seconds
max_truncerr
log_norm_delta_abs
reversibility_density_drift
```

Required conclusion:

- Does the one-site density trend across `n = 3,4,5,6` look size-stable?
- Does the D-improvement trend seen in `3x3` remain qualitatively consistent in larger PBC global density diagnostics?
- Where does increasing `D` saturate or stop helping?

## Campaign C: 7x7 Capacity And Short Dynamics

- [ ] **Step 1: Run 7x7 capacity boundary**

Run:

```bash
SQUAREPXP_LARGERD_N=7 \
SQUAREPXP_LARGERD_DT=0.01 \
SQUAREPXP_LARGERD_D=1 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.0 \
SQUAREPXP_LARGERD_EXACT_FINITE=false \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/7x7-capacity.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/7x7-capacity.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected: capacity artifact files exist, or the runtime/memory boundary is clearly recorded.

- [ ] **Step 2: Run short 7x7 dynamics if capacity is practical**

Run only if capacity is practical:

```bash
SQUAREPXP_LARGERD_N=7 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.04 \
SQUAREPXP_LARGERD_EXACT_FINITE=false \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/7x7-global-short.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/7x7-global-short.csv \
julia --project=. scripts/pxp_larger_d_ed_benchmark.jl
```

Expected: short global artifact files exist, or the runtime/memory boundary is clearly recorded.

- [ ] **Step 3: Analyze 7x7 feasibility**

Record:

```text
ed_basis_dimension
ed_constrained_dimension
ed_group_order
ed_runtime_seconds
ipeps_runtime_seconds
memory/runtime observations
```

Required conclusion:

- Is `7x7` feasible for repeated sweeps?
- If dynamics ran, are conclusions limited to global symmetric PBC observables?

## Optional Energy-Level Extension

The current M3 JSON/CSV schema does not expose ED energy levels. If the next agent finds an existing spectral helper, use it. If not, treat energy levels as a separate implementation task.

- [ ] **Step 1: Search for existing spectral helpers**

Run:

```bash
rg -n "eigen|spectrum|energy level|eigvals|eigs" src test scripts
```

- [ ] **Step 2: If no helper exists, record postponement**

Record in the results note:

```text
Energy-level comparison was not run because the current M3 ED benchmark reports
time-dependent density and return probability, not ED spectra. Adding spectra
requires a separate tested ED spectral-report path.
```

- [ ] **Step 3: If implementing spectra is explicitly requested**

Use TDD and add a separate small ED spectral path for sizes where exact diagonalization is feasible. Do not mix that source change into the benchmark campaign without tests.

## Results Note

- [ ] **Step 1: Create results note**

Create:

```text
docs/superpowers/notes/2026-05-17-m3-systematic-larger-d-results.md
```

The note must include:

1. Exact commands run.
2. Runtime/threading summary.
3. Artifact paths.
4. `3x3` exact-finite D-sweep conclusion.
5. `3x3` longer-time conclusion.
6. Size sweep conclusion across `3x3`, `4x4`, `5x5`, `6x6`, and any feasible `7x7` run.
7. Whether larger `D` improves, saturates, or worsens.
8. Whether smaller `dt` changes the conclusion.
9. Energy-level status: run, postponed, or requiring new implementation.
10. Clear valid and invalid conclusions.
11. Recommended next campaign.

- [ ] **Step 2: Optionally create compact machine-readable summaries**

Create only if useful and reasonably small:

```text
artifacts/m3-systematic/summary.csv
artifacts/m3-systematic/summary.json
```

Do not commit huge raw artifacts unless explicitly requested.

## Final Verification And Response

- [ ] **Step 1: Check final state**

Run:

```bash
git status --short --branch
git log --oneline --decorate -8
```

- [ ] **Step 2: Commit lightweight documentation**

If a results note or compact summary was created and is small enough to track:

```bash
git add docs/superpowers/notes/2026-05-17-m3-systematic-larger-d-results.md
git add artifacts/m3-systematic/summary.csv artifacts/m3-systematic/summary.json 2>/dev/null || true
git commit -m "docs: record systematic M3 larger-D benchmark results"
```

- [ ] **Step 3: Final response**

Report:

- Branch name and final commit SHA if committed.
- Which campaigns actually ran.
- Longest `3x3` time reached.
- Largest `D` tested for `3x3`.
- Largest ED size actually run.
- Whether increasing `D` improved exact finite `3x3` agreement over longer dynamics.
- Whether the one-site density trend is stable across `3x3` through larger PBC sizes.
- Whether `7x7` is feasible for repeated sweeps.
- Energy-level status.
- Runtime or memory bottlenecks.
- Artifact paths.
- Postponed work and why.
