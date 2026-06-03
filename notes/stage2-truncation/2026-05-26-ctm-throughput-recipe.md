# CTM throughput recipe

Purpose: pin down a default threading layout for `measure_ctm` / PEPSKit CTMRG
on this server so audit campaigns (Phases 2-5 of the ScarFinder reliability
plan) are tractable. Replaces the earlier ad-hoc guidance (the retired
`memory/short_term/next_steps.md`) once a recommended row is selected.

## Why this exists

PEPSKit CTMRG parallelism is **Julia-thread based**: setting `OPENBLAS_NUM_THREADS`
alone is insufficient. There are four in-process knobs that interact:

| Knob                          | Affects                                  | Set via                                                  |
|-------------------------------|------------------------------------------|----------------------------------------------------------|
| `JULIA_NUM_THREADS`           | PEPSKit `dtmap` parallel regions         | shell env (must be set before Julia starts)              |
| BLAS thread count             | LAPACK/BLAS work inside tensor ops       | `configure_ctm_threading!(blas_threads=...)`             |
| Strided thread count          | TensorOperations / Strided tensor work   | `configure_ctm_threading!(strided_threads=...)`          |
| Strided threaded `mul`        | Julia-threaded splitting of `*`          | `configure_ctm_threading!(strided_threaded_mul=...)`     |
| PEPSKit scheduler             | `dtmap` schedule (`:default`/`:dynamic`) | `configure_ctm_threading!(pepskit_scheduler=...)`        |

`JULIA_NUM_THREADS` cannot be changed inside a running session, so the timing
matrix script varies the other three within one session, and is re-invoked
once per thread count.

## Script

`scripts/run_ctm_timing_matrix.jl` builds a representative iPEPS state
(`Lx x Ly` cell, `D` bond dim, optional short real-time PXP evolution at
`evolve_time`), then for each `(blas_threads, strided_threads, threaded_mul,
scheduler, chi)` combination performs one warm-up `measure_ctm` call followed
by `reps` timed calls. Output is written to JSON and CSV; passing
`SQUAREPXP_TIMING_APPEND=true` with a shared CSV path concatenates rows from
multiple invocations.

### Environment variables

| Variable                               | Default                  | Notes                                       |
|----------------------------------------|--------------------------|---------------------------------------------|
| `SQUAREPXP_TIMING_CELL_LX`             | `3`                      | unit cell width                             |
| `SQUAREPXP_TIMING_CELL_LY`             | `3`                      | unit cell height                            |
| `SQUAREPXP_TIMING_D`                   | `2`                      | iPEPS bond dimension                        |
| `SQUAREPXP_TIMING_DT`                  | `0.02`                   | Trotter step                                |
| `SQUAREPXP_TIMING_EVOLVE_TIME`         | `0.02`                   | real-time PXP evolution before timing       |
| `SQUAREPXP_TIMING_CHI_VALUES`          | `16`                     | comma-separated CTMRG bond dims             |
| `SQUAREPXP_TIMING_BLAS_VALUES`         | `1`                      | comma-separated BLAS thread counts          |
| `SQUAREPXP_TIMING_STRIDED_VALUES`      | `Threads.nthreads()`     | comma-separated Strided thread counts       |
| `SQUAREPXP_TIMING_THREADED_MUL_VALUES` | `false,true`             | comma-separated booleans                    |
| `SQUAREPXP_TIMING_SCHEDULER_VALUES`    | `default,dynamic`        | comma-separated scheduler symbols           |
| `SQUAREPXP_TIMING_REPS`                | `3`                      | timed reps per config (excludes warmup)     |
| `SQUAREPXP_TIMING_LABEL`               | `""`                     | free-text tag, e.g. `j42-blas1`             |
| `SQUAREPXP_TIMING_OUTPUT_CSV`          | `artifacts/ctm_timing_matrix.csv` |                                       |
| `SQUAREPXP_TIMING_OUTPUT_JSON`         | `artifacts/ctm_timing_matrix.json` |                                       |
| `SQUAREPXP_TIMING_APPEND`              | `false`                  | append to CSV instead of overwriting        |

## Recommended sweep

Run the script three times to populate the matrix:

```bash
# 1) Single-threaded baseline (cache miss/compile floor)
SQUAREPXP_TIMING_LABEL=j1 \
SQUAREPXP_TIMING_OUTPUT_CSV=artifacts/ctm_timing_matrix.csv \
SQUAREPXP_TIMING_OUTPUT_JSON=artifacts/ctm_timing_matrix-j1.json \
SQUAREPXP_TIMING_APPEND=false \
JULIA_NUM_THREADS=1 julia --project=. scripts/run_ctm_timing_matrix.jl

# 2) Mid-range
SQUAREPXP_TIMING_LABEL=j8 \
SQUAREPXP_TIMING_OUTPUT_CSV=artifacts/ctm_timing_matrix.csv \
SQUAREPXP_TIMING_OUTPUT_JSON=artifacts/ctm_timing_matrix-j8.json \
SQUAREPXP_TIMING_APPEND=true \
JULIA_NUM_THREADS=8 julia --project=. scripts/run_ctm_timing_matrix.jl

# 3) Full-server target
SQUAREPXP_TIMING_LABEL=j42 \
SQUAREPXP_TIMING_OUTPUT_CSV=artifacts/ctm_timing_matrix.csv \
SQUAREPXP_TIMING_OUTPUT_JSON=artifacts/ctm_timing_matrix-j42.json \
SQUAREPXP_TIMING_APPEND=true \
JULIA_NUM_THREADS=42 julia --project=. scripts/run_ctm_timing_matrix.jl
```

For a serious calibration also sweep `SQUAREPXP_TIMING_CHI_VALUES=16,32,64`,
since the optimal split between Strided and PEPSKit threads is `chi`-dependent.

## How to read the output

For each row, compare `median_seconds` across the threading dimensions while
holding `chi`, `D`, and `cell_*` fixed. The winning recipe is the one with
the lowest `median_seconds` *and* `ctm_converged = true`. If a row fails to
converge under a given threading config, do not select it even if it is
faster — divergence under threaded `mul` would mean the throughput "win" came
from skipping work. Drop it and pick the next-best stable config.

Once chosen, record the recipe in this file under "Selected recipe" below and
update `memory/short_term/next_steps.md` / `current_state.md` to point at it.

## Smoke validation

A `JULIA_NUM_THREADS=4` smoke at `3x3, D=2, chi=16, two configs` is enough to
confirm the script runs end-to-end. The full server matrix (`j42` with
`chi=16,32,64`) takes hours and should be scheduled separately.

## Selected recipe

Pending — fill in after the matrix has been run. The current memory-recorded
guess (used only to bootstrap the calibration) is:

```text
JULIA_NUM_THREADS=42
SQUAREPXP_CTM_BLAS_THREADS=1
SQUAREPXP_CTM_STRIDED_THREADS=42
SQUAREPXP_CTM_STRIDED_THREADED_MUL=true
SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic
```

This block is the hypothesis being tested, not the answer.
