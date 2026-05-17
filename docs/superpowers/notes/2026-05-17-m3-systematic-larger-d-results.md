# M3 Systematic Larger-D Results

Date: 2026-05-17

## Scope

All runs used the all-zero/all-down initial state and translationally symmetric
PBC dynamics only. The direct ED observable is the translation-invariant
one-site excitation density. For `3x3`, exact finite iPEPS contraction was also
recorded. For larger sizes, `density_simple` is treated as an iPEPS diagnostic,
not an exact finite contraction.

Energy-level comparison was not run because the current M3 ED benchmark reports
time-dependent density and return probability, not ED spectra. Searching
`src`, `test`, and `scripts` found no ED spectral helper; spectra require a
separate tested ED spectral-report path.

## Threading And Preflight

Initial runtime check:

```bash
julia --project=. -e 'using LinearAlgebra; println("BLAS threads = ", BLAS.get_num_threads()); println("Julia threads = ", Threads.nthreads())'
```

Output: BLAS threads `32`, Julia threads `1`.

For single ED-focused runs, the newer git and explicit BLAS setting were used:

```bash
PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH
OPENBLAS_NUM_THREADS=42
```

Check output under that environment: BLAS threads `42`, Julia threads `1`.

Preflight:

```bash
julia --project=. test/runtests.jl test_pxp_larger_d_ed_benchmark.jl
PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH \
  julia --project=. test/runtests.jl test_pxp_ed_benchmark.jl test_pxp_validation.jl
```

Results:

- `test_pxp_larger_d_ed_benchmark.jl`: `58/58` passed in `2m00.4s`.
- ED/validation preflight initially failed because the system git rejected
  `git -C`; rerunning with the newer git passed `206/206` in `2m11.2s`.

## Commands And Artifacts

Campaign A was first started as one process, then stopped after observing that
it used one core. It was rerun as per-case workers with durable per-case
artifacts in `artifacts/m3-systematic/parts/` and logs in
`logs/m3-systematic/`.

```bash
mkdir -p artifacts/m3-systematic/parts logs/m3-systematic

printf '%s\n' \
  '0.02 1' '0.02 2' '0.02 3' '0.02 4' \
  '0.01 1' '0.01 2' '0.01 3' '0.01 4' \
  '0.005 1' '0.005 2' '0.005 3' '0.005 4' |
xargs -n 2 -P 12 bash -c 'dt="$0"; D="$1"; stem="3x3-t020-dt${dt}-D${D}";
  PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH \
  OPENBLAS_NUM_THREADS=1 \
  SQUAREPXP_LARGERD_N=3 \
  SQUAREPXP_LARGERD_DT="$dt" \
  SQUAREPXP_LARGERD_D="$D" \
  SQUAREPXP_LARGERD_CUTOFF=1e-12 \
  SQUAREPXP_LARGERD_TOTAL_TIME=0.20 \
  SQUAREPXP_LARGERD_EXACT_FINITE=true \
  SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
  SQUAREPXP_LARGERD_JSON="artifacts/m3-systematic/parts/${stem}.json" \
  SQUAREPXP_LARGERD_CSV="artifacts/m3-systematic/parts/${stem}.csv" \
  /usr/bin/time -p julia --project=. scripts/pxp_larger_d_ed_benchmark.jl \
    > "logs/m3-systematic/${stem}.log" 2>&1'
```

The optional `3x3`, `total_time=0.50` extension used the same per-case command
with `SQUAREPXP_LARGERD_TOTAL_TIME=0.50` and `OPENBLAS_NUM_THREADS=3`.

Per-case outputs were merged into:

- `artifacts/m3-systematic/3x3-exact-long.json`
- `artifacts/m3-systematic/3x3-exact-long.csv`
- `artifacts/m3-systematic/3x3-exact-longer.json`
- `artifacts/m3-systematic/3x3-exact-longer.csv`

Campaign B was run as one ED-focused process:

```bash
PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH \
OPENBLAS_NUM_THREADS=42 \
SQUAREPXP_LARGERD_N=3,4,5,6 \
SQUAREPXP_LARGERD_DT=0.02,0.01 \
SQUAREPXP_LARGERD_D=1,2,3 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.10 \
SQUAREPXP_LARGERD_EXACT_FINITE=true \
SQUAREPXP_LARGERD_EXACT_FINITE_MAX_SITES=9 \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/size-sweep-3to6.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/size-sweep-3to6.csv \
/usr/bin/time -p julia --project=. scripts/pxp_larger_d_ed_benchmark.jl \
  > logs/m3-systematic/size-sweep-3to6.log 2>&1
```

Runtime: `real 1639.71`, `user 12997.77`, `sys 1397.79`.

Campaign C capacity probe:

```bash
PATH=/usr/share/atom/resources/app/node_modules/dugite/git/bin:$PATH \
OPENBLAS_NUM_THREADS=42 \
SQUAREPXP_LARGERD_N=7 \
SQUAREPXP_LARGERD_DT=0.01 \
SQUAREPXP_LARGERD_D=1 \
SQUAREPXP_LARGERD_CUTOFF=1e-12 \
SQUAREPXP_LARGERD_TOTAL_TIME=0.0 \
SQUAREPXP_LARGERD_EXACT_FINITE=false \
SQUAREPXP_LARGERD_JSON=artifacts/m3-systematic/7x7-capacity.json \
SQUAREPXP_LARGERD_CSV=artifacts/m3-systematic/7x7-capacity.csv \
/usr/bin/time -p julia --project=. scripts/pxp_larger_d_ed_benchmark.jl \
  > logs/m3-systematic/7x7-capacity.log 2>&1
```

Observation: after about 20 minutes, the process was still in a single-core
setup path, using about `1.8 GB` RSS and writing no artifact. It was stopped.
Short `7x7` dynamics were not run.

Compact summaries:

- `artifacts/m3-systematic/summary.csv`
- `artifacts/m3-systematic/summary.json`

## 3x3 Exact-Finite Trends

At `total_time=0.20`, exact-finite density errors:

| dt | D=1 | D=2 | D=3 | D=4 |
| --- | ---: | ---: | ---: | ---: |
| 0.02 | `-2.109e-2` | `-2.496e-5` | `-6.448e-4` | `-9.702e-3` |
| 0.01 | `-2.634e-2` | `-1.883e-5` | `-1.973e-3` | `-1.156e-2` |
| 0.005 | `-3.024e-2` | `-1.743e-5` | `-3.541e-3` | `-4.986e-3` |

At `total_time=0.50`, exact-finite density errors:

| dt | D=1 | D=2 | D=3 | D=4 |
| --- | ---: | ---: | ---: | ---: |
| 0.02 | `-1.520e-1` | `-1.210e-3` | `-9.700e-3` | `-3.072e-2` |
| 0.01 | `-1.586e-1` | `-1.204e-3` | `-1.714e-2` | `-1.673e-2` |
| 0.005 | `-1.628e-1` | `-1.205e-3` | `-2.412e-2` | `-6.633e-3` |

Conclusion: increasing from `D=1` to `D=2` strongly improves exact-finite
`3x3` density and return-probability agreement at both `0.20` and `0.50`.
Increasing beyond `D=2` is not monotonic in this campaign. `D=3` and `D=4`
remain better than `D=1` in exact-finite density, but are usually worse than
`D=2`. The improvement is therefore saturated/degraded beyond `D=2`, not a
monotone D-scaling result.

Smaller `dt` does not change the qualitative conclusion that `D=2` is best in
the exact-finite density comparison. For `D=1`, smaller `dt` did not improve
the endpoint density error in these first-order serial-update runs.

Return probability tells the same broad story: at `total_time=0.50`, `D=2`
has return-probability error about `4.8e-3`, while `D=1` is about `0.80` to
`0.89`; `D=3`/`D=4` are intermediate and not monotone.

Important diagnostic: D>1 `density_simple` does not track the exact finite
density and should not be used as exact finite truth. The exact-finite iPEPS
density field is the relevant `3x3` comparison.

## Size Sweep

At `total_time=0.10`, ED excitation density is size-stable:

| n | dt=0.02 ED density | dt=0.01 ED density |
| --- | ---: | ---: |
| 3 | `0.0098352277` | `0.0098352277` |
| 4 | `0.0098356158` | `0.0098356158` |
| 5 | `0.0098356148` | `0.0098356148` |
| 6 | `0.0098356148` | `0.0098356148` |

The larger-size iPEPS simple-density diagnostic does not show the same
D-improvement as the `3x3` exact-finite observable. For example at `dt=0.02`:

| n | D=1 simple error | D=2 simple error | D=3 simple error |
| --- | ---: | ---: | ---: |
| 3 | `-1.578e-3` | `-5.984e-3` | `-6.173e-3` |
| 4 | `-1.579e-3` | `-6.352e-3` | `-6.704e-3` |
| 5 | `-1.579e-3` | `-6.565e-3` | `-6.723e-3` |
| 6 | `-1.579e-3` | `-6.703e-3` | `-6.202e-3` |

The valid conclusion is that ED one-site excitation density is stable across
`3x3` through `6x6` at this short time. It is not valid to infer from the
larger-size simple-density diagnostic that larger D worsens the actual finite
iPEPS state; exact finite contraction is only populated for `3x3` here.

The optional `D=4`, `total_time=0.20` size sweep was not run. The completed
single-process `D=1..3`, `3x3..6x6`, `total_time=0.10` sweep took about
27 minutes; extending both D and time would be substantially more expensive.

## 7x7 Feasibility

The `7x7` capacity probe did not complete within about 20 minutes and produced
no JSON/CSV artifact. Observed resource use during setup was about one CPU core
and `1.8 GB` RSS, indicating a non-BLAS setup bottleneck before any useful
capacity row was emitted.

Conclusion: `7x7` is not feasible for repeated sweeps with the current helper
and campaign settings. No short `7x7` dynamics were run.

## Valid Conclusions

- `3x3` exact-finite density supports a strong improvement from `D=1` to `D=2`.
- `D=2` is the best exact-finite `3x3` setting tested at `total_time=0.20` and
  `0.50`.
- `D=3` and `D=4` do not give monotone improvement over `D=2` in this campaign.
- ED one-site excitation density is stable across `n=3,4,5,6` at
  `total_time=0.10`.
- `density_simple` remains a diagnostic for D>1 and should be interpreted
  separately from exact finite density.
- Current `7x7` setup is not practical for repeated sweeps.

## Not Valid

- Claiming publication-grade physics from these runs.
- Treating D>1 `density_simple` as an exact finite observable.
- Inferring spatially resolved ED fields from the symmetry-reduced ED basis.
- Claiming monotone improvement with D beyond `D=2`.
- Claiming `7x7` dynamics feasibility from the stopped capacity probe.

## Recommended Next Campaign

Superseding update from later on 2026-05-17: for D>1 local-density comparisons,
use CTM/environment observables rather than simple/local observables. A direct
CTM probe at `3x3`, `t = 0.02`, `chi = 2` showed D=2 CTM density matching exact
finite density while `density_simple` did not. The active next campaign is
therefore CTM observable throughput and warmed iPEPS+CTM timing, not larger ED.

Use `3x3` exact finite contraction and CTM density to investigate why `D=2`
dominates `D=3` and `D=4`, with particular attention to log-normalization
growth, split/truncation diagnostics, and CTM finite-chi diagnostics. For larger
sizes, add a reusable ED-reference cache or snapshotable iPEPS-observable path
before repeating expensive sweeps.
