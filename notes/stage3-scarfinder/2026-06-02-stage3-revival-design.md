# Stage-3 staggered-mag revival search — design (2026-06-02)

Verified by probes. Fully ADDITIVE and backward compatible: no changes to
ScarFinderResult/Iteration/CandidateScore, CSV/JSON, CTM paths, or
IPEPSEvolution/Observables. Owner metric = staggered-magnetization REVIVAL.

## Key fact
The series `s(t_i)=density_even-density_odd` is already in
`result.iterations[i].observables` (a SimpleObservableSummary) — no struct change.
Only `s(t_0)` (pre-evolution) is missing; the search layer supplies it via
`measure_simple(seed)`.

## API to add
- `src/ScarFinder.jl`:
  - `staggered_magnetization_series(result; initial=nothing)` (after `_imbalance` L531):
    `[it.observables.density_even-it.observables.density_odd for it in result.iterations]`,
    prepend `_imbalance(initial)` (or Float64) when given. SIGNED series.
  - `revival_quality(series; collapse_frac=0.5)` (after `_revival_strength` L537):
    s0=abs(series[1]); thresh=collapse_frac*s0; i_c=first index with |s|<=thresh;
    return `max(|series[i_c+1:end]|)/s0`, else 0.0 (empty / s0==0 / no collapse /
    collapse only at last sample). Verified: revival 0.9 > decay 0.35; edge cases 0.
  - `RevivalObjective` (L233-245): add `collapse_frac::Float64` field + kwarg
    (default 0.5, validate [0,1)). Keep positional `observable,weight`. KEEP
    `_revival_strength` + per-iteration `_score_value` path UNCHANGED.
  - `_objective_parameters` (L582): add `:revival_collapse_frac`.
  - `revival_score(result, objective::RevivalObjective; initial=nothing)`: requires
    `observable===:sublattice_imbalance`; `revival_quality(series; collapse_frac)`.
  - `scarfinder_search(state_family, params; objective=RevivalObjective(),
    baseline=nothing, max_blockade_violation=Inf, max_truncerr=Inf, scarfinder_kwargs...)`:
    per `(label,builder)` run `scarfinder!(builder(), params; objective, ...)` UNCHANGED,
    score `revival_score(result, rev_obj; initial=measure_simple(state))`; baseline =
    Néel via `checkerboard_square_ipeps(family_cell; excited_on=:even, maxdim=params.trotter.maxdim)`
    scored once; return `ScarFinderSearchResult` with rows sorted by descending
    `margin = revival - baseline_revival`, feasible first. `margin>0` ⇔ better than Néel.
  - structs `SearchCandidateRow`(label,revival,margin,max_blockade,max_truncerr,feasible),
    `ScarFinderSearchResult`(baseline_revival,rows,objective,params), `best_margin`,
    `angle_state_family(cell; thetas, maxdim=1, symmetric=true)`, `_family_cell`.
- `src/SquareIPEPS.jl`: `sublattice_angle_square_ipeps(cell; theta_even,
  theta_odd=pi/2-theta_even, maxdim=1)` — `state_at(c)=[cos θ, sin θ]` (idx1=Rydberg,
  density=cos²θ); θ_even=0,θ_odd=π/2 == Néel (verified imb=1.0). Reuses private
  `_product_square_ipeps`. Blockade-safe for θ_even∈[0,π/4].
- Re-export all new public symbols from `src/SquarePXPDynamics.jl`; every new export
  needs a docstring (test_public_docs).

## Tests (test_scarfinder.jl; keep 4x4, maxdim<=2, iters<=8 — heavy runs time out)
revival_quality semantics incl. "synthetic revival outranks monotone decay";
series shape/length; RevivalObjective backward compat + collapse_frac;
sublattice_angle density convention (θ=0 == checkerboard); scarfinder_search e2e
(rows sorted, baseline_revival present, θ=0 margin≈0); deterministic
"better-than-Néel" via hand-built series; reject non-revival objective.

## Owner decisions / risks (flagged)
1. Whether a real post-collapse revival is observable at affordable (cell, D, time)
   is UNVERIFIED physics — simple-update truncation at small D may damp it. Pick
   params known to show a peak, or accept margin measured on decay/early-return.
2. collapse_frac default 0.5 is a free parameter (threaded through RevivalObjective).
3. "tallest post-collapse peak" vs strictly "first return peak" — 3-line change if first.
4. Series sampled every projection_time; set fine enough to resolve a peak.
5. Signed series + abs in scorer counts return to either polarity (sublattice swap).
6. revival_score/scarfinder_search reject :density (no collapse-revive meaning).
