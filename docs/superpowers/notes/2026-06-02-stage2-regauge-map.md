# Stage-2 code map + diagnosis: bond truncation / gauge / environment (2026-06-02)

Read-only map of the evolution truncation+gauge pipeline, to ground the Stage-2
"good bond truncation + recover good condition with proper regauging" goal.

## What the evolution actually does (it is correct simple-update)
`project_star!` (src/StarSimpleUpdate.jl:453) runs a QR-reduced five-site star
update per center:
1. `_absorb_star_weights` (:250) — each leaf absorbs the FULL λ on its
   center-link **and** its 3 external links (`absorb_link_weight`,
   SquareIPEPS.jl:386 → `T * diag(λ)`).
2. `_qr_reduce_leaves` (:268) — `factorize(...; ortho="left", which_decomp="qr")`
   splits each absorbed leaf into Q (external/environment legs) · R (physical).
3. gate · center · R-factors → reduced core `theta` (:480–487).
4. `_split_reduced_theta` (:307) — sequential SVDs along `split_order`; each:
   - `_svd_with_rel_floor` (:235) applies maxdim/cutoff + the rel_floor relative
     SV condition floor (keep σ ≥ rel_floor·σ_max; caps cond at 1/rel_floor).
   - `_normalized_link_weight_from_singular_values` (:186) → `scale=‖σ‖₂`,
     `λ=σ/scale`; core re-scaled `core = V * scale` (:342) to carry the norm.
5. `_reconstruct_leaf` (:364) — `Q·active`, then `deabsorb_link_weight` the 3
   external λ back out (SquareIPEPS.jl:402 → `T * diag(inv(λ))`, with atol guard).
6. `_commit_star_update!` (:392) — write Γ tensors + new λ, track log-norm.

Weights live in `SquareIPEPSState.link_weights::Dict{BondKey,Vector{Float64}}`
(SquareIPEPS.jl:37), normalized to unit ‖·‖₂, init `[1,0,…,0]`. Hot-path read
`_link_weight_view` (:283).

## Key diagnosis
- Absorbing the **full** λ on environment (external) legs is the textbook
  simple-update mean-field environment (super-orthogonalization: λ² is the
  approximate environment). This is NOT a bug; the evolution is correct.
- Measurement ALREADY does the correct symmetric √λ readout:
  PEPSKitMeasurements.jl:751 splits each bond λ as √λ⊗√λ onto endpoints before
  CTMRG. So the Γ-λ→PEPS conversion for observables is right.
- There is **no regauging code** in evolution (`gauge` field hardcoded `:simple`,
  SquareIPEPS.jl:39); no `regauge`/`canonical`/full-update path exists.
- The only conditioning knob today is `rel_floor` (default 1e-4) — a crude,
  regime-dependent cap that has a bad pocket at the revival.

=> The D=4-worse-than-D=2/3 non-monotonicity AT THE REVIVAL is the inherent
mean-field limitation of simple update on a loopy lattice (chi-sweep already
proved it is an evolution, not a measurement, limit). Closing it requires a
REAL environment inside the evolution loop, not a better knob.

## Stage-2 options, cheapest → most-correct
1. **rel_floor auto / Vidal canonicalization polish** (small): re-canonicalize
   the Γ-λ form between sweeps so λ are true Schmidt values (gauge fix on the
   tree-approx), and/or make rel_floor adaptive. Cheap, may smooth conditioning,
   will NOT remove the mean-field error at the revival.
2. **Loop/cluster-corrected simple update** (medium): include nearest-neighbour
   loop tensors (e.g. simple-update with loop gas / "gauge fixing" à la
   Evenbly, or nth-order cluster) — partial environment, no full CTM.
3. **Full update** (large): use the CTMRG environment (already built for
   measurement) inside the truncation — environment-aware SVD of the bond.
   This is the principled fix for the revival non-monotonicity but is a
   research-grade feature (env construction per bond/sweep, cost, stability).

## Exact hook points for an environment-aware truncation
- `project_star!`:477 — precompute env context before absorption.
- `_qr_reduce_leaves`:276 — where the mean-field (orthogonality) assumption
  enters; swap in an environment-weighted reduction.
- `_split_reduced_theta`:334 — reweight σ by the environment before
  normalize/truncate (the core of full-update truncation).
- `_split_reduced_theta`:342 — adjust `scale` to conserve norm under env weights.
- `_reconstruct_leaf`:375 — match the deabsorption to the chosen gauge.

CTM environment already available via `pepskit_ctmrg_context` /
`to_pepskit_infinitepeps` (PEPSKitMeasurements.jl) — the full-update path can
reuse it rather than building a new contractor.
