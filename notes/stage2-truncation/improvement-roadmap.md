# Stage-2 improvement roadmap

Distilled from a 6-lens brainstorm (`2026-06-02-pxp-improvement-brainstorm.md`)
and the trajectory correction (`../methodology/revival-validation.md`). Three
levers; only env attacks the truncation directly.

## Current status (corrected)

By the n(t) trajectory the iPEPS shows CLEAN MONOTONE D-convergence up to D=4
(RMS D2 2.6e-2 → D4 9.8e-3) — there is NO demonstrated ceiling and NO bad pocket.
The earlier "mean-field environment ceiling" was a t=2.6 endpoint artifact. The
evolution is exact-finite-correct; CTM flatters (use the exact oracle); regauging
(`canonicalize_simple!`) is a no-op for dynamics (gauge transform).

## Top priorities

1. **Boundary-MPS exact contractor — DONE (2026-06-03), perf follow-up open.**
   Built in `src/Observables.jl`: a double-layer column-ring boundary contraction
   (`exact_density_finite(...; method = :boundary)`), exact, with bounded memory
   (closes physical legs locally → no `2^N` factor; boundary kept factorized).
   Validated == the dense path to ~1e-9 on 3×3 and 4×4 (D=2,3,4), including the
   evolved-checkerboard benchmark; reproduces the recorded trajectory RMS
   (D2 2.62e-2, D3 1.82e-2, D4 9.79e-3). `method=:auto` uses dense for ≤16 sites,
   boundary for larger.
   - **Why it was needed (confirmed):** the dense single-layer contraction is
     memory-pathological for *entangled* D≥5 on a 4×4 torus (~`2^(N/2)·D^(2Lx)`;
     observed **~244 GB** RSS at D=5 on the shared host — i.e. the OOM is real, not
     just on small nodes). The boundary path stays bounded.
   - **Open perf follow-up:** the boundary path is slow at full-rank D≥5 (~one
     `D^12` seam SVD per sweep, and `_boundary_density_finite` does N+1 sweeps per
     density). Cheap wins: per-site **environment caching** (N+1 sweeps → ~2) and
     proper **periodic-MPS seam canonicalization** (avoid the bloated seam tensor).
     Needed to make the full D=5,6 trajectory fast; a background run is in flight.
2. **Re-settle rel_floor + the regression test by TRAJECTORY** (not t=2.6). See
   `rel-floor.md`. Move the `test_d_convergence.jl` revival gate off the endpoint
   to a max/RMS trajectory metric.
3. **Environment-aware (cluster/full) update** — truncate against a real bond
   environment instead of the single-site λ². Hooks (`StarSimpleUpdate.jl`): the
   mean-field assumption enters at `_absorb_star_weights` / `_qr_reduce_leaves`;
   reweight σ in `_split_reduced_theta` (:334/:342). The 2×2-cluster update (#5)
   contracts the nearest loops exactly and sidesteps environment inversion. Now
   framed as "tighten convergence / remove the revival-rise lag", not "break a
   wall". NOTE: the norm environment ⟨ψ|ψ⟩ is PSD even in real time — the
   "indefinite metric" worry is a numerical-regularization issue, not fundamental.

## Lower-priority levers (from the brainstorm)

- Exact 2-color checkerboard schedule + order-4 Trotter (cuts Trotter error only;
  see `../conventions/pxp-model-states-trotter.md`).
- ED-free validators for sizes beyond exact contraction: reverse-evolution echo
  R(t), energy-drift E(t)−E(0) (reuse `reverse_evolve!` / `exact_pxp_energy_density_finite`).
- isoTNS / BP-gauge / constraint-resolved bond basis (ansatz-level, larger effort).

Full ranked list + the ruled-out-flags analysis: `2026-06-02-pxp-improvement-brainstorm.md`.
