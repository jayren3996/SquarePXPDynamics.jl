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
   - **Open perf follow-up (gates the D=5,6 trajectory).** At full-rank D≥5 the
     boundary path is slow. Measured (4×4, D=5, one sweep): rows 3–4 are cheap
     (~6–8 s) but **row 2 dominates** — zip 46 s + recompress 116 s ≈ 162 s — because
     that is where the full-`D²` boundary first forms and a ~`D^12` (≈3.9 GB) seam
     intermediate appears in BOTH the column-1 zip and the seam two-site SVD.
     `_boundary_density_finite` then repeats a full sweep N+1 times (Z + one per
     site) ⇒ ~46 min/density at D=5. **Two optimizations were tried and FAILED for
     the entangled regime (2026-06-03):**
     - *Thread-parallelism over the N sweeps* — each concurrent sweep holds its own
       ~3.9 GB seam, so `-t 16` blew RSS to ~146 GB and stalled. (Kept the
       `Threads.@threads`; it only helps when per-sweep memory is small.)
     - *Per-row environment caching* (one bottom + one top sweep, then each `Z_c` a
       cheap ring combine of cached environments). Correct (== dense to ~1e-9 on
       4×4 D≤4) and N+1→2 sweeps in principle, BUT the combine does no
       recompression, so on an *entangled* state the bottom+row+top ring contraction
       forms ~`χ^4` intermediates and **OOMs at the D=5 revival peak** (t=2.6); even
       the product point was not faster (~700 s). Reverted.
     - **Root cause / real fix:** the per-sweep cost and these OOMs all trace to the
       same thing — the *entangled* double-layer boundary genuinely has large bonds
       (Schmidt rank up to ~`D^4`), and a `D²`-bonded seam contraction is ~`D^12`.
       A practical contractor needs recompression *interleaved with* the
       environment combine (a proper boundary-MPS-with-environments / variational
       contraction), or a different decomposition — non-trivial, still open.
     - **For now:** the committed boundary path is correct + memory-bounded but slow
       (~46 min/density at D=5); the dense path is fine for D≤6 *low-entanglement*
       but OOMs (~244 GB) at the entangled peak. So the full D=5,6 *trajectory*
       (which must reach the entangled revival) is still gated on the faster
       contraction above.
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
   **STATE SNAPSHOT (2026-06-04): `2026-06-04-ctm-aware-evolution-state.md`** —
   read this for exactly where the CTM-aware evolution stands. TL;DR: the per-bond
   CTM truncation primitive is validated (26× better than bare SVD); the full-star
   geometry passes horizontally but the vertical bond is blocked on an untested
   `rotl90(env)` fix; per-star CTMRG is trajectory-infeasible; the open fork is
   which loop-carrying engine to build (exact-finite-cluster env vs NTU patch), and
   the blowup is CUMULATIVE so only a full-trajectory run can demonstrate the fix.

## Lower-priority levers (from the brainstorm)

- Exact 2-color checkerboard schedule + order-4 Trotter (cuts Trotter error only;
  see `../conventions/pxp-model-states-trotter.md`).
- ED-free validators for sizes beyond exact contraction: reverse-evolution echo
  R(t), energy-drift E(t)−E(0) (reuse `reverse_evolve!` / `exact_pxp_energy_density_finite`).
- isoTNS / BP-gauge / constraint-resolved bond basis (ansatz-level, larger effort).

Full ranked list + the ruled-out-flags analysis: `2026-06-02-pxp-improvement-brainstorm.md`.
