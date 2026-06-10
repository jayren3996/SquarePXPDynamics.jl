# CTM-aware 2D PXP dynamics — status for external review (2026-06-09)

**Audience:** an external reviewer / collaborator. This is a self-contained summary; the blow-by-blow
running log (with every workflow, refutation, and number) is in
`notes/stage2-truncation/2026-06-07-workflow-fanout-review-verdict.md`. Design docs:
`docs/superpowers/plans/2026-06-04-ctm-aware-evolution.md`, `...2026-06-07-ctm-upgrade-workflow.md`.

---

## 1. Goal and setup

`SquarePXPDynamics.jl` does **2D square-lattice PXP real-time dynamics** on a PEPSKit-native iPEPS. Time
evolution is a Trotterized **5-site "star" gate** (center + 4 leaves) applied in a `:serial` order-2
palindrome; each gate application truncates the four center↔leaf virtual bonds. The scientific target is the
**4×4 Néel quench** — a collapse (t≈1.0–1.5) and **revival** (peak t≈2.6) — benchmarked against an **exact
diagonalization (ED) oracle** for the finite 4×4 torus (`artifacts/neel_to_revival_4x4.json`). The observable
is the local density n(t) via `exact_density_finite` (a single-layer finite contraction, ledger- and
scale-invariant).

**The problem we set out to fix:** at bond dimension **D=2** the bare (SVD) truncation makes the 4×4 quench
**catastrophically unstable** — n(t) blows up from ~0.14 to ~0.43 at t=1.6→1.7 and never recovers the revival
(dead by t=2.2). The hypothesis (standard for 2D tensor networks) is that the bare per-bond SVD uses a
**bad mean-field environment** (the severed-loop λ² weights), and that a **loop-carrying bond environment**
(CTM-aware truncation) would let a small bond dimension behave like a larger one.

The truncation primitive is PEPSKit's `bond_truncate(a, b, benv, alg)`, which takes a `{2,2}` bond environment
`benv` (a `TensorMap` on the two bond-end virtual "q-legs", double-layer bra-ket with physical legs traced).
Our work plugs **different environments** into this same primitive behind a flag
(`env = :bare | :exact_cluster | :ntu`), via `scripts/g0/star_env.jl` (`build_star_bond_env`,
`project_star_env!`).

---

## 2. Result that IS solid — "effect A": a loop env stabilizes D=2

We built the **exact-cluster environment**: the full 4×4-torus double-layer contraction with the active bond's
q-legs open (`build_star_bond_env(...; env=:exact_cluster)`). It is the **gold-standard finite-torus bond
environment** (validated: a transpose-oracle value check at 1e-11, no-op-at-maxdim 1e-16, fixgauge robust even
at cond(Z'Z)→∞).

Forking a bare-D2 trajectory at the disease onset (t=1.6) and switching only the truncation environment, the
**exact-cluster env at D=2 removes the blowup and reproduces the entire revival** (`env_window_dt02_to2.8.out`):

| t | env-D2 (exact-cluster) | exact ED | Δ |
|---|------------------------|----------|---|
| 1.8 | 0.2232 | 0.2437 | −0.021 |
| 2.0 | 0.3165 | 0.3374 | −0.021 |
| 2.4 | 0.4647 | 0.4670 | −0.002 |
| **2.6 (peak)** | **0.4928** | **0.4831** | **+0.010** |
| 2.8 | 0.4825 | 0.4623 | +0.020 |

Stable throughout (minkept=2, nfb=0, finite to t=2.8). **A loop-carrying environment turns the
catastrophically-unstable bare-D2 quench into a stable integrator accurate to ~±2% of exact ED across the whole
revival.** This is the core positive finding and we consider it well-established (independently reviewed; the
bare-vs-env control and ledger-bug fix are documented).

## 3. The residual ~2% is bond-dimension error, not Trotter (the lever is D)

A dt-convergence test (`env_dtconv_dt01_to1.8.out`) from the identical fork: env(t=1.8) is **0.223205 at
dt=0.02 and 0.223225 at dt=0.01** (a +2e-5 move). Richardson dt→0 ≈ 0.22323 — still −0.0205 below ED. So the
env converges (in dt) to its own D=2-truncated dynamics, ~2% below true physics; **closing that gap needs more
bond dimension D, not smaller dt.** (Native simple-update at rel_floor=0 confirms the direction: bare-D3 at
t=1.6 = 0.15554 vs ED 0.15558 — D=3 is near-exact where D=2 blew up.)

---

## 4. The obstruction: no feasible *and* accurate environment at the bond dimension that matters

To improve accuracy we need a loop env at **D≥3**. Three environment families, and **each fails a different
requirement**:

| engine | accuracy | feasible at D≥3 | status |
|--------|----------|-----------------|--------|
| **exact-cluster** (full 4×4 torus) | ~2% (gold) | **NO** | reference only |
| **NTU-P6** (local 6-cell patch) | **~14%** | yes (cheap) | validated, but crude |
| **CTMRG** (χ-controlled boundary) | unknown on 4×4 | yes | viable but unvalidatable here |

**(a) exact-cluster is infeasible at D≥3.** The full-torus contraction cost scales ~(D²)^width; a D=3
trajectory ran **32 h / 70 GB and did not finish a single 0.1-time-chunk** before we killed it
(`env_d3_control_to2.0.out`). It is a 4×4-D=2 reference, never a production engine — exactly the design doc's
reason for preferring a *local* environment.

**(b) NTU-P6 (the cheap local patch) is structurally too crude on 4×4.** We audited + value-validated the
4-direction NTU builder (vertical patch transpose-oracle 6e-14, no-op 1e-15; `ntu_vertical_gate.jl`). It is
**5.5× cheaper than exact-cluster (62 min vs 5.7 h) and stable** (removes the blowup, reproduces the revival
*shape* and timing). **But it systematically undershoots by ~14% at the peak** (vs exact-cluster's ~2%) —
`ntu_window_P6_dt02_to2.8.out`. The open-boundary P6 patch is a genuinely cruder environment (frozen-bond
rel-dev 0.949 vs exact, "nearly orthogonal"); the per-step bias accumulates. The two refinement levers are
both closed on a 4×4: the `boundary=:lambda` mean-field dressing is *worse* (it double-counts the SU weights,
660× worse at a frozen bond), and the larger patches P10/P12 **self-alias on the 4×4 wrap** (they need
linear size ≥5). So NTU-P6 is a cheap *stabilizer* with a ~14% accuracy floor that cannot be refined here.

**(c) CTMRG is viable but cannot be validated on the 4×4 testbed.** CTMRG (`leading_boundary` →
`bondenv_fu`) **converges on healthy states** (err~1e-10) and its CTM-aware truncation **beats bare-SVD 9–26×**
at a frozen bond. But every cheap validation route is blocked:
  - CTMRG computes the **infinite-tiling (thermodynamic-limit)** environment; the ED oracle and exact-cluster
    are the **finite 4×4 torus** — *different boundary conditions*, so CTMRG cannot be checked against the 4×4
    ED the way exact-cluster was.
  - It will **not converge on the bare-blown-up D=2 states** (the regime where the env matters), and the only
    *healthy* high-entanglement states need **D≥3, where the exact-cluster reference is infeasible** — so the
    feasible-reference regime and the healthy-high-entanglement regime **never overlap**.
  - The frozen-bond env-comparison is **gauge-confounded** (the CTM environment carries a gauge freedom not
    fixed by the objective; raw-tensor rel-dev to exact-cluster is erratic across χ even though CTMRG
    converged to one fixed point — `ctm_chi_ladder_t1.0.out`), and the single-bond `|dn|` is
    **non-discriminating** (all envs are near-lossless at one frozen bond).

**Net obstruction.** The 4×4 + exact-ED testbed has delivered its science (effect A; the accuracy lever is D),
but it **cannot cheaply validate an environment that is simultaneously accurate and feasible at D≥3.** The only
remaining route to settle CTM accuracy is a full **CTM-aware trajectory** — which is a real engineering build
(the CTMRG env must be amortized: rebuilt once per lattice sweep, not per-bond-per-star; a single CTMRG solve
was 75 s–1461 s at D=2), and which on 4×4 would model the *thermodynamic limit* (no finite-4×4 oracle to check
against — validation would have to be by χ-convergence + finite-size extrapolation + cross-observables).

---

## 5. Specific questions for external review

1. **Validation without an oracle.** What is the accepted standard for validating a CTM-aware *real-time*
   evolution in 2D when there is no ED reference at the target system size? Is χ-convergence + finite-size
   extrapolation + an independent observable (energy / Loschmidt) sufficient, and how is the non-unitary
   env-weighted norm bookkeeping (`add_log_norm!`) handled for energy/Loschmidt (density is a robust ratio;
   these are not)?
2. **Finite vs infinite, and the revival.** Is the t≈2.6 4×4 revival a genuine scar feature or a finite-size
   artifact? The CTMRG (infinite-tiling) vs exact-cluster (finite-torus) environment gap is exactly the
   finite-size signal — but we could not measure it (gauge-confounded). Is there a clean way to quantify it?
3. **A feasible accurate environment.** Is there an environment *between* the (accurate, infeasible)
   full-torus and the (cheap, ~14%) open-boundary P6 — e.g. a boundary-MPS / finite-χ environment for the
   *finite* torus, or a larger/lambda-corrected patch — that is accurate **and** feasible at D≥3 on a system
   that still has an ED oracle? (On 4×4, P10/P12 self-alias; a 5×5–6×6 torus would admit them but loses the
   cheap ED reference.)
4. **Is D=2-plus-good-env even the right target?** Effect A shows a loop env makes D=2 ~2%-accurate. Is the
   physically interesting program "small-D + good environment" (our framing), or should this be D≥3 simple
   update with a feasible environment, or a different ansatz entirely?
5. **NTU boundary closure.** The `:trace` (max-mixed) patch-boundary closure biases toward the bare-SVD
   subspace and undershoots; `:lambda` double-counts the SU weights and is worse. Is there a correct local
   boundary closure for a real-time SU/NTU patch environment (the literature is mostly imaginary-time/ground
   state)?

---

## 6. Code / data pointers

- **Engine:** `scripts/g0/star_env.jl` — `build_star_bond_env` (4-direction env: `:bare` / `:exact_cluster` /
  `:ntu`), `project_star_env!` (lossless-star-then-env-truncate), `_build_star_bond_env_ntu` (NTU patch).
- **NTU builder + patches:** `scripts/g0/ntu_bondenv.jl` (validated horizontal builder; P6/P10/P12/torus).
- **CTMRG primitive:** `scripts/g0/ctm_bond_proto.jl` (bondenv_fu glue, value-correct), `ctm_chi_ladder.jl`
  (the χ-ladder probe).
- **Gates / validation:** `scripts/g0/envaware_integration.jl` (G0–G4, D-parameterized),
  `ntu_vertical_gate.jl` (vertical-P6 transpose oracle).
- **Trajectory drivers:** `scripts/g0/env_fork_trajectory.jl` (parameterized: tstart,tend,dt,env,alg,patch,D),
  `env_d3_control.jl`, `env_bare_control.jl`, `onset_intervention_probe.jl`, `native_d34_trajectory.jl`.
- **Data (all results / .out files):** `notes/stage2-truncation/data/`.
- **Full running log of the review process:** `notes/stage2-truncation/2026-06-07-workflow-fanout-review-verdict.md`.

## 7. Honest one-paragraph summary

A loop-carrying bond environment **provably stabilizes** the 2D PXP D=2 quench and reproduces the 4×4 revival
to ~2% of exact ED (the catastrophic bare-D2 blowup is an environment artifact, decisively). The residual ~2%
is a bond-dimension ceiling (dt-converged), so improving accuracy needs D≥3. **The obstruction is that no
environment we have is both accurate and feasible at D≥3 on a system small enough to have an exact oracle:**
the exact full-torus env is accurate but explodes (32 h/70 GB at D=3); the cheap local NTU patch is feasible
but structurally ~14% off and unrefinable on 4×4; CTMRG is viable and scalable but targets the thermodynamic
limit, so it cannot be validated against the finite 4×4 ED, and the cheap frozen-bond shortcuts are
gauge-confounded / non-discriminating. We are pausing here for external guidance on the validation strategy
(esp. how to validate a CTM-aware real-time evolution without an oracle, and whether the revival is a
finite-size feature) before committing to the (substantial) amortized CTM-aware trajectory build.
