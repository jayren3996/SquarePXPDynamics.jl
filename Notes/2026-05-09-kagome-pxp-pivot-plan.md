# Kagome PXP Pivot Plan

Date: 2026-05-09. Author: agent-assisted planning session.

## Why this note exists

We started the project on the **triangular lattice** with the goal of running ScarFinder for the 2D triangular PXP model with iPEPS evolution. After implementing the dense 7-site star gate, projected PXP machinery, D=1 Simple Update oracle, and an initial ScarFinder loop, we hit two structural blockers when extending Simple Update to D > 1:

1. **Algorithmic mismatch on triangular.** The 7-site star gate spans 1 center + 6 spokes. Sequential SVD peel-split of the post-gate cluster does not preserve translational invariance: for a 3-site unit cell, three "sibling" spoke tensors at the same rep are NOT three versions of the same site tensor; they each update a *different* bond direction of that rep. Naive average-and-error-out is mathematically incoherent. Repair requires a per-bond merge solving a small consistency system, with no clean published recipe for triangular star gates.

2. **Cost wall.** Cluster materialization scales as `2^7 × D^18` worst case. Realistic peak intermediates (with smart contraction order) hit GBs at D=4 and TBs at D≥6. Even with HOSVD or BP-gauge methods the *application* of the multi-body gate stays expensive on triangular.

Decision: **switch the lattice to Kagome**. Kagome is non-bipartite (frustrated via corner-sharing triangles), has mature PEPS/PESS infrastructure, and has dramatically smaller multi-body gates (5-site instead of 7-site, `2^5 × D^k` cluster scaling).

This note records the literature scan that motivated the switch and lays out a concrete algorithmic plan for the Kagome PXP ScarFinder pipeline.

---

## Literature scan

### A. Kagome PXP / Rydberg scars

- **Samajdar, Ho, Pichler, Lukin, Sachdev**, "Quantum phases of Rydberg atoms on a kagome lattice", PNAS 2021. arXiv: [2011.12295](https://arxiv.org/abs/2011.12295). DMRG ground-state phase diagram. Identifies √7 phase (12-fold degenerate at appropriate blockade range), maps to triangular dimer model. Studies ground states, not dynamics or scars. **Useful baseline for which Kagome PXP phases exist.**

- **Sahay, Mohapatra, Vishwanath et al.**, "Systematic construction of quantum many-body scars in frustrated Rydberg arrays", arXiv: [2605.05297](https://arxiv.org/abs/2605.05297) (May 2026, days old). Graph-theoretic framework for QMBSs on arbitrary lattices including frustrated geometries. Two mechanisms: type-I (locally entangled states for mild frustration) and type-II (strong frustration pinning). Hexagonal lattice example shown. **Provides the conceptual hook for "kagome PXP scars are an open question worth pursuing."**

- **Kerschbaumer et al.**, "Quantum Many-Body Scars beyond the PXP model in Rydberg simulators", arXiv: [2410.18913](https://arxiv.org/abs/2410.18913) (PRL 134, 160401, 2025). 1D chains and triangular ladder (quasi-1D). Uses MPS/TEBD with bond dim χ=2-3. **Not directly applicable but confirms scars exist beyond standard PXP.**

- **No published "kagome PXP scars" paper exists.** This is a genuine research opening — both the physics question and the numerical algorithmics are open. 

- **Surace et al.**, "Quantum simulation of dynamical gauge theories in periodically driven Rydberg atom arrays", arXiv: [2408.02733](https://arxiv.org/abs/2408.02733) (2024). Uses Kagome PXP framework for U(1) lattice gauge simulation via Floquet engineering. Confirms kagome PXP is a meaningful target Hamiltonian.

### B. PESS / simplex-tensor methods for Kagome

- **Xie, Chen, Yang, Yao, Xiang**, "Tensor renormalization of quantum many-body systems using projected entangled simplex states", PRX 4, 011025 (2014). arXiv: [1307.5696](https://arxiv.org/abs/1307.5696). The original PESS paper. **Designed for Kagome.** Each up-triangle gets a rank-3 simplex tensor S (no physical leg); each kagome site has a 3-leg tensor T (1 phys + 2 simplex legs). Simple Update via HOSVD on simplex tensors. State-of-the-art kagome Heisenberg ground-state energies. Already in our literature folder as `pess_simple_update_1307.5696.pdf`.

- **Liao, Wang, Xie, et al.** various follow-ups — PESS is the canonical kagome ansatz; many ground-state benchmarks exist. PEPS without the simplex structure converges poorly on kagome because corner-sharing-triangle frustration is hard to capture with bond singular values alone.

- **Time-evolution work in PESS framework is mostly imaginary-time** (ground state). Real-time PESS evolution exists but is less mature; the framework supports it.

### C. BP-gauge SU and BP-based dynamics

- **Alkabetz & Arad**, "Tensor networks contraction and the belief propagation algorithm", PRR 3, 023073 (2021). arXiv: [2008.04433](https://arxiv.org/abs/2008.04433). Establishes that BP on a tensor network is **equivalent to the mean-field environment used in vanilla Simple Update**. Foundational reference for "SU is BP."

- **Tindall & Fishman**, "Gauging tensor networks with belief propagation", SciPost Phys. 15, 222 (2023). arXiv: [2306.17837](https://arxiv.org/abs/2306.17837). BP-based gauging algorithm. Demonstrated on square, hexagonal, random regular, cubic. **Not explicitly tested on kagome or triangular**, but the algorithm is lattice-agnostic. Confirms regauging via BP + simple update gate application is mathematically equivalent to using BP messages as environments.

- **Tindall, Mello, Fishman, Stoudenmire, Sels**, "Dynamics of disordered quantum systems with two- and three-dimensional tensor networks", arXiv: [2503.05693](https://arxiv.org/abs/2503.05693) (March 2025). **Real-time evolution with BP-SU at scale.** Lattices: cylindrical 2D, diamond cubic 3D, dimerized cubic 3D. Bond dims χ_BP = 6 to 32 during evolution. **Two-body gates only.** State-of-the-art accuracy on disordered Ising-type Hamiltonians. Not tested on kagome or with multi-body gates.

- **Tindall, Fishman, Stoudenmire, Sels**, "Efficient tensor network simulation of IBM's Eagle kicked Ising experiment", PRX Quantum 5, 010308 (2024). arXiv: [2306.14887](https://arxiv.org/abs/2306.14887). BP on heavy-hex lattice (tree-friendly geometry). High-accuracy simulation outperforming the quantum hardware.

- **Anonymous(?) authors**, "Belief Propagation and Tensor Network Expansions for Many-Body Quantum Systems: Rigorous Results and Fundamental Limits", arXiv: [2604.03228](https://arxiv.org/abs/2604.03228) (recent). **Rigorous bounds.** Loop-decay condition guarantees BP exponential accuracy in gapped phases; systematic failure near criticality. Tested on TFIM 2D/3D.

- **Beyond-BP / cluster-corrected** BP: arXiv: [2510.02290](https://arxiv.org/abs/2510.02290), [2604.24760](https://arxiv.org/abs/2604.24760), [2411.04957](https://arxiv.org/abs/2411.04957). Loop corrections for non-tree lattices, with the goal of improving accuracy on geometries with short loops (which Kagome has via corner-sharing triangles). Active development frontier.

### D. Key takeaways from the scan

1. **Kagome PXP scar dynamics is unstudied** in published work. There's a legitimate research project here.
2. **PESS is THE canonical Kagome ansatz**, and it's well-validated for ground states.
3. **BP-SU is well-validated for 2-body bond gates** on tree-like and moderately-loopy lattices, including 2D/3D dynamics. Not validated for multi-body gates or for kagome specifically.
4. **The multi-body PXP gate problem doesn't disappear by changing methods.** Whether we use cluster-and-split SU, BP-SU with cluster gate, or PESS-SU, applying the 5-site PXP star gate requires a 5-site cluster operation. The improvement from triangular to kagome is in cluster size (5 vs 7 sites, 8 vs 18 external bonds), not in the algorithmic approach to multi-body gates.
5. **Loop corrections may matter for kagome.** Corner-sharing triangles produce 3-cycles that vanilla BP doesn't see. Either accept the loop error for first-cut work or plan for cluster-corrected BP later.

---

## Kagome geometry primer (for the design that follows)

Kagome lattice: 2D non-bipartite lattice of corner-sharing triangles. Each site has **coordination 4**. The natural unit cell contains **3 sites** (one up-triangle, labeled A, B, C). The medial lattice of kagome is the triangular lattice; the dual is honeycomb-of-triangles.

Standard Bravais vectors: `a1 = (2, 0)`, `a2 = (1, √3)` (in some natural unit). Sublattice positions: `A = (0, 0)`, `B = (1, 0)`, `C = (1/2, √3/2)` within each unit cell.

Each site belongs to **2 triangles** (one up-triangle, one down-triangle). Each triangle contains 3 sites of distinct sublattices.

The 4 NNs of a center site `c` (of sublattice X) are:
- 2 NNs in the up-triangle containing `c` (sublattices Y and Z, say).
- 2 NNs in the down-triangle containing `c` (also sublattices Y and Z, but distinct lattice sites).

So in the natural 3-site UC, the 5-site star {center + 4 NNs} contains positions on sublattices {X, Y, Z, Y', Z'} which wrap to reps {X, Y, Z, Y, Z} under the 3-site UC. **Same sibling-aliasing issue we hit on triangular**, just smaller (2 spokes per non-center rep instead of 3).

To eliminate aliasing: use a **9-site enlarged UC** (3×3 block of the natural up-triangle). With 9 reps, the 5 star positions land on 5 distinct reps. This is the cleanest path for a cluster-and-split SU. Cost: 3× more reps to track but each rep tensor is the same size, so memory is 3× and per-step work is ~3× compared to the natural 3-site UC.

---

## Proposed algorithm: PESS-SU + BP-gauge cluster gate, 9-site UC

The plan combines three established ingredients:

### Component 1: PESS ansatz (Xie 2014)

State representation:
- One **simplex tensor `S_△`** per up-triangle: rank-3, no physical leg, 3 virtual legs (one per kagome site of the triangle).
- One **simplex tensor `S_▽`** per down-triangle: same shape.
- One **site tensor `T_v`** per kagome site: 1 phys leg + 2 simplex legs (one toward the up-triangle the site belongs to, one toward the down-triangle).
- Per-bond λ spectra on the 6 bonds of each triangle (3 bonds per triangle × 2 triangles per site = 6 unique bonds per site, exactly matching standard PESS bookkeeping).

Why PESS over plain iPEPS: PESS natively captures the 3-body correlations within each triangle via the simplex tensor `S`. On kagome this is essential because the relevant frustration unit IS the triangle. Plain iPEPS with a 4-leg site tensor wastes representation capacity trying to encode triangle correlations through bond-only entanglement.

### Component 2: BP-gauge maintenance

After each star-gate application we re-establish the BP fixed-point gauge on the PESS tensors. Operationally:

1. Run BP message-passing on the PESS tensor network until convergence (or fixed iteration cap).
2. Use BP messages as the local environment for any subsequent local truncation step.
3. The PESS λ spectra act as the "Vidal form" representation of the BP fixed point on each bond.

Tindall-Fishman 2023 ([2306.17837](https://arxiv.org/abs/2306.17837)) shows this is equivalent to using BP-message environments for local truncation. On kagome, the corner-sharing-triangle structure gives 3-cycles (the triangles themselves) which introduce loop corrections. For first-cut work we accept the BP-tree approximation and validate empirically; loop corrections (cluster-BP) are a follow-up.

### Component 3: 5-site cluster gate application

The PXP star gate `G_eff = P_blockade · exp(-i dt · X_c P_1 P_2 P_3 P_4)` is intrinsically 5-body. Algorithm per gate application:

1. **Identify the 5 star positions** on the lattice. With the 9-site enlarged UC these are 5 distinct reps; no aliasing.
2. **Absorb λ on all bonds of each star position.** Each kagome site has 4 bonds (2 to adjacent triangles' simplex tensors, 2 to within-triangle partners — wait, actually in PESS each site has 2 simplex-tensor legs, NOT 4 bond legs to other sites).

Let me restate the PESS bond structure precisely:
- Each kagome site `v` has tensor `T_v` with 3 legs: phys + simplex-up-leg + simplex-down-leg.
- Each up-triangle has tensor `S_△` with 3 legs (one per site of the triangle).
- Each down-triangle has tensor `S_▽` with 3 legs.
- λ spectra live on the bonds between site tensors and simplex tensors (6 such bonds per kagome site... no wait, 2 per site, since each site connects to 2 simplices).

So site `v` has 2 bond legs (one to the up-triangle simplex, one to the down-triangle simplex), each with its own λ spectrum.

For the 5-site star centered at `c`:
- `c` connects to 2 simplices: up-triangle `△_c` and down-triangle `▽_c`.
- Each NN of `c` is in either `△_c` or `▽_c`. Specifically: 2 NNs in `△_c` (the other two sites of `△_c`) and 2 NNs in `▽_c`.
- So the 5-site star, as a PESS network, consists of: `T_c, T_{NN1_△}, T_{NN2_△}, T_{NN1_▽}, T_{NN2_▽}, S_△, S_▽` = **5 site tensors + 2 simplex tensors = 7 tensors** in the cluster.

3. **Contract the gate.** Apply `G_eff` (a 5-phys-leg ITensor with 5 in-phys + 5 out-phys) into the cluster's 5 phys legs. The resulting cluster has 5 fresh out-phys legs + external bond legs (those of NN sites that go to other simplices outside the star = 2 per NN × 4 NNs = 8 external simplex-bond legs).

4. **Decompose back to PESS tensors.** This is the multi-step part:
   - **Refit each site tensor** by SVD across the cut "this site's phys + this site's outgoing simplex leg | rest". Standard simple-update factorization.
   - **Refit each simplex tensor** by HOSVD across the 3 legs of the simplex.
   - Use BP messages as the environment for each truncation step (to make the truncation locally-optimal in the BP-environment sense).
   - Update the affected λ spectra; renormalize.

5. **Re-run BP** to refresh the gauge before the next gate application.

### Cluster scaling estimate

With 5 site tensors + 2 simplex tensors in the cluster:
- Site tensors: phys=2, simplex bonds dim=D. Each `T_v` has size `2 × D^2`.
- Simplex tensors: each `S_△` has size `D^3`.
- Worst-case fully-contracted cluster: `2^5 × D^8` external × `D^?` internal. Actually since simplex tensors have rank 3 and absorb the within-triangle correlations, the effective "external bond dimension" is dominated by the 8 dangling simplex legs (2 NNs × 2 simplices × 2 per NN = 8 net external).
- D=4: `32 × 65,536 ≈ 2M entries × 16 bytes = 32 MB`. Comfortable.
- D=8: `32 × 16M ≈ 0.5G entries × 16 = 8 GB`. Tight but workstation-feasible.
- D=12: `32 × ~430M ≈ 14G entries × 16 = 220 GB`. Not feasible on a workstation; doable on a server.
- D=16: `32 × 4.3G ≈ 137G entries × 16 = 2.2 TB`. Not feasible.

**Practical operating range: D = 4 to D = 8.** Above D=8 we'd need either (a) cluster-and-split with smarter intermediate contraction order, (b) cluster-BP for loop corrections, or (c) NTU/full-update infrastructure. ScarFinder for kagome PXP at D=4-8 should already produce meaningful results; expansion to D>8 is a research milestone not a starting requirement.

---

## Validation strategy

Same three-layer structure as the (now-defunct) triangular plan, adapted for kagome:

### Layer 1: Kernel ED tests on the 5-site cluster

- Construct a hand-built PESS state with known D=2 entries.
- Compute the 5-site cluster vector by direct contraction (32-dim vector).
- Apply a Haar-random unitary `G` directly: `psi_after_ED = G · psi_before`.
- Apply via `apply_star_gate_simple_update!` on the PESS state.
- Reconstruct the cluster vector from the post-gate PESS state.
- Assert agreement up to global phase, tested at:
  - large `maxdim` (no truncation): exact match within `1e-10`.
  - small `maxdim` (force truncation): bond dims `<= maxdim`, `discarded_weight > 0`, finite norms.

### Layer 2: Finite-torus integration test

- 12-site or 24-site kagome torus matched to the 9-site UC.
- Hilbert space ≤ `2^24 = 16M`, ED-tractable on a workstation.
- Initial state: all-down product (blockade-allowed). Run 3-5 projected PXP Trotter steps on both PESS-SU iPEPS and torus ED. Compare per-sublattice `<Z>`.
- Tolerance: `1e-3` absolute on local Z observables at `dt=0.01`, second-order Trotter, `D=4`.

### Layer 3: Solvable kagome stabilizer benchmark

- Kagome cluster-state Hamiltonian `H_cluster = -∑_v X_v ∏_{u ∈ NN(v)} Z_u`. All terms commute → no Trotter error. Initial state: all-Z+ product.
- Closed-form `<Z_c>(t)` evolution (analogous to `cluster_center_z_expectation_exact` in `src/SolvableModels.jl`).
- Run real-time evolution at `D=4`, measure drift from analytic = pure truncation error. Target: `<1e-6` for short times.

---

## What carries over from the triangular work

Reusable directly:

- **Spin operators** ([src/SpinOps.jl](src/SpinOps.jl)) — basis-independent.
- **Schedules** ([src/Schedules.jl](src/Schedules.jl)) — first/second-order Trotter weights apply to any color partition.
- **SolvableModels** structure ([src/SolvableModels.jl](src/SolvableModels.jl)) — analytic helpers, model template stays.
- **Models** ([src/Models.jl](src/Models.jl)) — `pxp_star_hamiltonian` parametrized by gate dimension + projector + flip; will need a 32×32 kagome variant alongside the 128×128 triangular one.
- **Gates** ([src/Gates.jl](src/Gates.jl)) — `dense_gate` / `projected_gate` are agnostic to the gate dimension; small refactor to remove the hard-coded `(128, 128)` size check.
- **ScarFinder loop structure** ([src/ScarFinder.jl](src/ScarFinder.jl)) — the search/seed/iterate/rank machinery is geometry-agnostic.
- **ITensors.jl backend, root Project.toml, test harness** — unchanged.

Replaced:

- **Geometry** ([src/Geometry.jl](src/Geometry.jl)) — triangular axial coords → kagome (Bravais coords + sublattice index).
- **States** ([src/States.jl](src/States.jl)) — `TriangularIPEPS` → `KagomePESS` with simplex tensors.
- **Observables** ([src/Observables.jl](src/Observables.jl)) — local-expectation logic adapted for PESS structure.
- **SimpleUpdate** ([src/SimpleUpdate.jl](src/SimpleUpdate.jl)) — full rewrite around PESS HOSVD updates + BP-gauge maintenance. The 7 dead helpers from the triangular rank-1 placeholder are deleted; the 4 new helpers (absorb / build cluster / peel-or-HOSVD-split / writeback) get reimagined for the PESS structure.
- **Evolution** ([src/Evolution.jl](src/Evolution.jl)) — drives the kagome scheduler over the kagome coloring (smaller, since kagome has a cleaner 3-coloring structure than triangular's 7-coloring).

New code:

- **BP message-passing module** — `src/BP.jl`. Iterative message updates on the PESS tensor network. Convergence detection. Gauge restoration.
- **Finite-cluster ED helpers** — `test/util_kagome_finite_ed.jl`. Torus Hamiltonian builder, ED time-evolution reference, local observables on the torus.

---

## Phased plan (sketch — formal task plan goes in `docs/superpowers/plans/` after re-brainstorming)

1. **Geometry rewrite.** Kagome Bravais lattice, 3-site natural UC and 9-site enlarged UC, sublattice partitioning, color schedule for parallel star updates. Tests for neighbor relations and color disjointness.

2. **PESS state container.** `KagomePESS` with site + simplex tensors, λ bookkeeping per simplex bond, Vidal-form invariants. Product-state and random initializers. Tests for shape, opposite-bond sharing, and norm.

3. **5-site projected PXP gate.** 32×32 dense gate, blockade projector for the 5-site cluster (5 internal star edges to check), tests against hand-computed reference matrix elements.

4. **BP module.** Message-passing, convergence, BP-fixed-point gauge restoration. Tests on simple product-state BP fixed points.

5. **PESS Simple Update for site-product gates.** Carries over from triangular path with cluster-size adaptations. The product-gate fast path validates the absorb/split/writeback skeleton without needing the multi-body machinery yet.

6. **PESS Simple Update for general 5-site gates with BP-gauge truncation.** This is the meat. Cluster-and-split with BP-environment-aware truncation. Resolves the sibling-aliasing problem by construction (9-site UC) and uses BP messages instead of vanilla λ for environment.

7. **Validation layers 1-3.**

8. **ScarFinder driver kagome adaptation.** Mostly re-wiring; the core search loop is unchanged.

9. **End-to-end smoke test.** Kagome PXP ScarFinder at D=4, ThreeSiteUnitCell, modest niterations. Verify bond growth and truncation activate, candidate diagnostics finite, deterministic ranking.

---

## Open questions to decide before writing the formal spec

1. **PESS or plain iPEPS for the first kagome implementation?** PESS is the right end state but plain iPEPS is closer to the existing triangular code. Trade-off: PESS rewrite is bigger upfront; plain iPEPS has worse algorithmic fit and may need to be rewritten anyway. **Lean: PESS from the start.**

2. **9-site UC always, or 3-site UC with sibling-merge?** 9-site avoids the symmetry problem entirely at 3× cost. 3-site is half-baked unless we figure out the sibling merge — and we already know that's hard. **Lean: 9-site.**

3. **BP-gauge from day one, or vanilla SU first?** Vanilla SU is simpler to validate; BP adds a layer. But BP is the right truncation environment, and rewriting SU twice is wasteful. **Lean: vanilla λ environment for first cut, BP integration as a follow-on PR after validation layers 1-3 are passing.**

4. **Real-time evolution validated against torus ED only, or also against published kagome dynamics work?** No published kagome PXP dynamics exists, so torus ED is the only oracle. **Confirmed: torus ED is the validation target.**

5. **Bond dimension target.** First-cut: D=4. Production: D=8. Aspirational: D=12 (would need server-class memory). **Lean: D=4 in CI tests, D=8 in benchmarks, D=12 documented as out-of-scope for the first PR.**

---

## Risks and explicit non-goals for the kagome pivot

Risks:

- **No precedent for kagome PXP dynamics in the literature.** Validation depends entirely on internal ED references (torus + stabilizer benchmark). If physics intuition disagrees with results, no published numbers to cross-check.
- **PESS rewrite cost.** Bigger than the triangular SU we drafted; 1-3 weeks of focused work.
- **BP loop corrections.** Corner-sharing triangles → 3-cycles. Vanilla BP gauges them as trees; loop error is real but bounded for gapped phases. May need cluster-BP follow-up if scarring physics turns out to be at-or-near a critical regime.
- **D=8+ memory.** Cluster scaling improves a lot from triangular but still hits a wall around D=12-16. Not a blocker for first scientific results; is a blocker for benchmark-comparable accuracy.

Explicit non-goals (kept as future work):

- Loop-corrected BP / cluster-BP. Vanilla BP first, cluster corrections only if validation says we need them.
- NTU / full update on kagome. PESS-SU + BP-gauge first; NTU is the next-tier accuracy upgrade.
- Energy expectation at D>1 via boundary MPS or CTMRG. PESS local energy estimators are approximate; precise expectation values can wait.
- Imaginary-time energy correction in ScarFinder. Depends on (#3) being ready.
- 4-site, larger UCs beyond the 9-site enlargement.
- Triangular PXP. Triangular work is parked; if kagome scarfinder produces interesting results we may revisit triangular with a PESS-on-triangular ansatz later.

---

## Next concrete step

Re-run the brainstorming → spec → plan flow on this kagome pivot. The brainstorming will scope the PESS/BP/UC choices above; the spec will pin the kagome PESS data layout and SU update kernel; the plan will produce task-by-task implementation steps analogous to the triangular plan but adapted.

Tentative target file paths:

- `docs/superpowers/specs/2026-05-09-kagome-pess-pxp-dynamics-design.md`
- `docs/superpowers/plans/2026-05-09-kagome-pess-pxp-dynamics.md`

The triangular plan and spec stay in the repo as historical artifacts; they are not active.
