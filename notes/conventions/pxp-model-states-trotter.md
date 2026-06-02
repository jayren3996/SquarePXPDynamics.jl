# PXP model, states, and the checkerboard Trotter split

See also `definitions_and_conventions.md`, `physics_context.md` in this folder.

## Model

Square-lattice PXP star term: `h_c = X_c · P_down(right) · P_down(up) ·
P_down(left) · P_down(down)` — Pauli-X on the center site times projectors onto
|down⟩ on its 4 neighbors. `H = Σ_c h_c`. Basis order `1 = :up`, `2 = :down`;
star site order `(center, right, up, left, down)`. The dense square-star
Hamiltonian (`PXPModel.jl`) is the source of truth for the energy operator.

## Initial-state constraint (verified, with a legible error since 2026-06-02)

- `:z_up` (all-up) ALWAYS crashes the projected star update — every center has
  excited neighbors, so the projected gate is blockade-forbidden.
- `:checkerboard_*` / Néel WORKS on even×even cells with `schedule = :serial`;
  it only crashes on ODD×ODD cells (the checkerboard has a seam).
- `:down` (all-down) is the safe quench state on any cell.
- Practical rule: use Néel on even×even cells (`schedule = :serial`), or `:down`.
  The crash now throws a legible error naming the fix.

## Checkerboard = an EXACT 2-term Trotter split

`[h_c, h_c'] = 0` unless centers c, c' are nearest neighbors. Two terms collide
ONLY through the center X (`X_c` hits `h_c'` iff c ∈ {c'} ∪ N(c')). When only
their P-supports overlap (diagonal / dist-2 / knight), they STILL commute —
verified to machine zero (`scripts/verify_pxp_commutation.jl`: non-NN ‖[h,h']‖=0,
NN = 4√2). So the checkerboard 2-coloring A={(x+y) even}, B={(x+y) odd} gives
`exp(-i H_A t) = ∏_{c∈A} exp(-i h_c t)` EXACTLY (verified to 1e-15 on 4×2, 6×2
tori, dt up to 1.3). PXP is a clean 2-term split H = H_A + H_B, like 1D TEBD.

**Use:** a `:checkerboard` schedule + order-4 (Yoshida S4 on the exact 2-term
split) is cleaner than 5-color for higher-order Trotter and cuts Trotter error at
fixed dt. CAVEAT: this only reduces TROTTER error — at the revival the error is the
truncation/rise lag, not Trotter, so order-4 helps short time, not the revival.
