# Prototype: validate the PEPSKit full-update bond-truncation pipeline on OUR
# native PXP state type, for ONE horizontal center-right bond. This is the core
# primitive of a CTM-aware star evolution; get the leg bookkeeping right here
# (by running) before integrating into project_star.
#
# Pipeline (mirrors PEPSKit's test/bondenv/*.jl idiom):
#   lossless gate on one star  ->  CTMRG env on the enlarged state  ->
#   _qr_bond(center,right) -> bondenv_fu -> positive_approx (PSD check) ->
#   bond_truncate(ALS, truncrank D) -> _qr_bond_undo -> measure density.
# Compare CTM-full-update truncation vs a bare SVD truncation of the SAME bond,
# both against the lossless (exact, untruncated) density.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

const D = 2
const BIG = 64
const CHI = 24
const T_FREEZE = 1.0

cell = PeriodicSquareUnitCell(4, 4)
psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = D)
evolve_pepskit!(psi, T_FREEZE; dt = 0.02, order = 2, schedule = :serial,
                maxdim = D, cutoff = 1e-12, rel_floor = 1e-3)
@printf("frozen density at t=%.2f: %.6f\n", T_FREEZE, exact_density_finite(psi; max_sites = 16))

# pick a center; its right leaf shares center.E <-> right.W (horizontal bond)
c0 = cell.reps[1]
cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
rc, cc = cidx[1], cidx[2]
Nr, Nc = size(psi.peps)
cp1 = mod1(cc + 1, Nc)
@printf("center PEPSKit (r,c)=(%d,%d), right leaf (r,c)=(%d,%d)\n", rc, cc, rc, cp1)

# --- lossless single-star gate (exact one-gate ceiling) ---
lossless = deepcopy(psi)
project_star_pepskit!(lossless, c0; dt = 0.02, maxdim = BIG,
                      cutoff = 0.0, rel_floor = 0.0, evolution = :real, projected = true)
n_lossless = exact_density_finite(lossless; max_sites = 16)
@printf("n_lossless (no truncation)              = %.8f\n", n_lossless)
A = lossless.peps.A[rc, cc]
B = lossless.peps.A[rc, cp1]
@printf("center-right bond dim after gate: %s\n", string(space(A)))

# --- CTMRG environment on the enlarged measurement PEPS ---
mp = measurement_peps(lossless)
env0 = PEPSKit.CTMRGEnv(randn, ComplexF64, mp, ComplexSpace(CHI))
env, info = PEPSKit.leading_boundary(env0, mp; alg = :simultaneous,
                                     tol = 1e-9, miniter = 2, maxiter = 100, verbosity = 1,
                                     trunc = PEPSKit.truncrank(CHI))
println("CTMRG done.")

# --- _qr_bond -> bondenv_fu -> positive_approx ---
X, a, b, Y = PEPSKit._qr_bond(A, B)
@printf("a space (from _qr_bond): %s\n", string(space(a)))
@printf("b space (from _qr_bond): %s\n", string(space(b)))
benv = PEPSKit.bondenv_fu(rc, cc, X, Y, env)
@printf("benv space: %s ; dual flags: %s\n", string(space(benv)),
        string([isdual(space(benv, ax)) for ax in 1:numind(benv)]))
Z = PEPSKit.positive_approx(benv)
@printf("positive_approx OK; cond(Z'Z)=%.3e\n", cond(Z' * Z))

# --- glue _qr_bond output to bond_truncate's expected {2,1}/{1,2}, leg2=phys ---
# bond_truncate wants a[env,phys; bond] (={2,1}), b[bond; phys,env] (={1,2}).
# _qr_bond gives a (q <- P,E) = {1,2} legs (q,P,E); b permuted to (?,?,?).
# Permute to the expected forms (verified by running):
a_bt = permute(a, ((1, 2), (3,)))   # (q,P ; E)  -> env=q, phys=P, bond=E
b_bt = permute(b, ((1,), (2, 3)))   # (E ; P,q?) -> guess; fix from error if wrong
@printf("a_bt space: %s\nb_bt space: %s\n", string(space(a_bt)), string(space(b_bt)))

alg = ALSTruncation(; trunc = PEPSKit.truncrank(D) & PEPSKit.truncerror(; atol = 1e-12),
                    maxiter = 50, tol = 1e-12)
a1, s1, b1, tinfo = PEPSKit.bond_truncate(a_bt, b_bt, benv, alg)
@printf("bond_truncate fid (CTM-aware) = %.8e ; kept dim = %d\n", tinfo.fid, dim(space(s1, 1)))

# --- reattach and measure density ---
a1f, b1f = PEPSKit.absorb_s(a1, s1, b1)
A2, B2 = PEPSKit._qr_bond_undo(X, a1f, b1f, Y)
ctm = deepcopy(lossless)
ctm.peps.A[rc, cc] = A2
ctm.peps.A[rc, cp1] = B2
n_ctm = exact_density_finite(ctm; max_sites = 16)

# --- bare SVD truncation of the SAME bond (no environment), for comparison ---
ab = PEPSKit._combine_ab(a_bt, b_bt)
a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc = PEPSKit.truncrank(D))
a0f, b0f = PEPSKit.absorb_s(a0, s0, b0)
# undo permute back to _qr_bond_undo's expected a/b forms (match a_bt/b_bt layout)
A0, B0 = PEPSKit._qr_bond_undo(X, permute(a0f, ((1,2),(3,))), permute(b0f, ((1,),(2,3))), Y)
bare = deepcopy(lossless)
bare.peps.A[rc, cc] = A0
bare.peps.A[rc, cp1] = B0
n_bare = exact_density_finite(bare; max_sites = 16)

println("\n===== RESULT (one horizontal bond) =====")
@printf("n_lossless         = %.8f\n", n_lossless)
@printf("n_CTM-aware trunc  = %.8f  (|dn| = %.4e)\n", n_ctm, abs(n_ctm - n_lossless))
@printf("n_bare-SVD trunc   = %.8f  (|dn| = %.4e)\n", n_bare, abs(n_bare - n_lossless))
println("DONE")
