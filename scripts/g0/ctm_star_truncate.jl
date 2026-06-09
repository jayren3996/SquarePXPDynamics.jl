# Full per-star CTM-aware truncation of all 4 center-leaf bonds, with no-op
# self-checks. Architecture: lossless gate -> CTMRG on enlarged state ->
# truncate right/left (horizontal) directly, up/down via rotl90 (vertical->horizontal).
#
# Self-check: at maxdim=BIG every bond truncation is a no-op, so the density must
# reproduce n_lossless EXACTLY. This validates the leg/rotation geometry before any
# real (maxdim=D) reduction is trusted. Then compare CTM-aware vs bare SVD at maxdim=D.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

const T_FREEZE = 1.0
const Dtest = 2
const BIG = 64
const CHI = 24

# Truncate the horizontal bond (r,c)-(r,c+1) of `peps` (InfinitePEPS) using CTM `env`.
# Mutates peps.A[r,c] and peps.A[r,cp1]. trunc = D (Int) for truncrank, or BIG.
function truncate_bond_h!(peps, env, r, c, D; use_ctm = true)
    Nr, Nc = size(peps)
    cp1 = mod1(c + 1, Nc)
    A, B = peps.A[r, c], peps.A[r, cp1]
    X, a, b, Y = PEPSKit._qr_bond(A, B)
    a_bt = permute(a, ((1, 2), (3,)))   # (env,phys ; bond)
    b_bt = permute(b, ((1,), (2, 3)))   # (bond ; phys,env)
    alg = ALSTruncation(; trunc = PEPSKit.truncrank(D) & PEPSKit.truncerror(; atol = 1e-12),
                        maxiter = 50, tol = 1e-12)
    if use_ctm
        benv = PEPSKit.bondenv_fu(r, c, X, Y, env)
        a1, s1, b1, _ = PEPSKit.bond_truncate(a_bt, b_bt, benv, alg)
    else
        ab = PEPSKit._combine_ab(a_bt, b_bt)
        a1, s1, b1 = PEPSKit.svd_trunc(permute(ab, ((1, 3), (4, 2))); trunc = PEPSKit.truncrank(D))
    end
    a1f, b1f = PEPSKit.absorb_s(a1, s1, b1)
    A2, B2 = PEPSKit._qr_bond_undo(X, permute(a1f, ((1, 2), (3,))), permute(b1f, ((1,), (2, 3))), Y)
    peps.A[r, c] = A2
    peps.A[r, cp1] = B2
    return nothing
end

function ctmrg_env(peps)
    env0 = PEPSKit.CTMRGEnv(randn, ComplexF64, peps, ComplexSpace(CHI))
    env, _ = PEPSKit.leading_boundary(env0, peps; alg = :simultaneous, tol = 1e-9,
                                      miniter = 2, maxiter = 100, verbosity = 0,
                                      trunc = PEPSKit.truncrank(CHI))
    return env
end

# Build a frozen lossless-gated state + its CTM env. Returns (state, env, n_lossless, rc, cc).
function setup()
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = Dtest)
    evolve_pepskit!(psi, T_FREEZE; dt = 0.02, order = 2, schedule = :serial,
                    maxdim = Dtest, cutoff = 1e-12, rel_floor = 1e-3)
    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]
    project_star_pepskit!(psi, c0; dt = 0.02, maxdim = BIG, cutoff = 0.0, rel_floor = 0.0,
                          evolution = :real, projected = true)
    n_lossless = exact_density_finite(psi; max_sites = 16)
    env = ctmrg_env(measurement_peps(psi))
    return psi, env, n_lossless, rc, cc
end

# Density after replacing state.peps with `peps` (InfinitePEPS).
function density_with(state, peps)
    s2 = deepcopy(state)
    for i in eachindex(peps.A)
        s2.peps.A[i] = peps.A[i]
    end
    return exact_density_finite(s2; max_sites = 16)
end

println("=== NO-OP geometry self-checks (maxdim=BIG; each must == n_lossless) ===")
psi, env, n_lossless, rc, cc = setup()
Nr, Nc = size(measurement_peps(psi))
@printf("n_lossless = %.8f ; center (r,c)=(%d,%d) ; grid %dx%d\n", n_lossless, rc, cc, Nr, Nc)

# right
let peps = deepcopy(measurement_peps(psi))
    truncate_bond_h!(peps, env, rc, cc, BIG; use_ctm = true)
    @printf("right (H) no-op: %.8f  dev=%.2e\n", density_with(psi, peps), abs(density_with(psi, peps) - n_lossless))
end
# left
let peps = deepcopy(measurement_peps(psi))
    truncate_bond_h!(peps, env, rc, mod1(cc - 1, Nc), BIG; use_ctm = true)
    @printf("left  (H) no-op: %.8f  dev=%.2e\n", density_with(psi, peps), abs(density_with(psi, peps) - n_lossless))
end
# up/down via rotl90: the original up/down bonds are the HORIZONTAL bonds carrying
# the enlarged dim (>Dtest) in the rotated frame. Find them rather than guess coords.
function enlarged_h_bonds(peps; thresh = Dtest)
    Nr2, Nc2 = size(peps)
    bonds = Tuple{Int, Int}[]
    for r in 1:Nr2, c in 1:Nc2
        dim(space(peps.A[r, c], 3)) > thresh && push!(bonds, (r, c))  # E (leg 3) bond enlarged
    end
    return bonds
end
let peps0 = deepcopy(measurement_peps(psi))
    pr = rotl90(peps0); er = rotl90(env)   # rotate the SAME fixed-point env (su_iter idiom)
    hb = enlarged_h_bonds(pr)
    @printf("  (rotated enlarged horizontal bonds: %s)\n", string(hb))
    # DIAGNOSTIC: dump spaces at each enlarged bond + its E-neighbor + env edges
    Nr2, Nc2 = size(pr)
    for (r, c) in hb
        cp = mod1(c + 1, Nc2)
        @printf("  bond (%d,%d)-(%d,%d): A=%s  B=%s\n", r, c, r, cp,
                string(space(pr.A[r, c])), string(space(pr.A[r, cp])))
        X, a, b, Y = PEPSKit._qr_bond(pr.A[r, c], pr.A[r, cp])
        @printf("    X=%s  Y=%s\n", string(space(X)), string(space(Y)))
        @printf("    env edges t4[%d,%d]=%s t1X[%d,%d]=%s\n",
                r, mod1(c-1,Nc2), string(space(er.edges[4, r, mod1(c-1,Nc2)])),
                mod1(r-1,Nr2), c, string(space(er.edges[1, mod1(r-1,Nr2), c])))
    end
    for (r, c) in hb
        truncate_bond_h!(pr, er, r, c, BIG; use_ctm = true)
    end
    peps_back = rotr90(pr)
    nb = density_with(psi, peps_back)
    @printf("up/down (rot) no-op: %.8f  dev=%.2e\n", nb, abs(nb - n_lossless))
end

println("\n=== ALL FOUR bonds, maxdim=D=$Dtest : CTM-aware vs bare SVD ===")
function truncate_all(use_ctm)
    peps = deepcopy(measurement_peps(psi))
    truncate_bond_h!(peps, env, rc, cc, Dtest; use_ctm)                 # right
    truncate_bond_h!(peps, env, rc, mod1(cc - 1, Nc), Dtest; use_ctm)   # left
    pr = rotl90(peps); er = use_ctm ? rotl90(env) : env                 # rotate the SAME fixed-point env
    for (r, c) in enlarged_h_bonds(pr)                                  # up/down (now horizontal)
        truncate_bond_h!(pr, er, r, c, Dtest; use_ctm)
    end
    return density_with(psi, rotr90(pr))
end
n_ctm = truncate_all(true)
n_bare = truncate_all(false)
@printf("n_lossless        = %.8f\n", n_lossless)
@printf("n_CTM-aware (D=%d) = %.8f  |dn|=%.4e\n", Dtest, n_ctm, abs(n_ctm - n_lossless))
@printf("n_bare-SVD  (D=%d) = %.8f  |dn|=%.4e\n", Dtest, n_bare, abs(n_bare - n_lossless))
println("DONE")
