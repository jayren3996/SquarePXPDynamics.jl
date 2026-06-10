# Probe: is the G3 full-star no-op-at-BIG gap (4.97e-9) an ALS-convergence floor
# or a routing/composition error? Test by (a) tightening ALS maxiter/tol, and
# (b) substituting FullEnvTruncation (exact, non-iterative) for the BIG no-op.

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "star_env.jl"))

# env_truncate_bond! variant taking an explicit alg (so we can swap ALS params /
# FullEnvTruncation). Mirrors star_env.jl:env_truncate_bond! exactly.
function env_truncate_bond_alg!(state, rc, cc, leaf_rc, leaf_cc, dir, alg; env=:exact_cluster)
    A = state.peps.A[rc, cc]; B = state.peps.A[leaf_rc, leaf_cc]
    r = _rot_for_dir(dir)
    Ar = rotate_peps(A, r); Br = rotate_peps(B, r)
    X, a, b, Y = _qr_bond(Ar, Br)
    benv = build_star_bond_env(state, X, Y, dir, rc, cc, leaf_rc, leaf_cc; env=env)
    Z = positive_approx(benv)
    Z2, a2, b2, (Linv, Rinv) = fixgauge_benv(Z, a, b)
    benv2 = Z2' * Z2; normalize!(benv2, Inf)
    X2, Y2 = _fixgauge_benvXY(X, Y, Linv, Rinv)
    a2_bt = permute(a2, ((1, 2), (3,))); b2_bt = permute(b2, ((1,), (2, 3)))
    aD, sD, bD, info = bond_truncate(a2_bt, b2_bt, benv2, alg)
    aDf, bDf = absorb_s(aD, sD, bD)
    A_new, B_new = _qr_bond_undo(X2, permute(aDf, ((1, 2), (3,))),
                                 permute(bDf, ((1,), (2, 3))), Y2)
    A_new = rotate_peps(A_new, mod(4 - r, 4)); B_new = rotate_peps(B_new, mod(4 - r, 4))
    state.peps.A[rc, cc] = A_new; state.peps.A[leaf_rc, leaf_cc] = B_new
    return info.fid
end

function full_star_noop(psi, c0, rc, cc, Nr, Nc, alg; label="")
    leafcell = Dict(:right=>(rc,mod1(cc+1,Nc)), :up=>(mod1(rc-1,Nr),cc),
                    :left=>(rc,mod1(cc-1,Nc)), :down=>(mod1(rc+1,Nr),cc))
    lossless = deepcopy(psi)
    project_star_pepskit!(lossless, c0; dt=0.02, maxdim=64, cutoff=0.0,
                          rel_floor=0.0, evolution=:real, projected=true)
    n_lossless = exact_density_finite(lossless; max_sites=16)
    st = deepcopy(psi)
    project_star_pepskit!(st, c0; dt=0.02, maxdim=10^6, cutoff=0.0,
                          rel_floor=0.0, evolution=:real, projected=true)
    fids = Float64[]
    for dir in (:right,:up,:left,:down)
        (lr,lc) = leafcell[dir]
        push!(fids, env_truncate_bond_alg!(st, rc, cc, lr, lc, dir, alg))
    end
    n_env = exact_density_finite(st; max_sites=16)
    @printf("  [%s] dev = %.4e ; fids = %s\n", label, abs(n_env-n_lossless),
            string(round.(1 .- fids, sigdigits=3))); flush(stdout)
    return abs(n_env - n_lossless)
end

function main()
    cell = PeriodicSquareUnitCell(4,4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on=:even, D=2)
    evolve_pepskit!(psi, 1.8; dt=0.02, order=2, schedule=:serial, maxdim=2, cutoff=1e-12, rel_floor=1e-3)
    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]; Nr, Nc = size(psi.peps)
    @printf("=== G3 probe: full-star no-op-at-BIG, ALS-floor vs routing ===\n"); flush(stdout)

    als50 = ALSTruncation(; trunc=truncrank(64)&truncerror(;atol=1e-12), maxiter=50, tol=1e-12)
    als300 = ALSTruncation(; trunc=truncrank(64)&truncerror(;atol=1e-14), maxiter=300, tol=1e-14)
    full_star_noop(psi, c0, rc, cc, Nr, Nc, als50;  label="ALS maxiter=50  tol=1e-12")
    full_star_noop(psi, c0, rc, cc, Nr, Nc, als300; label="ALS maxiter=300 tol=1e-14")
    fet = FullEnvTruncation(; trunc=truncrank(64)&truncerror(;atol=1e-14))
    try
        full_star_noop(psi, c0, rc, cc, Nr, Nc, fet; label="FullEnvTruncation (exact)")
    catch err
        @printf("  [FullEnvTruncation] THREW: %s\n", string(err)); flush(stdout)
    end
    println("DONE")
end

main()
