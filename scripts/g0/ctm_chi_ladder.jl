# ctm_chi_ladder.jl — de-risk the CTMRG pivot (cheap, frozen-bond).
#
# On the SAME frozen D=2 rise state used by reviewer_ntu_probe (lossless star ->
# one horizontal bond), build the CTMRG env at chi in {8,16,24,32} and compare its
# {2,2} bond env (bondenv_fu) to the EXACT-CLUSTER finite-torus reference
# (build_exact_cluster_bondenv). Two questions in one run:
#   (1) FINITE-vs-INFINITE gap: CTMRG gives the infinite-tiling env of the 4x4 cell;
#       exact-cluster is the finite 4x4 torus (matches the ED oracle). How far apart?
#   (2) chi-CONVERGENCE: does the CTM env settle as chi grows (its own internal limit)?
# Metrics (the rel-dev is the SAME one that gave NTU-P6 = 0.949, directly comparable):
#   reldev_inf = ||benv_ctm - benv_exact|| / ||benv_exact||   (both normalize!(Inf))
#   1 - overlap = 1 - |<benv_ctm,benv_exact>_F| / (||.||_F ||.||_F)   (scale-invariant)
#   |dn|        = single-bond density deviation vs lossless (weak discriminator, for record)
# vs NTU-P6 reldev 0.949 ("nearly orthogonal"): if CTM reldev drops well below that
# (-> small as chi grows), CTMRG's env is a genuine accuracy lever (pivot GO).

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "ntu_bondenv.jl"))   # -> build_exact_cluster_bondenv (+ spike include)

const D = 2
const BIG = 64
const T_FREEZE = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.0   # 1.0 = HEALTHY (pre-blowup)
const CHIS = [8, 16, 24, 32, 48]

frob_overlap(A, B) = real(dot(A, B)) / (norm(A) * norm(B))   # in [-1,1]

# truncate the bond with a given {2,2} benv, return |n_trunc - n_lossless|
function bond_dn(lossless, benv, X, a, b, Y, rc, cc, cp1, n_lossless)
    alg = ALSTruncation(; trunc=PEPSKit.truncrank(D) & PEPSKit.truncerror(; atol=1e-12),
                        maxiter=50, tol=1e-12)
    Z = PEPSKit.positive_approx(benv)
    Z2, a2, b2, (Linv, Rinv) = PEPSKit.fixgauge_benv(Z, a, b)
    benv2 = Z2' * Z2; normalize!(benv2, Inf)
    X2, Y2 = PEPSKit._fixgauge_benvXY(X, Y, Linv, Rinv)
    a1, s1, b1, tinfo = PEPSKit.bond_truncate(permute(a2, ((1,2),(3,))),
                                              permute(b2, ((1,),(2,3))), benv2, alg)
    a1f, b1f = PEPSKit.absorb_s(a1, s1, b1)
    A2, B2 = PEPSKit._qr_bond_undo(X2, permute(a1f, ((1,2),(3,))), permute(b1f, ((1,),(2,3))), Y2)
    st = deepcopy(lossless); st.peps.A[rc, cc] = A2; st.peps.A[rc, cp1] = B2
    return abs(exact_density_finite(st; max_sites=16) - n_lossless), tinfo.fid
end

function main()
    @printf("=== ctm_chi_ladder : D=%d frozen t=%.2f, chi=%s ===\n", D, T_FREEZE, string(CHIS))
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on=:even, D=D)
    evolve_pepskit!(psi, T_FREEZE; dt=0.02, order=2, schedule=:serial,
                    maxdim=D, cutoff=1e-12, rel_floor=1e-3)
    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]; Nr, Nc = size(psi.peps); cp1 = mod1(cc+1, Nc)

    lossless = deepcopy(psi)
    project_star_pepskit!(lossless, c0; dt=0.02, maxdim=BIG, cutoff=0.0, rel_floor=0.0,
                          evolution=:real, projected=true)
    n_lossless = exact_density_finite(lossless; max_sites=16)
    @printf("frozen n(t=%.2f)=%.6f ; n_lossless(1 star)=%.8f ; bond=%s\n",
            T_FREEZE, exact_density_finite(psi; max_sites=16), n_lossless,
            string(space(lossless.peps.A[rc, cc]))); flush(stdout)

    A = lossless.peps.A[rc, cc]; B = lossless.peps.A[rc, cp1]
    X, a, b, Y = PEPSKit._qr_bond(A, B)

    # exact-cluster finite-torus reference env + its single-bond |dn|
    benv_exact = build_exact_cluster_bondenv(lossless, X, Y, rc, cc, cp1)
    normalize!(benv_exact, Inf)
    dn_exact, fid_exact = bond_dn(lossless, benv_exact, X, a, b, Y, rc, cc, cp1, n_lossless)
    @printf("exact-cluster ref: |dn|=%.4e fid=%.8f ; benv space=%s\n",
            dn_exact, fid_exact, string(space(benv_exact))); flush(stdout)

    # bare-SVD single-bond |dn| (no env), for the ladder floor
    ab = PEPSKit._combine_ab(permute(a, ((1,2),(3,))), permute(b, ((1,),(2,3))))
    a0, s0, b0 = PEPSKit.svd_trunc(permute(ab, ((1,3),(4,2))); trunc=PEPSKit.truncrank(D))
    a0f, b0f = PEPSKit.absorb_s(a0, s0, b0)
    A0, B0 = PEPSKit._qr_bond_undo(X, permute(a0f, ((1,2),(3,))), permute(b0f, ((1,),(2,3))), Y)
    stb = deepcopy(lossless); stb.peps.A[rc, cc] = A0; stb.peps.A[rc, cp1] = B0
    dn_bare = abs(exact_density_finite(stb; max_sites=16) - n_lossless)
    @printf("bare-SVD (no env): |dn|=%.4e\n", dn_bare); flush(stdout)

    mp = measurement_peps(lossless)
    @printf("\n%-5s %-9s %-13s %-13s %-11s %-11s %-9s\n",
            "chi", "ctm_conv", "reldev_inf", "1-overlap", "|dn|_ctm", "fid_ctm", "wall_s")
    prev = nothing
    for chi in CHIS
        t0 = time()
        local reldev_inf = NaN; local one_minus_ov = NaN; local dn_ctm = NaN; local fid_ctm = NaN
        local conv = "?"
        try
            env0 = PEPSKit.CTMRGEnv(randn, ComplexF64, mp, ComplexSpace(chi))
            env, ctminfo = PEPSKit.leading_boundary(env0, mp; alg=:simultaneous,
                                tol=1e-9, miniter=4, maxiter=200, verbosity=2,
                                trunc=PEPSKit.truncrank(chi))
            conv = try
                ci = ctminfo
                hasproperty(ci, :converged) ? string(ci.converged) :
                    (ci isa Number ? string(ci) : "done")
            catch; "done" end
            benv_ctm = PEPSKit.bondenv_fu(rc, cc, X, Y, env)
            normalize!(benv_ctm, Inf)
            if space(benv_ctm) == space(benv_exact)
                reldev_inf = norm(benv_ctm - benv_exact) / norm(benv_exact)
                one_minus_ov = 1 - abs(frob_overlap(benv_ctm, benv_exact))
            else
                @printf("  chi=%d SPACE MISMATCH ctm=%s\n", chi, string(space(benv_ctm)))
            end
            dn_ctm, fid_ctm = bond_dn(lossless, benv_ctm, X, a, b, Y, rc, cc, cp1, n_lossless)
        catch err
            @printf("  chi=%d THREW: %s\n", chi, string(err))
        end
        @printf("%-5d %-9s %-13.4e %-13.4e %-11.4e %-11.6f %-9.1f\n",
                chi, conv, reldev_inf, one_minus_ov, dn_ctm, fid_ctm, time()-t0); flush(stdout)
        prev = chi
    end

    @printf("\n===== READ =====\n")
    @printf("NTU-P6 reldev_inf vs exact = 0.949 (nearly orthogonal). If CTM reldev_inf drops\n")
    @printf("well below that and SHRINKS with chi -> CTM env is an accuracy lever -> pivot GO.\n")
    @printf("If reldev_inf PLATEAUS at a nonzero floor -> that floor is the finite(4x4)-vs-\n")
    @printf("infinite(CTM) BC gap (bears on the revival finite-size question). exact-cluster\n")
    @printf("ref |dn|=%.4e ; bare |dn|=%.4e.\n", dn_exact, dn_bare)
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
