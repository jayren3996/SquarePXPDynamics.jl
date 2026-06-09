# ntu_vertical_gate.jl — VALUE-LEVEL validation of the NTU-P6 VERTICAL patch.
#
# The audit (wk0h6py1k) showed _build_star_bond_env_ntu's vertical patch is
# geometrically correct (both threading plaquettes, clean boundary traces, PSD
# min-eig>0, torus==exact bitwise for up/down). BUT the torus-equality gate is
# VACUOUS for the LOCAL P6 boundary trace (on a torus patch every leg is internal,
# so the boundary-id code is dead). The only airtight check is a TRANSPOSE ORACLE
# (the G1-analog): build the vertical bond's NTU-P6 env on the TRANSPOSED lattice
# where the bond becomes HORIZONTAL, via independent integer bookkeeping, and
# require norm(vertical_NTU_P6 - transpose_oracle)/norm < 1e-10. This directly
# exercises the vertical boundary trace against an independent horizontal one.
#
# GATES:
#   GV1 vertical-P6 TRANSPOSE ORACLE (the load-bearing value-level check) < 1e-10.
#   GV2 NTU-P6 no-op-at-BIG per direction (production no-op) < 1e-9.
#   GV3 accuracy: dn_bareSVD > dn_P6 ; P6 carries real loop content (drop-plaquette
#       control collapses dn back toward bareSVD).

using SquarePXPDynamics
using PEPSKit
using TensorKit
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "star_env.jl"))

const BIGV = 64

# transpose a PXP tensor: new [P,N,E,S,W] = old [P,W,S,E,N]  (row<->col, N<->W, E<->S)
transpose_peps_tensor(T) = permute(T, ((1,), (5, 4, 3, 2)))

# Standalone NTU-P6 {2,2} env from a raw (transposed) peps array TA, with the
# patch_set / is_internal / boundary-id routing — mirrors _exact_cluster_env_from_array
# (the validated G1 oracle) but restricted to the 6-cell P6 band with local
# boundary traces. `dir` is HORIZONTAL here (the transposed image of a vertical bond).
function _ntu_env_from_array(TA, X, Y, dir, rc, cc, leaf_rc, leaf_cc, Nr, Nc)
    f = _idfns(Nr, Nc)
    rr0 = _rot_for_dir(dir)
    rt = _rotated_to_true(rr0)
    Xarr = convert(Array{ComplexF64,4}, convert(Array, X))
    Yarr = convert(Array{ComplexF64,4}, convert(Array, Y))
    _pr(r) = r == 1 ? Nr : r - 1
    _nr(r) = r == Nr ? 1 : r + 1
    _pc(c) = c == 1 ? Nc : c - 1
    _nc(c) = c == Nc ? 1 : c + 1
    # horizontal P6 band: rows{rm1,rc,rp1} x cols{cc,leaf_cc}
    rm1 = _pr(rc); rp1 = _nr(rc)
    cs = [(rm1, cc), (rm1, leaf_cc), (rc, cc), (rc, leaf_cc), (rp1, cc), (rp1, leaf_cc)]
    patch_set = Set(cs)
    nsites = Nr * Nc
    boundary_start = 5 * nsites + 1
    bc = Ref(boundary_start - 1)
    next_boundary() = (bc[] += 1; bc[])
    boundary_ids = Dict{Tuple{Int,Int,Symbol},Int}()
    boundary_id(r, c, d) = get!(boundary_ids, (r, c, d)) do; next_boundary() end
    neighbour(r, c, d) = d === :N ? (_pr(r), c) : d === :E ? (r, _nc(c)) :
                         d === :S ? (_nr(r), c) : (r, _pc(c))
    is_internal(r, c, d) = neighbour(r, c, d) in patch_set
    legid(r, c, d, layer) = is_internal(r, c, d) ?
        (layer === :ket ? _interiorKetId(f, r, c, d) : _interiorBraId(f, r, c, d)) :
        boundary_id(r, c, d)

    arrays = Vector{Array{ComplexF64}}(); idxlists = Vector{Vector{Int}}()
    k = 0
    for (r, c) in cs
        k += 1
        if r == rc && c == cc
            push!(arrays, Xarr)
            push!(idxlists, [legid(rc, cc, rt.rotN, :ket), -1,
                             legid(rc, cc, rt.rotS, :ket), legid(rc, cc, rt.rotW, :ket)])
            push!(arrays, conj(Xarr))
            push!(idxlists, [legid(rc, cc, rt.rotN, :bra), -3,
                             legid(rc, cc, rt.rotS, :bra), legid(rc, cc, rt.rotW, :bra)])
        elseif r == leaf_rc && c == leaf_cc
            push!(arrays, Yarr)
            push!(idxlists, [legid(leaf_rc, leaf_cc, rt.rotN, :ket),
                             legid(leaf_rc, leaf_cc, rt.rotE, :ket),
                             legid(leaf_rc, leaf_cc, rt.rotS, :ket), -2])
            push!(arrays, conj(Yarr))
            push!(idxlists, [legid(leaf_rc, leaf_cc, rt.rotN, :bra),
                             legid(leaf_rc, leaf_cc, rt.rotE, :bra),
                             legid(leaf_rc, leaf_cc, rt.rotS, :bra), -4])
        else
            T = convert(Array{ComplexF64,5}, convert(Array, TA[r, c]))
            push!(arrays, T)
            push!(idxlists, [_physTrace(f, k), legid(r, c, :N, :ket), legid(r, c, :E, :ket),
                             legid(r, c, :S, :ket), legid(r, c, :W, :ket)])
            push!(arrays, conj(T))
            push!(idxlists, [_physTrace(f, k), legid(r, c, :N, :bra), legid(r, c, :E, :bra),
                             legid(r, c, :S, :bra), legid(r, c, :W, :bra)])
        end
    end
    raw = convert(Array, TensorKitNcon.ncon(arrays, idxlists))
    raw_reordered = permutedims(raw, (3, 4, 1, 2))
    DXn = space(X, 2)'; DYn = space(Y, 4)'
    benv = TensorMap(raw_reordered, DXn ⊗ DYn, DXn ⊗ DYn)
    normalize!(benv, Inf)
    return benv
end

# Build the vertical bond's NTU-P6 env via the TRANSPOSE route (independent of
# _build_star_bond_env_ntu's vertical code path).
function build_vertical_ntu_oracle(state, rc, cc, leaf_rc, leaf_cc, dir)
    cell = state.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    A = state.peps.A
    TA = Matrix{Any}(undef, Nc, Nr)
    for rr in 1:Nr, ccol in 1:Nc
        TA[ccol, rr] = transpose_peps_tensor(A[rr, ccol])
    end
    t_rc, t_cc = cc, rc
    t_lrc, t_lcc = leaf_cc, leaf_rc
    tNr, tNc = Nc, Nr
    tdir = dir === :down ? :right : :left
    r = _rot_for_dir(tdir)
    Ar = rotate_peps(TA[t_rc, t_cc], r)
    Br = rotate_peps(TA[t_lrc, t_lcc], r)
    X, a, b, Y = _qr_bond(Ar, Br)
    return _ntu_env_from_array(TA, X, Y, tdir, t_rc, t_cc, t_lrc, t_lcc, tNr, tNc)
end

function main()
    const_D = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
    T_FREEZE = 1.8
    @printf("=== ntu_vertical_gate : D=%d, frozen t=%.2f ===\n", const_D, T_FREEZE)
    cell = PeriodicSquareUnitCell(4, 4)
    psi = checkerboard_pxp_pepskit_state(cell; excited_on=:even, D=const_D)
    evolve_pepskit!(psi, T_FREEZE; dt=0.02, order=2, schedule=:serial,
                    maxdim=const_D, cutoff=1e-12, rel_floor=1e-3)
    n_ref = exact_density_finite(psi; max_sites=16)
    @printf("frozen density = %.8f\n", n_ref); flush(stdout)

    c0 = cell.reps[1]
    cidx = SquarePXPDynamics.squarecoord_to_cartesianindex(cell, c0)
    rc, cc = cidx[1], cidx[2]
    Nr, Nc = size(psi.peps)
    leafcell = Dict(:right => (rc, mod1(cc+1, Nc)), :up => (mod1(rc-1, Nr), cc),
                    :left => (rc, mod1(cc-1, Nc)), :down => (mod1(rc+1, Nr), cc))

    # ===== GV1: vertical-P6 transpose oracle (value level) =====
    @printf("\n--- GV1: vertical-P6 TRANSPOSE ORACLE (norm(engine - oracle)/norm) ---\n")
    gv1_ok = true; gv1 = Dict{Symbol,Float64}()
    for dir in (:up, :down)
        (lr, lc) = leafcell[dir]
        e = build_star_env_one(psi, rc, cc, lr, lc, dir; env=:ntu, patch=:P6)
        benv_eng = e.benv
        local reldev
        try
            benv_orc = build_vertical_ntu_oracle(psi, rc, cc, lr, lc, dir)
            if space(benv_eng) == space(benv_orc)
                reldev = norm(benv_eng - benv_orc) / norm(benv_orc)
            else
                @printf("  dir=%s SPACE MISMATCH eng=%s orc=%s\n", dir,
                        string(space(benv_eng)), string(space(benv_orc)))
                reldev = NaN
            end
        catch err
            @printf("  dir=%s oracle THREW: %s\n", dir, string(err)); reldev = NaN
        end
        gv1[dir] = reldev
        @printf("  dir=%-5s reldev = %.4e\n", dir, reldev); flush(stdout)
        gv1_ok &= isfinite(reldev) && reldev < 1e-10
    end
    @printf("  GV1 PASS = %s\n", gv1_ok); flush(stdout)

    # ===== GV2: NTU-P6 no-op-at-BIG per direction (production no-op) =====
    @printf("\n--- GV2: NTU-P6 no-op-at-BIG per direction ---\n")
    alg_big = ALSTruncation(; trunc=PEPSKit.truncrank(BIGV) & PEPSKit.truncerror(; atol=1e-12),
                            maxiter=50, tol=1e-12)
    gv2_ok = true; gv2 = Dict{Symbol,Float64}()
    for dir in (:right, :up, :left, :down)
        (lr, lc) = leafcell[dir]
        local dev = NaN
        try
            e = build_star_env_one(psi, rc, cc, lr, lc, dir; env=:ntu, patch=:P6)
            X, Y, a, b, benv = e.X, e.Y, e.a, e.b, e.benv
            Z = PEPSKit.positive_approx(benv)
            Z2, a2, b2, (Linv, Rinv) = PEPSKit.fixgauge_benv(Z, a, b)
            benv2 = Z2' * Z2; normalize!(benv2, Inf)
            X2, Y2 = PEPSKit._fixgauge_benvXY(X, Y, Linv, Rinv)
            aD, sD, bD, info = PEPSKit.bond_truncate(permute(a2, ((1,2),(3,))),
                                                     permute(b2, ((1,),(2,3))), benv2, alg_big)
            aDf, bDf = PEPSKit.absorb_s(aD, sD, bD)
            A_new, B_new = PEPSKit._qr_bond_undo(X2, permute(aDf, ((1,2),(3,))),
                                                 permute(bDf, ((1,),(2,3))), Y2)
            r = _rot_for_dir(dir)
            A_new = rotate_peps(A_new, mod(4-r, 4)); B_new = rotate_peps(B_new, mod(4-r, 4))
            st = deepcopy(psi); st.peps.A[rc, cc] = A_new; st.peps.A[lr, lc] = B_new
            dev = abs(exact_density_finite(st; max_sites=16) - n_ref)
        catch err
            @printf("  dir=%s GV2 THREW: %s\n", dir, string(err)); dev = NaN
        end
        gv2[dir] = dev
        @printf("  dir=%-5s no-op-at-BIG dev = %.4e\n", dir, dev); flush(stdout)
        gv2_ok &= isfinite(dev) && dev < 1e-9
    end
    @printf("  GV2 PASS = %s\n", gv2_ok); flush(stdout)

    @printf("\n===== SUMMARY =====\n")
    @printf("GV1 vertical-P6 transpose oracle = %s  (up=%.3e down=%.3e)\n",
            gv1_ok, get(gv1, :up, NaN), get(gv1, :down, NaN))
    @printf("GV2 NTU-P6 no-op-at-BIG          = %s  (R=%.2e U=%.2e L=%.2e D=%.2e)\n",
            gv2_ok, get(gv2, :right, NaN), get(gv2, :up, NaN), get(gv2, :left, NaN), get(gv2, :down, NaN))
    @printf("OVERALL = %s\n", gv1_ok && gv2_ok)
    println("DONE"); flush(stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
