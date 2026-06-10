# Exact constrained-sector finite-torus density for a PEPSKit-native PXP state
# (PXPIPEPSState), include-able after `using SquarePXPDynamics`.
#
# Idea (meet-in-the-middle constrained-sector contraction): the PXP state lives
# (up to tiny truncation leakage) in the hard-square-constrained sector, and for
# a FIXED physical configuration the network is single-layer with bond dim D.
#   1. Enumerate constrained row configs (no two adjacent excited on a ring of
#      Nc columns; bit c-1 set = column c EXCITED = P index 1).
#   2. For each PEPS row r (peps[r,:], legs [P,N,E,S,W]) and config m build
#      Row_r[m]: a (prod N-dims x prod S-dims) matrix; P fixed by m, E-W ring
#      contracted (wrap leg traced), N legs fused into the row index and S legs
#      into the column index WITH THE SAME column ordering (col 1 fastest), so
#      Row_r[m] * Row_{r+1}[m'] sums S(r,c) against N(r+1,c) per column.
#   3. Rows m, m' vertically stackable iff (m & m') == 0.
#   4./5. Halves: B_top[p] = Row_1[c1]*...*Row_h[c_h] (h = Nr/2), maps
#      [w = N legs of row 1 (the Nr<->1 seam)] x [v = S legs of row h (middle
#      seam)]; B_bot[s] = Row_{h+1}*...*Row_Nr maps [v] x [w]. Both built
#      tree-shared by DFS (each depth-k product is one GEMM, reused by all its
#      depth-(k+1) children).
#   6. Amplitudes for ALL (p,s) in one BLAS GEMM via tr(A*B) = vec(A^T).vec(B):
#      Psi = transpose(Tmat) * Bmat, Tmat[:,p] = vec(transpose(B_top[p])),
#      Bmat[:,s] = vec(B_bot[s]).
#   7. Pair (p,s) valid iff compat(c_h, c_{h+1}) AND compat(c_Nr, c1); invalid
#      pairs still get amplitudes — they measure constraint leakage at the two
#      seams.
#   8. n = sum_valid wgt*nex / (Nr*Nc*sum_valid wgt); leak = sum_invalid wgt /
#      sum_all wgt, with wgt = abs2(Psi).
#
# Per-row matrices are rescaled by a single per-row factor (same factor for all
# configs of that row), which multiplies EVERY amplitude by the same constant
# exp(logscale): n and leak are invariant; physical amplitudes are
# Psi * exp(logscale) (used by the 4x4 gate's per-amplitude cross-check).
#
# Cost at 6x6 D=3: row cache 6*18 x 8.5 MB ~ 1 GB; halves 2 x 2309 x 531441
# complex ~ 39 GB; final GEMM 2309 x 531441 x 2309 ~ minutes at 16 BLAS threads.

using SquarePXPDynamics
using LinearAlgebra

"""
    sector_row_masks(n)

All row patterns (UInt64, n bits) with no two adjacent excited sites on a
periodic ring of n; bit value 1 means "excited". Logic copied verbatim from
`constrained_row_masks` in scripts/run_neel_ed_krylov.jl.
"""
function sector_row_masks(n::Int)
    masks = UInt64[]
    last = UInt64(1) << (n - 1)
    for m::UInt64 = 0:(UInt64(1) << n - UInt64(1))
        has_linear = (m & (m << 1)) & ((UInt64(1) << n) - UInt64(1)) != 0
        has_wrap = (m & UInt64(1)) != 0 && (m & last) != 0
        (has_linear || has_wrap) || push!(masks, m)
    end
    return masks
end

# Site tensors of measurement_peps(psi), sliced by physical value and permuted
# to leg order (W, N, S, E). sites[r, c, p] for p in 1:2 (1 = EXCITED).
function _sector_site_slices(peps, Nr::Int, Nc::Int)
    sites = Array{Array{ComplexF64,4},3}(undef, Nr, Nc, 2)
    for r = 1:Nr, c = 1:Nc
        A = convert(Array{ComplexF64,5}, convert(Array, peps[r, c]))  # [P,N,E,S,W]
        size(A, 1) == 2 || error("physical dim must be 2 at ($r,$c), got $(size(A,1))")
        for p = 1:2
            sites[r, c, p] = permutedims(A[p, :, :, :, :], (4, 1, 3, 2))  # (W,N,S,E)
        end
    end
    # bond-dim consistency on the torus: E(r,c) ~ W(r,c+1), S(r,c) ~ N(r+1,c)
    for r = 1:Nr, c = 1:Nc
        A = sites[r, c, 1]
        cn = c == Nc ? 1 : c + 1
        rn = r == Nr ? 1 : r + 1
        size(A, 4) == size(sites[r, cn, 1], 1) || error("E/W bond dim mismatch at ($r,$c)")
        size(A, 3) == size(sites[rn, c, 1], 2) || error("S/N bond dim mismatch at ($r,$c)")
    end
    return sites
end

# Row transfer matrix of PEPS row r for config mask m: maps fused N legs (row
# index, column 1 fastest) -> fused S legs (col index, column 1 fastest).
# Ring contraction carries the wrap leg w = W of site (r,1) and traces it
# against E of site (r,Nc) at the end. Interleaved leg order is maintained
# automatically: T flat order is [w, N1,S1, N2,S2, ..., e].
function _sector_row_matrix(sites, r::Int, m::UInt64, Nc::Int)
    pof(c) = ((m >> (c - 1)) & 0x1) == 0x1 ? 1 : 2   # bit set -> P index 1 (excited)
    A1 = sites[r, 1, pof(1)]
    wdim = size(A1, 1)
    ndims_ = Int[size(A1, 2)]
    sdims_ = Int[size(A1, 3)]
    T = reshape(A1, :, size(A1, 4))                   # [w*N1*S1, e1]
    for c = 2:Nc
        Ac = sites[r, c, pof(c)]
        size(T, 2) == size(Ac, 1) || error("E/W mismatch in row $r at col $c")
        T = T * reshape(Ac, size(Ac, 1), :)           # append (N_c, S_c, e_c)
        push!(ndims_, size(Ac, 2))
        push!(sdims_, size(Ac, 3))
        T = reshape(T, :, size(Ac, 4))                # split trailing e_c
    end
    edim = size(T, 2)
    edim == wdim || error("wrap E/W mismatch in row $r")
    mid = prod(ndims_) * prod(sdims_)
    T3 = reshape(T, wdim, mid, edim)
    Rint = zeros(ComplexF64, mid)                     # trace wrap leg
    for d = 1:wdim
        @views Rint .+= T3[d, :, d]
    end
    dims = Vector{Int}(undef, 2Nc)                    # (n1,s1,n2,s2,...)
    for c = 1:Nc
        dims[2c-1] = ndims_[c]
        dims[2c] = sdims_[c]
    end
    Rt = reshape(Rint, dims...)
    perm = vcat(1:2:2Nc, 2:2:2Nc)                     # -> (n1..nNc, s1..sNc)
    return reshape(permutedims(Rt, perm), prod(ndims_), prod(sdims_))
end

# Cache all Nr x nmask row matrices; rescale each row family by one shared
# factor (max Frobenius norm) for over/underflow safety. Returns (rowm,
# logscale) with every amplitude scaled by exp(-logscale) relative to physical.
function _sector_row_cache(sites, masks::Vector{UInt64}, Nr::Int, Nc::Int)
    nm = length(masks)
    rowm = [Vector{Matrix{ComplexF64}}(undef, nm) for _ = 1:Nr]
    logscale = 0.0
    for r = 1:Nr
        for (i, m) in enumerate(masks)
            rowm[r][i] = _sector_row_matrix(sites, r, m, Nc)
        end
        s = maximum(norm(M) for M in rowm[r])
        (isfinite(s) && s > 0) || error("row $r: zero/non-finite transfer matrices")
        for M in rowm[r]
            M ./= s
        end
        logscale += log(s)
    end
    return rowm, logscale
end

# counts[k] = number of vertically compatible chains of length k (k = 1..h)
function _sector_chain_counts(compat::Matrix{Bool}, h::Int)
    nm = size(compat, 1)
    ways = ones(Int, nm)
    counts = Vector{Int}(undef, h)
    counts[1] = sum(ways)
    for k = 2:h
        ways = [sum(ways[i] for i = 1:nm if compat[i, j]; init = 0) for j = 1:nm]
        counts[k] = sum(ways)
    end
    return counts
end

# Build one half: DFS over compatible mask chains of length h = length(rows),
# tree-shared products (one GEMM per node at depth >= 2, reused by children).
# Writes column p of OUT as vec(transpose(prod)) (transpose_out, for the top
# half) or vec(prod) (bottom half). Returns per-stack (first mask idx, last
# mask idx, nex, chain) in enumeration order.
function _sector_build_half!(
    OUT::Matrix{ComplexF64},
    rows::Vector{Vector{Matrix{ComplexF64}}},
    masks::Vector{UInt64},
    compat::Matrix{Bool};
    transpose_out::Bool,
)
    h = length(rows)
    nm = length(masks)
    nstk = size(OUT, 2)
    firsts = Vector{Int}(undef, nstk)
    lasts = Vector{Int}(undef, nstk)
    nexs = Vector{Int}(undef, nstk)
    stacks = Vector{Vector{Int}}(undef, nstk)
    rdim = size(rows[1][1], 1)
    bufs = [Matrix{ComplexF64}(undef, rdim, size(rows[k][1], 2)) for k = 1:h]
    tbuf = Matrix{ComplexF64}(undef, size(rows[h][1], 2), rdim)
    path = Vector{Int}(undef, h)
    p = Ref(0)
    function rec(k::Int, prodmat::Matrix{ComplexF64}, first_mi::Int, nex::Int)
        for mi = 1:nm
            (k == 1 || compat[path[k-1], mi]) || continue
            path[k] = mi
            nexk = nex + count_ones(masks[mi])
            local cur::Matrix{ComplexF64}
            if k == 1
                cur = rows[1][mi]
            else
                mul!(bufs[k], prodmat, rows[k][mi])
                cur = bufs[k]
            end
            fmi = k == 1 ? mi : first_mi
            if k == h
                p[] += 1
                col = p[]
                if transpose_out
                    permutedims!(tbuf, cur, (2, 1))
                    @views OUT[:, col] .= vec(tbuf)
                else
                    @views OUT[:, col] .= vec(cur)
                end
                firsts[col] = fmi
                lasts[col] = mi
                nexs[col] = nexk
                stacks[col] = copy(path)
            else
                rec(k + 1, cur, fmi, nexk)
            end
        end
        return nothing
    end
    rec(1, bufs[1], 0, 0)
    p[] == nstk || error("stack count mismatch: built $(p[]), expected $nstk")
    return firsts, lasts, nexs, stacks
end

"""
    sector_density(psi::PXPIPEPSState; kwargs...) -> NamedTuple

Exact constrained-sector finite-torus density of `measurement_peps(psi)` via
meet-in-the-middle contraction over hard-square-valid configurations. Requires
even Nr. Returns `(n, leak, secs, ...)` where `n` is the mean excitation
density over the valid (fully constrained) sector, `leak` the weight fraction
of seam-violating configurations among all enumerated pairs, and `secs` the
wall time of this call. Extra diagnostic fields: ntop, nbot, K, npairs,
nvalid_pairs, est_gb, logscale, S_valid, S_all, masks, top_/bot_ stack data,
and (if `keep_amplitudes`) the full amplitude matrix `amps` (scaled by
exp(-logscale) relative to physical amplitudes).

Guards: throws before allocating if the estimated memory exceeds
`ram_guard_gb` or the estimated time exceeds `time_guard_s`.
"""
function sector_density(
    psi::PXPIPEPSState;
    blas_threads::Int = 16,
    ram_guard_gb::Float64 = 150.0,
    time_guard_s::Float64 = 1800.0,
    keep_amplitudes::Bool = false,
    verbose::Bool = true,
)
    t0 = time()
    LinearAlgebra.BLAS.set_num_threads(blas_threads)
    cell = psi.unitcell
    Nr, Nc = cell.Ly, cell.Lx
    iseven(Nr) || throw(ArgumentError("sector_density requires even Nr; got $Nr"))
    (Nr >= 2 && Nc >= 2) || throw(ArgumentError("cell too small"))
    h = Nr ÷ 2

    sites = _sector_site_slices(measurement_peps(psi), Nr, Nc)
    masks = sector_row_masks(Nc)
    nm = length(masks)
    compat = [(masks[i] & masks[j]) == 0 for i = 1:nm, j = 1:nm]

    rowm, logscale = _sector_row_cache(sites, masks, Nr, Nc)
    Mw = size(rowm[1][1], 1)        # fused N of row 1  (Nr<->1 seam)
    Mv = size(rowm[h][1], 2)        # fused S of row h  (middle seam)
    size(rowm[h+1][1], 1) == Mv || error("middle seam dim mismatch")
    size(rowm[Nr][1], 2) == Mw || error("wrap seam dim mismatch")
    K = Mw * Mv

    counts = _sector_chain_counts(compat, h)
    ntop = counts[h]
    nbot = ntop                      # same masks/compat for the bottom half
    ngemm_half = h >= 2 ? sum(counts[2:h]) : 0
    est_gb =
        (Float64(K) * (ntop + nbot) + Float64(ntop) * nbot + Float64(Nr) * nm * K) *
        16 / 2^30
    est_flops =
        2 * ngemm_half * 8.0 * Mw * Float64(Mv)^2 + 8.0 * ntop * Float64(K) * nbot
    est_secs = est_flops / 5.0e10    # very conservative 50 GFLOP/s
    if verbose
        println(
            "  [sector] Nr=$Nr Nc=$Nc nmask=$nm h=$h ntop=$ntop K=$K " *
            "npairs=$(ntop*nbot) est_gb=$(round(est_gb; digits=2)) " *
            "est_secs(cons.)=$(round(est_secs; digits=1))",
        )
        flush(stdout)
    end
    est_gb <= ram_guard_gb ||
        error("sector_density RAM guard: est $(round(est_gb;digits=1)) GB > $(ram_guard_gb) GB")
    est_secs <= time_guard_s ||
        error("sector_density time guard: est $(round(est_secs;digits=0)) s > $(time_guard_s) s")

    Tmat = Matrix{ComplexF64}(undef, K, ntop)
    tf, tl, tnex, tstacks =
        _sector_build_half!(Tmat, rowm[1:h], masks, compat; transpose_out = true)
    Bmat = Matrix{ComplexF64}(undef, K, nbot)
    bf, bl, bnex, bstacks =
        _sector_build_half!(Bmat, rowm[h+1:Nr], masks, compat; transpose_out = false)
    if verbose
        println("  [sector] halves built at t+$(round(time()-t0; digits=1)) s; GEMM...")
        flush(stdout)
    end
    Psi = transpose(Tmat) * Bmat     # Psi[p,s] = tr(B_top[p] * B_bot[s])

    S_all = 0.0
    S_valid = 0.0
    S_valid_n = 0.0
    nvalid = 0
    @inbounds for s = 1:nbot
        bfs = bf[s]
        bls = bl[s]
        bn = bnex[s]
        for p = 1:ntop
            w = abs2(Psi[p, s])
            S_all += w
            if compat[tl[p], bfs] && compat[bls, tf[p]]   # mid seam AND wrap seam
                S_valid += w
                S_valid_n += w * (tnex[p] + bn)
                nvalid += 1
            end
        end
    end
    S_valid > 0 || error("zero weight in valid sector")
    n = S_valid_n / S_valid / (Nr * Nc)
    leak = (S_all - S_valid) / S_all
    secs = time() - t0
    if verbose
        println(
            "  [sector] n=$(n) leak=$(leak) nvalid_pairs=$nvalid secs=$(round(secs; digits=1))",
        )
        flush(stdout)
    end
    return (
        n = n,
        leak = leak,
        secs = secs,
        ntop = ntop,
        nbot = nbot,
        K = K,
        npairs = ntop * nbot,
        nvalid_pairs = nvalid,
        est_gb = est_gb,
        logscale = logscale,
        S_valid = S_valid,
        S_all = S_all,
        masks = masks,
        top_first = tf,
        top_last = tl,
        top_nex = tnex,
        bot_first = bf,
        bot_last = bl,
        bot_nex = bnex,
        top_stacks = tstacks,
        bot_stacks = bstacks,
        amps = keep_amplitudes ? Psi : nothing,
    )
end

# In-script verification of the trace-flattening identity used for the big
# GEMM (step 6): tr(A*B) == transpose(vec(transpose(A))) * vec(B), in exactly
# the production layout Psi = transpose(Tmat) * Bmat. Non-square on purpose.
let
    A = ComplexF64[1+2im 3-1im 0.5im; -2.0 1im 4.0]            # 2x3
    B = ComplexF64[0.3 1im; 2-1im 0.7; 1.0 1+1im]              # 3x2
    direct = tr(A * B)
    Tmat = reshape(vec(permutedims(A)), :, 1)                  # vec(transpose(A))
    Bmat = reshape(vec(B), :, 1)
    viaflat = (transpose(Tmat)*Bmat)[1, 1]
    isapprox(direct, viaflat; rtol = 1e-13) ||
        error("trace-flatten identity check failed: $direct vs $viaflat")
end
