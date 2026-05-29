#!/usr/bin/env julia
#
# Single-Neel PXP ED dynamics on n×n PBC for n that's too large for brute-force
# enumeration (2^(n^2) too big). Builds the constrained basis by row-by-row
# composition: enumerate constrained rows, then compose into n-row stacks with
# vertical PXP constraint, then close PBC at the top↔bottom wraparound.
#
# Initial state: literal Neel-A (even sublattice excited; bit 0 in our
# convention means |up>/excited). Time evolution via Lanczos/Arnoldi Krylov on
# the sparse H = sum_i X_i restricted to constrained subspace. Observable:
# global density n(t) (fraction of excited sites).
#
# Saves a JSON with time grid + density trajectory.
#
# Environment variables:
#   NEEL_KRYLOV_N           lattice size (default 6)
#   NEEL_KRYLOV_TOTAL_TIME  (default 3.0)
#   NEEL_KRYLOV_DT_MEAS     (default 0.1)
#   NEEL_KRYLOV_DT_STEP     internal Krylov step (default 0.05)
#   NEEL_KRYLOV_M           Krylov subspace size (default 30)
#   NEEL_KRYLOV_OUTPUT_JSON (default artifacts/neel_ed_krylov_N.json)

using Pkg
project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using LinearAlgebra
using SparseArrays
using JSON3
using Printf

say(args...) = (println(args...); flush(stdout))

# ----------- Constrained basis via row-composition -----------

site_index(n::Int, x::Int, y::Int) = x + n * y + 1   # 1-based

"""
    constrained_row_masks(n)

Return all row patterns (UInt64 of width n bits) with no two adjacent excited
sites under periodic boundary conditions in 1D. Bit 0 == "up" (excited).
We use LSB == column 0, so the row is excited at column x iff bit x is 0
(after we encode "excited == 0" → mask bit value 1 means "excited").
For internal row enumeration we use the OPPOSITE convention: mask bit 1 means
"excited" — easier to work with. Translation back is done where needed.
"""
function constrained_row_masks(n::Int)
    masks = UInt64[]
    last = UInt64(1) << (n - 1)
    for m::UInt64 = 0:(UInt64(1) << n - UInt64(1))
        has_linear = (m & (m << 1)) & ((UInt64(1) << n) - UInt64(1)) != 0
        has_wrap = (m & UInt64(1)) != 0 && (m & last) != 0
        (has_linear || has_wrap) || push!(masks, m)
    end
    return masks
end

"""
    enumerate_constrained_states(n)

Return all 2D PXP constrained states on an n×n torus as Vector{UInt64} where
each state packs `n*n` bits: site (x, y) excited iff bit `y*n + x` is 1.
(We use the "excited == 1" mask convention internally.)

The enumeration composes constrained row masks vertically with the PXP
constraint that no two vertically adjacent sites are both excited, including
the y == n-1 <-> y == 0 wrap-around.
"""
function enumerate_constrained_states(n::Int)
    n <= 8 || error("n=$n too large; consider symmetry reduction or sector-restricted Krylov")
    row_masks = constrained_row_masks(n)
    say("  row patterns: $(length(row_masks))")
    # DFS over rows. State during DFS: list of chosen row masks so far.
    states = UInt64[]
    function dfs(rows::Vector{UInt64})
        if length(rows) == n
            # PBC closure: rows[n-1] and rows[0] must be vertically compatible.
            (rows[1] & rows[end]) == 0 || return
            packed = UInt64(0)
            for y = 1:n
                packed |= UInt64(rows[y]) << ((y - 1) * n)
            end
            push!(states, packed)
            return
        end
        last_row = isempty(rows) ? UInt64(0) : rows[end]
        for r in row_masks
            (r & last_row) == 0 || continue   # no two stacked excited
            push!(rows, r)
            dfs(rows)
            pop!(rows)
        end
    end
    # Optimization: also enforce row 1 vs row n PBC constraint by trying each
    # first row, and for each remember it so DFS can check the wrap at the end.
    # The DFS above already handles wraparound at the end of the recursion.
    dfs(UInt64[])
    return states
end

function neel_state(n::Int; excited_on::Symbol = :even)
    s = UInt64(0)
    for y = 0:(n - 1), x = 0:(n - 1)
        even = iseven(x + y)
        excited = excited_on === :even ? even : !even
        if excited
            s |= UInt64(1) << (y * n + x)
        end
    end
    return s
end

"""
    build_hamiltonian(states, n)

Sparse PXP H = sum_i X_i restricted to the supplied constrained basis.
"""
function build_hamiltonian(states::Vector{UInt64}, n::Int)
    nsites = n * n
    idx = Dict{UInt64,Int}()
    for (i, s) in enumerate(states)
        idx[s] = i
    end
    rows = Int[]
    cols = Int[]
    sizehint!(rows, length(states) * nsites ÷ 2)
    sizehint!(cols, length(states) * nsites ÷ 2)
    for (i, s) in enumerate(states), bit = 0:(nsites - 1)
        flipped = s ⊻ (UInt64(1) << bit)
        j = get(idx, flipped, 0)
        if j != 0
            push!(rows, i)
            push!(cols, j)
        end
    end
    return sparse(rows, cols, ones(Float64, length(rows)), length(states), length(states))
end

"""
    site_density_diagonal(states, n)

Return a length-N diagonal vector with entries equal to the average excitation
density (count of excited sites / nsites) per basis state.
"""
function global_density_diagonal(states::Vector{UInt64}, n::Int)
    nsites = n * n
    return Float64[count_ones(s) / nsites for s in states]
end

# ----------- Lanczos Krylov time evolution -----------

"""
    krylov_step!(psi, H, dt, m)

Advance `psi` by one Krylov step of duration `dt` using a Lanczos subspace of
dimension `m`. `H` must be Hermitian. Returns the updated psi (in-place).

Standard Lanczos for `exp(-i*dt*H) |psi>`: builds tridiagonal T, diagonalizes,
maps amplitudes back. Cheap for sparse H with O(N) nonzeros and modest m.
"""
function krylov_step!(psi::Vector{ComplexF64}, H::SparseMatrixCSC{Float64,Int},
                     dt::Float64, m::Int)
    N = length(psi)
    beta0 = norm(psi)
    @assert beta0 > 0
    V = Matrix{ComplexF64}(undef, N, m + 1)
    alpha = Vector{Float64}(undef, m)
    beta = Vector{Float64}(undef, m)
    V[:, 1] = psi / beta0
    j_max = m
    for j = 1:m
        w = H * @view(V[:, j])
        if j > 1
            w .-= beta[j - 1] .* @view(V[:, j - 1])
        end
        alpha[j] = real(dot(@view(V[:, j]), w))
        w .-= alpha[j] .* @view(V[:, j])
        # Reorthogonalize once for stability
        for k = 1:j
            c = dot(@view(V[:, k]), w)
            w .-= c .* @view(V[:, k])
        end
        beta[j] = norm(w)
        if beta[j] < 1e-12
            j_max = j
            break
        end
        if j < m
            V[:, j + 1] = w / beta[j]
        end
    end
    # Tridiagonal in real arithmetic
    T = SymTridiagonal(alpha[1:j_max], beta[1:(j_max - 1)])
    F = eigen(T)
    # exp(-i*dt*T) v_1 where v_1 = e_1*beta0
    coeffs = F.vectors' * (ComplexF64[i == 1 ? beta0 : 0.0 for i = 1:j_max])
    phased = coeffs .* exp.(-im .* dt .* F.values)
    a = F.vectors * phased
    fill!(psi, 0)
    @views for j = 1:j_max
        psi .+= a[j] .* V[:, j]
    end
    return psi
end

function main()
    n = parse(Int, get(ENV, "NEEL_KRYLOV_N", "6"))
    total_time = parse(Float64, get(ENV, "NEEL_KRYLOV_TOTAL_TIME", "3.0"))
    dt_meas = parse(Float64, get(ENV, "NEEL_KRYLOV_DT_MEAS", "0.1"))
    dt_step = parse(Float64, get(ENV, "NEEL_KRYLOV_DT_STEP", "0.05"))
    m = parse(Int, get(ENV, "NEEL_KRYLOV_M", "30"))
    output_path = get(ENV, "NEEL_KRYLOV_OUTPUT_JSON",
        joinpath(project_root, "artifacts", "neel_ed_krylov_$(n).json"))

    say("PXP single-Neel Krylov ED, n=$n, T=$total_time, dt_meas=$dt_meas, dt_step=$dt_step, m=$m")
    say("Enumerating constrained basis...")
    t_enum = @elapsed states = enumerate_constrained_states(n)
    say("  basis size = $(length(states)) (enumerated in $(round(t_enum, digits=2))s)")

    say("Building sparse H...")
    t_h = @elapsed H = build_hamiltonian(states, n)
    say("  nnz = $(nnz(H)) ($(round(nnz(H) / length(states), digits=1)) per row), built in $(round(t_h, digits=2))s")

    neel = neel_state(n; excited_on = :even)
    idx_neel = findfirst(==(neel), states)
    idx_neel === nothing && error("Neel state not in constrained basis (sanity check failed)")
    say("  Neel state at basis index $idx_neel of $(length(states))")

    dens_diag = global_density_diagonal(states, n)

    psi = zeros(ComplexF64, length(states))
    psi[idx_neel] = 1.0

    ts = Float64[0.0]
    densities = Float64[sum(abs2.(psi) .* dens_diag)]
    say(@sprintf("  t=%.2f n=%.6f", 0.0, densities[1]))

    save() = open(output_path, "w") do io
        JSON3.pretty(io, Dict(
            "cell_n" => n,
            "total_time" => total_time,
            "dt_meas" => dt_meas,
            "dt_step" => dt_step,
            "m" => m,
            "basis_size" => length(states),
            "nnz" => nnz(H),
            "ts" => ts,
            "density" => densities,
        ))
    end
    mkpath(dirname(output_path))
    save()

    # Time evolve up to total_time, measuring every dt_meas.
    t_now = 0.0
    next_meas = dt_meas
    while next_meas <= total_time + 1e-12
        # Step in dt_step chunks until we reach next_meas
        while t_now < next_meas - 1e-12
            step = min(dt_step, next_meas - t_now)
            krylov_step!(psi, H, step, m)
            t_now += step
        end
        push!(ts, next_meas)
        push!(densities, sum(abs2.(psi) .* dens_diag))
        say(@sprintf("  t=%.2f n=%.6f |psi|=%.6f", next_meas, densities[end], norm(psi)))
        save()
        next_meas += dt_meas
    end
    say("Done. Wrote $output_path")
end

main()
