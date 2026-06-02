using LinearAlgebra
using SparseArrays
using Printf

# ----------------------------------------------------------------------
# Single-site 2x2 operators. Convention: |up>=|0> (index 1), |down>=|1> (index 2).
# h_c = X_center * P_down(right) * P_down(up) * P_down(left) * P_down(down),
# matching src/PXPModel.jl square_pxp_star_hamiltonian().
# ----------------------------------------------------------------------
const I2 = sparse(ComplexF64[1 0; 0 1])
const X  = sparse(ComplexF64[0 1; 1 0])
const Pd = sparse(ComplexF64[0 0; 0 1])   # projector onto |down>

# Embed a Dict{site=>2x2 op} into an N-site Hilbert space (site 1 = most significant).
function embed(ops::Dict{Int,SparseMatrixCSC{ComplexF64,Int}}, N::Int)
    out = haskey(ops, 1) ? ops[1] : I2
    for s in 2:N
        out = kron(out, haskey(ops, s) ? ops[s] : I2)
    end
    return out
end

fnorm(A) = sqrt(real(sum(abs2, A)))
opnorm_dense(A) = opnorm(Matrix(A))   # spectral norm

# ======================================================================
# OPEN patch (no wrap). Used for commutator norms: place BOTH centers as
# interior cells (all 4 neighbors in-bounds) so each h_c is a full bulk
# 5-site star. The commutator is a purely local property, so an open
# patch with interior centers reproduces the bulk operator exactly.
# ======================================================================
struct Patch; Lx::Int; Ly::Int; end
pnsites(p::Patch) = p.Lx * p.Ly
plin(p::Patch, x::Int, y::Int) = (y * p.Lx + x) + 1   # zero-based coords
inb(p::Patch, x, y) = (0 <= x < p.Lx) && (0 <= y < p.Ly)

function interior(p::Patch, x::Int, y::Int)
    return inb(p,x+1,y) && inb(p,x-1,y) && inb(p,x,y+1) && inb(p,x,y-1)
end

function h_center_open(p::Patch, x::Int, y::Int)
    N = pnsites(p)
    @assert interior(p, x, y) "center ($x,$y) must be interior so all 4 neighbors exist"
    c = plin(p, x, y)
    ops = Dict{Int,SparseMatrixCSC{ComplexF64,Int}}()
    ops[c] = X
    for (nx, ny) in ((x+1,y),(x,y+1),(x-1,y),(x,y-1))
        s = plin(p, nx, ny)
        ops[s] = haskey(ops,s) ? ops[s]*Pd : Pd
    end
    return embed(ops, N)
end

comm(A,B) = A*B - B*A

println("="^72)
println("PXP commutation: h_c = X_c * P_down(right)*P_down(up)*P_down(left)*P_down(down)")
println("="^72)

# Helper: build a patch big enough so both centers are interior, return ||[h1,h2]||.
function commnorm(ca::Tuple{Int,Int}, cb::Tuple{Int,Int})
    # bounding box of both stars, pad by 1 so both centers are interior
    xs = [ca[1],cb[1]]; ys=[ca[2],cb[2]]
    minx=minimum(xs)-1; maxx=maximum(xs)+1
    miny=minimum(ys)-1; maxy=maximum(ys)+1
    Lx = maxx-minx+1; Ly = maxy-miny+1
    p = Patch(Lx,Ly)
    a=(ca[1]-minx, ca[2]-miny); b=(cb[1]-minx, cb[2]-miny)
    h1=h_center_open(p,a...); h2=h_center_open(p,b...)
    return fnorm(comm(h1,h2)), pnsites(p)
end

# ---- Part 1: commutator norms at each separation class ----
println("\n--- Part 1: ||[h_c, h_c']|| (Frobenius) by center separation ---\n")

cases = [
    ("identical (same center)",                (3,3), (3,3)),
    ("nearest-neighbor dist 1 (horizontal)",   (3,3), (4,3)),
    ("nearest-neighbor dist 1 (vertical)",     (3,3), (3,4)),
    ("diagonal NNN dist sqrt2 (shares 2 proj)",(3,3), (4,4)),
    ("diagonal NNN dist sqrt2 (other diag)",   (3,3), (4,2)),
    ("straight dist 2 (shares 1 proj, horiz)", (3,3), (5,3)),
    ("straight dist 2 (shares 1 proj, vert)",  (3,3), (3,5)),
    ("knight (dx=2,dy=1) dist sqrt5",          (3,3), (5,4)),
    ("dist 3 straight, disjoint support",      (3,3), (6,3)),
]
# NOTE: the meaningful cases are the OVERLAPPING-P-support pairs (diagonal NNN
# shares 2 projectors, straight-dist-2 shares 1) — they show commutation holds
# even when P-supports overlap, only failing at NN where the X collides. Far
# (dist>=3) pairs trivially commute (disjoint support); we keep one cheap dist-3
# check. Pairs with dist>=3 in BOTH axes would pad to a >=30-site patch (2^30+
# sparse, infeasible) and prove nothing new, so they are omitted.

results = Tuple{String,Float64,Int}[]
for (label, ca, cb) in cases
    nrm, N = commnorm(ca, cb)
    push!(results, (label, nrm, N))
    @printf "%-44s  ||[h,h']||_F = %.3e   (N=%d sites)\n" label nrm N
end

# Summary verdict
nn = [r[2] for r in results if occursin("nearest-neighbor", r[1])]
notnn = [r[2] for r in results if !occursin("nearest-neighbor", r[1]) && !occursin("identical", r[1])]
tol = 1e-12
all_nn_nonzero = all(>(1e-6), nn)
all_notnn_zero = all(<(tol), notnn)
ident_zero = results[1][2] < tol
println()
@printf "nearest-neighbor commutators: min=%.3e  (expect > 0)\n" minimum(nn)
@printf "non-NN (incl diag/dist2/dist>=3) commutators: max=%.3e (expect ~0)\n" maximum(notnn)
@printf "identical center commutator: %.3e (trivially 0)\n" results[1][2]
println("commute iff NOT nearest-neighbor: ", all_nn_nonzero && all_notnn_zero && ident_zero)

# ======================================================================
# Part 2: factorization on a PERIODIC lattice.
# exp(-i H_A dt) == prod_{c in A} exp(-i h_c dt) for checkerboard A.
# We need a torus where A-sublattice centers are mutually non-NN AND each
# h_c is a full periodic star. Use a periodic lattice and pick A = {(x+y) even}.
# Dense exp limits N: use a 3x4 = 12-site torus (2^12 = 4096) for the dense check.
# ======================================================================
println("\n" * "="^72)
println("--- Part 2: factorization exp(-i H_A dt) =?= prod_c exp(-i h_c dt) ---")
println("="^72)

struct Torus; Lx::Int; Ly::Int; end
tnsites(t::Torus) = t.Lx*t.Ly
tlin(t::Torus, x::Int, y::Int) = (mod(y,t.Ly)*t.Lx + mod(x,t.Lx)) + 1

function h_center_torus(t::Torus, x::Int, y::Int)
    N = tnsites(t)
    c = tlin(t,x,y)
    ops = Dict{Int,SparseMatrixCSC{ComplexF64,Int}}()
    # Build neighbor list; on small tori a neighbor may coincide with center or another
    # neighbor -> accumulate by multiplication to be exact.
    push_op!(s, op) = (ops[s] = haskey(ops,s) ? ops[s]*op : op)
    push_op!(c, X)
    for (nx,ny) in ((x+1,y),(x,y+1),(x-1,y),(x,y-1))
        push_op!(tlin(t,nx,ny), Pd)
    end
    return embed(ops, N)
end

function run_factorization(t::Torus, dt::Float64)
    N = tnsites(t)
    # checkerboard A: (x+y) even
    A = [(x,y) for y in 0:t.Ly-1 for x in 0:t.Lx-1 if iseven(x+y)]
    # verify A-centers are pairwise non-NN on THIS torus (so the claim's premise holds)
    function torus_NN(a,b)
        dx = mod(a[1]-b[1]+t.Lx, t.Lx); dx=min(dx,t.Lx-dx)
        dy = mod(a[2]-b[2]+t.Ly, t.Ly); dy=min(dy,t.Ly-dy)
        return (dx,dy)==(1,0) || (dx,dy)==(0,1)
    end
    bad = [(a,b) for a in A for b in A if a!=b && torus_NN(a,b)]
    HA = spzeros(ComplexF64, 2^N, 2^N)
    hs = Matrix{ComplexF64}[]
    for c in A
        h = h_center_torus(t, c...)
        HA += h
        push!(hs, Matrix(h))
    end
    # exp(-i HA dt)
    U_sum = exp(-im*dt*Matrix(HA))
    # prod_c exp(-i h_c dt) (order should not matter since they commute; use list order)
    U_prod = Matrix{ComplexF64}(I, 2^N, 2^N)
    for h in hs
        U_prod = U_prod * exp(-im*dt*h)
    end
    disc = opnorm(U_sum - U_prod)
    # also check pairwise commutation within A directly
    maxcomm = 0.0
    for i in 1:length(hs), j in i+1:length(hs)
        maxcomm = max(maxcomm, fnorm(hs[i]*hs[j]-hs[j]*hs[i]))
    end
    return A, bad, disc, maxcomm
end

using Printf
# Dense exp(-iHt) needs a 2^N x 2^N matrix; N=16 (4x4 torus) is 65536^2 ~ 68 GB,
# infeasible. The factorization is EXACTLY equivalent to pairwise within-A
# commutation, which Part 1 already proves for every separation class that occurs
# in a 4x4 checkerboard. Here we directly confirm the identity on small even-width
# tori (clean, seam-free checkerboard) that fit a dense exp.
for (Lx,Ly) in ((4,2),(6,2))
    t = Torus(Lx,Ly)
    N = tnsites(t)
    if N > 14
        println("skip $(Lx)x$(Ly) torus (N=$N too large for dense exp)")
        continue
    end
    for dt in (0.1, 0.7, 1.3)
        A, bad, disc, maxcomm = run_factorization(t, dt)
        @printf "Torus %dx%d (N=%2d) dt=%.2f: |A|=%2d  A-NN-pairs(should be 0)=%d  max within-A ||[h,h']||=%.2e  ||exp(-iHAdt)-prod||_op=%.3e\n" Lx Ly N dt length(A) length(bad)÷2 maxcomm disc
    end
end

println("\nDONE")
