# Test-only utilities for finite-cluster reference computations.
# Not part of the package; kept under test/ and included from runtests.jl.

using ITensors
using LinearAlgebra
using TriangularPEPSDynamics

"""
    cluster_vector_from_state(state, center) -> Vector{ComplexF64}

Contract the 7 star tensors at `center` (center first, then 6 directional
neighbors) into a single 128-dim vector indexed in the same site order as
the dense 7-site Hamiltonian (center physical index first, then directions
1..6 in `TRIANGULAR_DIRECTIONS` order).

Implementation note: this contracts the cluster as if it were finite (no
wrap-around). Each external bond of a star site is contracted with the
neighboring iPEPS tensor *only if* that neighbor is one of the other star
sites; external bonds dangling out of the cluster are summed over with the
all-ones vector — i.e., we trace the lambdas absorbed in the actual iPEPS
into the dangling legs. For tests that compare cluster vectors before and
after a gate, this is consistent because we use the same convention on both
sides.
"""
function cluster_vector_from_state(state::TriangularIPEPS, center::Coord)
    star = star_sites(center)
    reps = [wrap_coord(state.unitcell, sc) for sc in star]
    site_tensors = [state.tensors[rep] for rep in reps]
    phys = [state.phys_inds[rep] for rep in reps]

    # Build per-position physical index renames so each star position has a
    # unique physical leg even when reps are shared across positions.
    fresh_phys = [Index(2, "phys_pos_$(i)") for i in 1:7]
    renamed = [replaceind(site_tensors[i], phys[i], fresh_phys[i]) for i in 1:7]

    # Sum over all bond legs that are NOT shared between two star positions
    # by contracting with a vector of ones. Shared bonds contract internally.
    star_bond_keys = Set{Tuple{Coord,Int}}()
    for (pos_i, sc_i) in enumerate(star), d in 1:6
        nbr = neighbor(sc_i, d)
        for (pos_j, sc_j) in enumerate(star)
            if pos_j != pos_i && nbr == sc_j
                push!(star_bond_keys, (reps[pos_i], d))
            end
        end
    end

    for (pos_i, sc_i) in enumerate(star), d in 1:6
        if (reps[pos_i], d) in star_bond_keys
            continue
        end
        bind = bond_index(state, sc_i, d)
        ones_vec = ITensor(ones(ComplexF64, dim(bind)), bind)
        renamed[pos_i] = renamed[pos_i] * ones_vec
    end

    cluster = renamed[1]
    for k in 2:7
        cluster = cluster * renamed[k]
    end

    # cluster now has 7 fresh physical legs and zero virtual legs. Reshape
    # to a 128-dim vector in the canonical site order.
    return ComplexF64.(reshape(array(cluster, fresh_phys...), 128))
end
