using ITensors

function _dense_square_star_index(values)
    idx = 1
    for (site, value) in enumerate(values)
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

function _itensor_gate_entry(G, sites, out_values, in_values)
    out = prime.(sites)
    return G[(out[i] => out_values[i] for i in eachindex(sites))...,
             (sites[i] => in_values[i] for i in eachindex(sites))...]
end

function _itensor_basis_ket(sites, values)
    ket = ITensor(ComplexF64, sites...)
    ket[(sites[i] => values[i] for i in eachindex(sites))...] = 1.0
    return ket
end

function _dense_from_itensor_gate(G, sites)
    dense = zeros(ComplexF64, 2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES)
    for out_values in Iterators.product((1:2 for _ in 1:SQUARE_STAR_SITES)...)
        out_idx = _dense_square_star_index(out_values)
        for in_values in Iterators.product((1:2 for _ in 1:SQUARE_STAR_SITES)...)
            in_idx = _dense_square_star_index(in_values)
            dense[out_idx, in_idx] = _itensor_gate_entry(G, sites, out_values, in_values)
        end
    end
    return dense
end

function _opposite_dir(dir)
    dir === :right && return :left
    dir === :up && return :down
    dir === :left && return :right
    dir === :down && return :up
    error("unreachable")
end

function _filled_site_tensor(psi, c)
    p = physical_index(psi, c)
    left = link_index(psi, c, :left)
    right = link_index(psi, c, :right)
    up = link_index(psi, c, :up)
    down = link_index(psi, c, :down)
    T = ITensor(ComplexF64, p, left, right, up, down)
    for pv in 1:dim(p), lv in 1:dim(left), rv in 1:dim(right), uv in 1:dim(up), dv in 1:dim(down)
        T[p => pv, left => lv, right => rv, up => uv, down => dv] =
            complex(100pv + 10lv + rv, 10uv + dv)
    end
    return T
end

@testset "safe iPEPS link weight helpers" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    c = SquareCoord(1, 1)

    lambda = link_weight(psi, c, :right)
    @test lambda == [1.0, 0.0]
    lambda[1] = 0.25
    @test link_weight(psi, c, :right) == [1.0, 0.0]

    values_by_dir = Dict(
        :right => [0.6, 0.8],
        :up => [0.3, 0.7],
        :left => [0.2, 0.9],
        :down => [0.4, 0.5],
    )
    for dir in (:right, :up, :left, :down)
        values = values_by_dir[dir]
        set_link_weight!(psi, c, dir, values)
        @test link_weight(psi, c, dir) == values
        @test link_weight(psi, neighbor(cell, c, dir), _opposite_dir(dir)) == values

        link = link_index(psi, c, dir)
        Lambda = link_weight_tensor(psi, c, dir)
        @test inds(Lambda) == (link, prime(link))
        for i in 1:dim(link), j in 1:dim(link)
            @test Lambda[link => i, prime(link) => j] == (i == j ? values[i] : 0.0)
        end
    end

    @test_throws ArgumentError link_weight(psi, c, :diagonal)
    @test_throws ArgumentError set_link_weight!(psi, c, :right, [1.0])
    @test_throws ArgumentError set_link_weight!(psi, c, :right, [1.0, -0.1])
    @test_throws ArgumentError set_link_weight!(psi, c, :right, [1.0, Inf])

    psi.link_weights[bondkey(cell, c, :right)] = [1.0]
    @test_throws ArgumentError link_weight_tensor(psi, c, :right)
end

@testset "S3 link weight absorption helpers" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    c = SquareCoord(2, 2)
    T = _filled_site_tensor(psi, c)

    for dir in (:right, :up, :left, :down)
        set_link_weight!(psi, c, dir, [2.0, 4.0])
        absorbed = absorb_link_weight(T, psi, c, dir)
        restored = deabsorb_link_weight(absorbed, psi, c, dir)
        @test inds(absorbed) == inds(T)
        @test norm(restored - T) < 1e-12
    end

    set_link_weight!(psi, c, :right, [2.0, 0.0])
    deabsorbed = deabsorb_link_weight(T, psi, c, :right)
    expected = ITensor(ComplexF64, inds(T)...)
    p = physical_index(psi, c)
    left = link_index(psi, c, :left)
    right = link_index(psi, c, :right)
    up = link_index(psi, c, :up)
    down = link_index(psi, c, :down)
    for pv in 1:dim(p), lv in 1:dim(left), rv in 1:dim(right), uv in 1:dim(up), dv in 1:dim(down)
        expected[p => pv, left => lv, right => rv, up => uv, down => dv] =
            rv == 1 ? T[p => pv, left => lv, right => rv, up => uv, down => dv] / 2 : 0.0
    end
    @test norm(deabsorbed - expected) < 1e-12

    @test_throws ArgumentError absorb_link_weight(ITensor(physical_index(psi, c)), psi, c, :right)
end

@testset "link weight entropy diagnostics" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    c = SquareCoord(1, 1)

    @test weight_entropy([1.0, 0.0]) ≈ 0.0
    @test weight_entropy([1.0, 1.0]) ≈ log(2)
    @test_throws ArgumentError weight_entropy([-1.0, 1.0])
    @test_throws ArgumentError weight_entropy([0.0, 0.0])

    set_link_weight!(psi, c, :up, [1 / sqrt(2), 1 / sqrt(2)])
    @test bond_entropy(psi, c, :up) ≈ log(2)
    @test bond_entropy(psi, neighbor(cell, c, :up), :down) ≈ log(2)

    entropies = all_bond_entropies(psi)
    @test length(entropies) == 2 * cell.Lx * cell.Ly
    @test entropies[bondkey(cell, c, :up)] ≈ log(2)
end

@testset "ITensor square-star PXP gates match dense gates" begin
    sites = [Index(2, tag) for tag in ("center", "right", "up", "left", "down")]
    t = 0.17

    for (G, dense) in (
        (square_pxp_gate_itensor(sites, t; evolution = :real),
         square_pxp_gate(t; evolution = :real)),
        (projected_square_pxp_gate_itensor(sites, t; evolution = :imaginary),
         projected_square_pxp_gate(t; evolution = :imaginary)),
    )
        @test inds(G) == (prime.(sites)..., sites...)
        @test _dense_from_itensor_gate(G, sites) ≈ dense
        for in_values in Iterators.product((1:2 for _ in 1:SQUARE_STAR_SITES)...)
            in_idx = _dense_square_star_index(in_values)
            out_tensor = G * _itensor_basis_ket(sites, in_values)
            for out_values in Iterators.product((1:2 for _ in 1:SQUARE_STAR_SITES)...)
                out_idx = _dense_square_star_index(out_values)
                out = prime.(sites)
                @test out_tensor[(out[i] => out_values[i] for i in eachindex(sites))...] ≈ dense[out_idx, in_idx]
            end
        end
    end

    @test_throws ArgumentError square_pxp_gate_itensor(sites[1:4], t)
    @test_throws ArgumentError square_pxp_gate_itensor([Index(3, "bad"), sites[2:end]...], t)
    @test_throws ArgumentError square_pxp_gate_itensor(sites, t; evolution = :thermal)
end
