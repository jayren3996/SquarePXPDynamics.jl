using ITensors
using LinearAlgebra

function _d2_dense_index(values)
    idx = 1
    nsites = length(values)
    for (site, value) in enumerate(values)
        1 <= value <= 2 || throw(ArgumentError("basis values must be 1 or 2"))
        idx += (value - 1) * 2^(nsites - site)
    end
    return idx
end

function _d2_star_coords(cell, center)
    c = wrap(cell, center)
    return (
        center = c,
        right = neighbor(cell, c, :right),
        up = neighbor(cell, c, :up),
        left = neighbor(cell, c, :left),
        down = neighbor(cell, c, :down),
    )
end

function _d2_dense_all_down(cell)
    values = ntuple(_ -> 2, length(cell.reps))
    state = zeros(ComplexF64, 2^length(cell.reps))
    state[_d2_dense_index(values)] = 1
    return state
end

function _d2_apply_dense_star(state, cell, center, step)
    sites = Tuple(_d2_star_coords(cell, center))
    positions = ntuple(site -> findfirst(==(sites[site]), cell.reps), length(sites))
    all(!isnothing, positions) || throw(ArgumentError("star sites must be in cell reps"))
    gate = projected_square_pxp_gate(step; evolution = :real)
    next_state = zeros(ComplexF64, length(state))

    for in_values in Iterators.product((1:2 for _ = 1:length(cell.reps))...)
        in_idx = _d2_dense_index(in_values)
        amplitude = state[in_idx]
        iszero(amplitude) && continue
        star_in = ntuple(site -> in_values[positions[site]], length(sites))
        star_in_idx = _d2_dense_index(star_in)
        for star_out in Iterators.product((1:2 for _ = 1:length(sites))...)
            star_out_idx = _d2_dense_index(star_out)
            out_values = collect(in_values)
            for site = 1:length(sites)
                out_values[positions[site]] = star_out[site]
            end
            out_idx = _d2_dense_index(Tuple(out_values))
            next_state[out_idx] += gate[star_out_idx, star_in_idx] * amplitude
        end
    end
    return next_state
end

function _d2_dense_density(state, cell)
    normsq = sum(abs2, state)
    normsq > 0 || throw(ArgumentError("dense state has zero norm"))
    density = 0.0
    for values in Iterators.product((1:2 for _ = 1:length(cell.reps))...)
        idx = _d2_dense_index(values)
        density += count(==(1), values) * abs2(state[idx])
    end
    return density / (length(cell.reps) * normsq)
end

function _d2_exact_serial_star_densities(cell, step)
    state = _d2_dense_all_down(cell)
    densities = Float64[]
    for center in cell.reps
        state = _d2_apply_dense_star(state, cell, center, step)
        push!(densities, _d2_dense_density(state, cell))
    end
    return densities, state
end

function _d2_absorb_all_weights_once(psi)
    tensors = Dict(c => copy(T) for (c, T) in psi.tensors)
    for key in keys(psi.link_weights)
        tensors[key.site] = absorb_link_weight(tensors[key.site], psi, key.site, key.dir)
    end
    return tensors
end

function _d2_dense_contract_ipeps(psi)
    cell = psi.unitcell
    length(cell.reps) == 9 || throw(ArgumentError("test helper is limited to 3x3 cells"))
    tensors = _d2_absorb_all_weights_once(psi)
    theta = ITensor()
    started = false
    for c in cell.reps
        if !started
            theta = tensors[c]
            started = true
        else
            theta = @disable_warn_order theta * tensors[c]
        end
    end
    phys = Tuple(physical_index(psi, c) for c in cell.reps)
    for p in phys
        hasind(theta, p) || throw(ArgumentError("dense contraction lost physical index $p"))
    end
    state = zeros(ComplexF64, 2^length(cell.reps))
    for values in Iterators.product((1:2 for _ = 1:length(cell.reps))...)
        state[_d2_dense_index(values)] = theta[(phys[i] => values[i] for i = 1:length(phys))...]
    end
    return state
end

function _d2_trace_serial_star(; maxdim_d2 = 2, step = 0.02, cutoff = 1e-12)
    cell = PeriodicSquareUnitCell(3, 3)
    reference = _d2_dense_all_down(cell)
    psi_d1 = product_square_ipeps(cell; state = :down, maxdim = 1)
    psi_d2 = product_square_ipeps(cell; state = :down, maxdim = 1)
    records = NamedTuple[]

    for (star_index, center) in enumerate(cell.reps)
        reference = _d2_apply_dense_star(reference, cell, center, step)
        info_d1 = project_star!(
            psi_d1,
            center,
            step;
            evolution = :real,
            projected = true,
            maxdim = 1,
            cutoff,
        )
        info_d2 = project_star!(
            psi_d2,
            center,
            step;
            evolution = :real,
            projected = true,
            maxdim = maxdim_d2,
            cutoff,
        )
        dense_d1 = _d2_dense_contract_ipeps(psi_d1)
        dense_d2 = _d2_dense_contract_ipeps(psi_d2)
        simple_d1 = measure_simple(psi_d1)
        simple_d2 = measure_simple(psi_d2)
        push!(
            records,
            (
                star_index,
                center,
                reference_density = _d2_dense_density(reference, cell),
                dense_d1_density = _d2_dense_density(dense_d1, cell),
                dense_d2_density = _d2_dense_density(dense_d2, cell),
                simple_d1_density = simple_d1.density,
                simple_d2_density = simple_d2.density,
                d2_log_norm = log_norm(psi_d2),
                d2_norm_factors = copy(info_d2.norm_factors),
                d2_keptdims = copy(info_d2.keptdims),
                d2_min_lambda = copy(info_d2.min_lambda),
                d2_max_truncerr = info_d2.max_truncerr,
                d1_info = info_d1,
                d2_info = info_d2,
            ),
        )
    end
    return records
end

function _d2_ipeps_after_serial_stars(nstars; maxdim = 2, step = 0.02, cutoff = 1e-12)
    cell = PeriodicSquareUnitCell(3, 3)
    0 <= nstars <= length(cell.reps) || throw(ArgumentError("nstars must be in 0:9"))
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    last_info = nothing
    for center in cell.reps[1:nstars]
        last_info = project_star!(
            psi,
            center,
            step;
            evolution = :real,
            projected = true,
            maxdim,
            cutoff,
        )
    end
    return psi, last_info
end

@testset "3x3 PXP dense serial-star reference" begin
    cell = PeriodicSquareUnitCell(3, 3)
    @test cell.reps == [SquareCoord(x, y) for y = 1:3 for x = 1:3]

    densities, final_state = _d2_exact_serial_star_densities(cell, 0.02)

    @test length(densities) == 9
    @test _d2_dense_density(final_state, cell) ≈ 0.00039962698926202146 atol = 1e-15
end

@testset "3x3 product iPEPS dense contraction" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test _d2_dense_density(_d2_dense_contract_ipeps(down), cell) ≈ 0.0 atol = 1e-15
    @test _d2_dense_density(_d2_dense_contract_ipeps(up), cell) ≈ 1.0 atol = 1e-15
end

@testset "D2 serial PXP first-divergence localization" begin
    records = _d2_trace_serial_star()

    @test records[1].dense_d1_density ≈ records[1].reference_density atol = 5e-8
    @test records[1].dense_d2_density ≈ records[1].reference_density atol = 5e-8
    @test records[1].simple_d2_density ≈ records[1].dense_d2_density atol = 5e-8
    @test records[end].reference_density ≈ 0.00039962698926202146 atol = 1e-15
    @test records[end].dense_d2_density ≈ records[end].reference_density atol = 5e-15
    @test records[end].simple_d2_density ≈ 0.00021094978264193 atol = 5e-15

    first_dense_divergence = findfirst(
        r -> !isapprox(r.dense_d2_density, r.reference_density; atol = 5e-8, rtol = 0),
        records,
    )
    first_simple_divergence = findfirst(
        r -> !isapprox(r.simple_d2_density, r.dense_d2_density; atol = 5e-8, rtol = 0),
        records,
    )

    @test first_dense_divergence === nothing
    @test first_simple_divergence == 3
    divergent = records[first_simple_divergence]
    @test divergent.center == SquareCoord(3, 1)
    @test divergent.reference_density ≈ 0.00013326224449912612 atol = 5e-15
    @test divergent.dense_d2_density ≈ divergent.reference_density atol = 5e-15
    @test divergent.simple_d2_density ≈ 0.000111054099003352 atol = 5e-15
    @test divergent.d2_log_norm ≈ 1.7328679513998648 atol = 5e-15
    @test_broken divergent.simple_d2_density ≈ divergent.dense_d2_density atol = 5e-8
end

@testset "D2 first divergent state measurement-only regression" begin
    psi, info = _d2_ipeps_after_serial_stars(3)
    dense_density = _d2_dense_density(_d2_dense_contract_ipeps(psi), psi.unitcell)
    simple_density = density_simple(psi)

    @test info.keptdims[:left] == 2
    @test info.keptdims[:right] == 2
    @test info.max_truncerr < 1e-24
    @test dense_density ≈ 0.00013326224449912612 atol = 5e-15
    @test simple_density ≈ 0.000111054099003352 atol = 5e-15
    @test_broken simple_density ≈ dense_density atol = 5e-8
end
