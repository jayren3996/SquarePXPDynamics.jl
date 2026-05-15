using ITensors

@testset "square PEPS product state" begin
    psi = product_square_peps(2, 2; state = :down, maxdim = 3)

    @test length(psi.sites) == 4
    @test psi.maxdim == 3
    @test haskey(psi.tensors, SquareCoord(1, 1))
    @test dim(physical_index(psi, SquareCoord(1, 1))) == 2
    @test link_index(psi, SquareCoord(1, 1), :right) ==
          link_index(psi, SquareCoord(2, 1), :left)
    @test dim(link_index(psi, SquareCoord(1, 1), :right)) == 3

    T = site_tensor(psi, SquareCoord(1, 1))
    p = physical_index(psi, SquareCoord(1, 1))
    other = filter(!=(p), inds(T))
    @test T[p=>1, (i=>1 for i in other)...] == 0
    @test T[p=>2, (i=>1 for i in other)...] == 1

    up = product_square_peps(1, 1; state = :up)
    Tup = site_tensor(up, SquareCoord(1, 1))
    pup = physical_index(up, SquareCoord(1, 1))
    boundary = filter(i -> dim(i) == 1, inds(Tup))
    @test Tup[pup=>1, boundary[1]=>1, boundary[2]=>1, boundary[3]=>1, boundary[4]=>1] == 1
    @test Tup[pup=>2, boundary[1]=>1, boundary[2]=>1, boundary[3]=>1, boundary[4]=>1] == 0

    @test_throws ArgumentError product_square_peps(0, 1)
    @test_throws ArgumentError product_square_peps(1, 0)
    @test_throws ArgumentError product_square_peps(1, 1; state = :plus)
    @test_throws ArgumentError product_square_peps(1, 1; maxdim = 0)
end
