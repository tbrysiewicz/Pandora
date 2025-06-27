using Test

@testset "Galois group generators count" begin
    T = TwentySevenLines()
    G = galois_group(T; n_monodromy_loops=3)
    @test length(gens(G)) == 3
end

@testset "Degree bounds for TwentySevenLines" begin
    T = TwentySevenLines()
    @test bkk_bound(T) == 45
    @test bezout_bound(T) == 81
    @test degree(T) == 27
end