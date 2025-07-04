using Test

@testset "Automated Algorithm Tests" verbose=true begin


    for AD in keys(Pandora.ALGORITHM_DATA)
        @testset "Algorithm Data Tests for $(AD)" begin
            D = Pandora.ALGORITHM_DATA[AD]
            
            #Type Tests
            @test isa(D, AlgorithmDatum)
            @test isa(input_properties(D), Vector{EnumerativeProperty})
            @test isa(output_property(D), EnumerativeProperty)
            @test isa(reliability(D), Symbol)

            if D.automated
                #Smoke Tests
                @test try
                    T = TwentySevenLines()
                    Pandora.learn!(T, output_property(D); algorithm = AD)
                    true
                catch
                    false
                end
            end
            
        end
    end


end
 
@testset "Subtle Core Code Tests" verbose=true begin
    #Checks that kwargs get passed through the knowledge tree
    @testset "Galois group generators count" begin
        T = TwentySevenLines()
        G = galois_group(T; n_monodromy_loops=3)
        @test length(monodromy_sample(T)) == 3
    end    
end

#Unit Tests
@testset "Unit Tests" verbose=true begin 
    @testset "Degree bounds for TwentySevenLines" begin
        T = TwentySevenLines()
        @test bkk_bound(T) == 45
        @test bezout_bound(T) == 81
        @test degree(T) == 27
    end
    @testset "Galois Group for TwentySevenLines" begin
        T = TwentySevenLines()
        G = galois_group(T)
        @test is_transitive(G)
        @test order(G) == 51840
        @test is_primitive(G)
    end
end


@testset "Draw Newton Polygons" begin
    # Smoke test for drawing Newton polygons of an enumerative problem
    @test try
        # Variable and parameter setup
        @var x y
        @var a b c

        # Define polynomials
        f1 = a + b*x + c*y + b*c*x*y
        f2 = b + c*x^2

        # Create system and enumerative problem
        F = System([f1, f2], variables=[x, y], parameters=[a, b, c])
        E = EnumerativeProblem(F)

        # Compute Newton polytopes
        NP = newton_polytopes(E)

        # Attempt to visualize the first Newton polytope
        visualize(NP[1])
        true
    catch
        false
    end
end


#=
#todo Tests
@testset "Todo Tests" verbose=true begin
    @testset "Monodromy todo for TwentySevenLines" begin
        T = TwentySevenLines()
        @test try
            T = TwentySevenLines()
            @test monodromy_basis(T)
            true
        catch
            false
        end
        @test try
            T = TwentySevenLines()
            @test monodromy_reorder!(T)
            true
        catch
            false
        end
        @test try
            T = TwentySevenLines()
            @test monodromy_coordinates(T)
            true
        catch
            false
        end
    end
end
=#