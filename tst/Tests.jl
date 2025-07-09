using Test

@testset "All algorithms" verbose=true begin 

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


@testset "Readme tests" begin
    # 1. Define variables and parameters
    @var x y
    @var a[1:4] b[1:4]

    # 2. Define your system
    f1 = a[1]*x + a[2]*y + a[3]*x^2*y + a[4]*x*y^2
    f2 = b[1]*x + b[2]*y + b[3]*x^2*y + b[4]*x*y^2
    E = EnumerativeProblem([f1, f2], variables = [x, y], parameters = vcat(a, b), torus_only=true)

    # 3. Compute invariants and knowledge
    @test degree(E) == 4
    @test bkk_bound(E) == 4
    @test bezout_bound(E) == 9
    K = knowledge(E)
    @test length(K) > 0
    # Output of last(knowledge(E)) is a KnowledgeNode, structure may vary #[unknown output]

    # 4. Explore group structure
    G = monodromy_group(E)
    @test is_primitive(G) == false
    @test order(G) == 8
    gensG = gens(G)
    @test length(gensG) == 2
    @test is_decomposable(E) == true
    @test is_lacunary(E) == true

    # 5. Sample reality and optimize
    ER = explore_reality(E; n_samples=1000)
    @test isa(ER, Vector)
    # Output of ER is a Dict with integer keys and values #[unknown output]

    O = maximize_n_real_solutions(E)
    @test !isnothing(O)
    C = certify(record_fibre(O), E)
    @test C isa Pandora.CertificationResult

    # 6. Visualize and save discriminant
    V, P = visualize(E; near = record_parameters(O), strategy = :quadtree)
    @test !isnothing(V)
    @test !isnothing(P)
    save(P, "OutputFiles/MyDiscriminant.png")
    @test isfile("OutputFiles/MyDiscriminant.png") #[unknown output if path changes]

    # 7. Refine and save improved visualization
    refine!(V)
    refine!(V)
    P2 = visualize(V)
    save(P2, "OutputFiles/MyBetterDiscriminant.png")
    @test isfile("OutputFiles/MyBetterDiscriminant.png") #[unknown output if path changes]

    # 8. Visualize supports and Newton polytopes
    SupportVisualization = visualize_support(E)
    save(SupportVisualization[1], "OutputFiles/support1.png")
    save(SupportVisualization[2], "OutputFiles/support2.png")
    @test isfile("OutputFiles/support1.png") #[unknown output if path changes]
    @test isfile("OutputFiles/support2.png") #[unknown output if path changes]

    NP = newton_polytopes(E)
    @test length(NP) == 2
    mv = volume(sum(NP)) - volume(NP[1]) - volume(NP[2])
    @test mv == 4

    # 9. Automate knowledge and summarize
    T = TwentySevenLines()
    automate!(T)
    summarize(T)
    @test isfile("OutputFiles/latex_summary.pdf") #[unknown output if path changes]
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