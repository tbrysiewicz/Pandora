using Test

macro smoketest(expr)
    return quote
        try
            $(esc(expr))
           true
        catch err
           false
        end
    end
end

@testset "Algorithm Tests" verbose=true begin


    for AD in keys(Pandora.ALGORITHM_DATA)
        @testset "Algorithm Data Tests for $(AD)" begin
            D = Pandora.ALGORITHM_DATA[AD]

            #Type Tests
            @test isa(D, AlgorithmDatum)
            @test isa(input_properties(D), Vector{EnumerativeProperty})
            @test isa(output_property(D), EnumerativeProperty)
            @test isa(reliability(D), Symbol)
            @test isa(name(D), String)
            @test isa(default_kwargs(D), Dict{Symbol,Any})

            #Smoke Tests
            T = TwentySevenLines()
            @smoketest(
                Pandora.learn!(T, output_property(D); algorithm = AD)
            )
            
        end
    end


end