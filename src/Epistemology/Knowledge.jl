mutable struct Knowledge
    explanations :: Dict{Symbol,String}

    function Knowledge()
        Kn = new()
        Kn.explanations = Dict{Symbol,String}()
        return(Kn)
    end
end
