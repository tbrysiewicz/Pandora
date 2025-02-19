mutable struct Knowledge
    explanations :: Dict{Symbol,String}

    function Knowledge()
        Kn = new()
        Kn.explanations = Dict{Symbol,String}()
        return(Kn)
    end
end


function justify(EP::EnumerativeProblem,fact::Symbol,justification::String)
    (knowledge(EP).explanations)[fact] = justification
end