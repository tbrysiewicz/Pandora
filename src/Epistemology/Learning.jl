
function knows(EP::EnumerativeProblem, EProp::EnumerativeProperty)

    K = knowledge(EP)
    prop_indices = findall(p->p.property==EProp,K)
    if length(prop_indices)==0
        return(false)
    else
        return(true)
    end
end

function know!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}, value::T) where T
    EA = user_given(EProp,value)
    learn!(EP,EProp; algorithm = EA)
    return()
end

function learn!(EP::EnumerativeProblem, EProp::EnumerativeProperty{T}; algorithm = nothing, user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}()) where T
    if algorithm == nothing
        error("Taylor needs to code how to find best algorithm")
    end
    EA = algorithm
    f = core_function(EA)
    for i in input_properties(EA)
        if knows(EP,i)==false
            error("Property \n    ["*string(i)*"]\nis required for algorithm \n    ["*name(EA)*"]\nbut is not known")
        end
    end
    input_knowledge = [get_knowledge(i,EP) for i in input_properties(EA)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    input_parameters = parameters(EA)
    for k in keys(user_parameters)
        if haskey(k,input_parameters)
            input_parameters[k] = user_parameters[k]
        end
    end
    flattened_dict = [input_parameters[k] for k in keys(input_parameters)]
	o = f(input_knowledge_values...;flattened_dict...)
    new_knowledge = Knowledge{T}(EProp,o,input_knowledge,input_parameters,EA)
	return(discover!(EP,new_knowledge))
end

function discover!(EP::EnumerativeProblem,K::Knowledge)
    push!(EP.knowledge,K)
    return(K)
end