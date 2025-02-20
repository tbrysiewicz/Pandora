function get_knowledge!(EProp::EnumerativeProperty,EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    gk = get_knowledge(EProp,EP;user_parameters=user_parameters)
    if gk != nothing #Check if knowledge of this property is already known
        return(gk)   #If so, return it
    else             #Otherwise,
        algorithm = find_algorithm(EProp,EP)  #find and algorithm which will obtain it
        if algorithm != nothing
                                              #and learn that knowledge
            return(learn!(EP,EProp;algorithm = algorithm, user_parameters=user_parameters))
        end
    end
end


function get_knowledge_value!(EProp::EnumerativeProperty,EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    K = get_knowledge!(EProp,EP; user_parameters=user_parameters)
    return(value(K))
end


function get_knowledge(EProp::EnumerativeProperty,EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    if knows(EP,EProp)==false
        return(nothing)
    end
    K = knowledge(EP)
    prop_indices = findall(p->p.property==EProp,K)
    if  length(prop_indices)==1
        return(K[prop_indices[1]])
    else
        println("Warning: the property [",EProp,"] has more than one knowledge node.")
        return(K[prop_indices[1]])
    end
end

function get_knowledge_value(EProp::EnumerativeProperty,EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    K = get_knowledge(EProp,EP; user_parameters=user_parameters)
    if K != nothing
        return(value(K))
    else
        return(nothing)
    end
end

function (EProp::EnumerativeProperty)(EP::EnumerativeProblem; user_parameters :: Dict{Symbol,Any} = Dict{Symbol,Any}())
    get_knowledge_value!(EProp,EP; user_parameters=user_parameters)
end


function degree(EP::EnumerativeProblem)
    enumerative_degree(EP)
end
