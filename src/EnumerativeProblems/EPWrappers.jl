
#A wrapper is something that directs how something is computed based on what it 
#   already knows. e.g. if we've computed something before, a wrapper will navigate toward
#   looking it up instead of recomputation. If a computation of something is specialized
#   for specific objects, the wrapper will assess the best algorithm for computing the 
#   desired feature of the object. 


function degree(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute
        delete!(data(EP),:degree)
        println(AlterWarning)
        println(RecomputationWarning)
        populate!(EP)
    end
    if haskey(data(EP),:degree)
        return(data(EP)[:degree])
    else
        if !is_populated(EP) 
            println(AlterWarning)
            populate!(EP)
        end
        d = length(base_solutions(EP))
        data(EP)[:degree]=d
        return(d)
    end
end



function monodromy_group(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:monodromy_group) 
        println(RecomputationWarning)
    end
    get!(data(EP),:monodromy_group,compute_monodromy_group(EP))
end

function monodromy_dictionary(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:monodromy_group) 
        println(RecomputationWarning)
    end
    get!(data(EP),:monodromy_dictionary,compute_monodromy_dictionary(EP))
end


function galois_group(EP::EnumerativeProblem; force_recompute=false)
    monodromy_group(EP; force_recompute = force_recompute)
end


function components(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:components) 
        println(RecomputationWarning)
    end
    get!(data(EP),:components,compute_components(EP))
end