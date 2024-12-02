

function degree(EP::EnumerativeProblem)
    if !is_populated(EP)
        println(AlterWarning)
        populate!(EP)
    end
    return(length(base_solutions(EP)))
end



function monodromy_group(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:monodromy_group) 
        println(RecomputationWarning)
    end
    get!(data(EP),:monodromy_group,compute_monodromy_group(EP))
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