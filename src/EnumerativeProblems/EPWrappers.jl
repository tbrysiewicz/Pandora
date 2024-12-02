

function degree(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute
        delete!(data(EP),:degree)
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