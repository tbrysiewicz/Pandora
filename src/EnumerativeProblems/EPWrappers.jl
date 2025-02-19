
#A wrapper is something that directs how something is computed based on what it 
#   already knows. e.g. if we've computed something before, a wrapper will navigate toward
#   looking it up instead of recomputation. If a computation of something is specialized
#   for specific objects, the wrapper will assess the best algorithm for computing the 
#   desired feature of the object. 

#Think of wrappers for the things that the user touches: they wrap together the stuff the
#  user shouldn't see/have to worry about.


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
    get!(data(EP),:monodromy_group) do
        justify(EP,:monodromy_group,"Computed from monodromy dictionary.")
        compute_monodromy_group(EP)
    end
end

function monodromy_dictionary(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:monodromy_dictionary) 
        println(RecomputationWarning)
    end
    get!(data(EP),:monodromy_dictionary) do
         justify(EP,:monodromy_dictionary,"Computed numerically.")
         compute_monodromy_dictionary(EP)
    end
end


function galois_group(EP::EnumerativeProblem; force_recompute=false)
    monodromy_group(EP; force_recompute = force_recompute)
end
 
function transitivity_basis(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:transitivity_basis) 
        println(RecomputationWarning)
    end
    get!(data(EP),:transitivity_basis) do
        transitivity_basis(galois_group(EP))
    end
end   

function components(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:components) 
        println(RecomputationWarning)
    end
    get!(data(EP),:components) do
        compute_components(EP)
    end
end




function bezout(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:bezout) 
        println(RecomputationWarning)
    end
    get!(data(EP),:bezout) do
        compute_bezout(EP)
    end
end


function bkk(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:bkk) 
        println(RecomputationWarning)
    end
    get!(data(EP),:bkk) do
        compute_bkk(EP)
    end
end


function affine_bkk(EP::EnumerativeProblem; force_recompute=false)
    if force_recompute 
        delete!(data(EP),:affine_bkk) 
        println(RecomputationWarning)
    end
    get!(data(EP),:affine_bkk) do
        compute_affine_bkk(EP)
    end
end