

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