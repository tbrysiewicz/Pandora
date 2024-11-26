################################################################
################################################################
##############       Fundamental Getters         ###############
################################################################
################################################################


#If one changes their mind about what to call EnumerativeProblem.base_fibre, or EnumerativeProblem.system
#   one should only need to change the fundamental getters and fundamental populators (see EPPopulators.jl)
#   for enumerative problems. Everything esle should work fine. 

function base_fibre(EP::EnumerativeProblem) :: Fibre
    EP.base_fibre
end

function system(EP::EnumerativeProblem) :: System
    EP.system
end

function is_populated(EP::EnumerativeProblem) :: Bool
    isdefined(EP,:base_fibre)
end


################################################################
################################################################
##############          Simple Getters           ###############
################################################################
################################################################



function base_solutions(EP::EnumerativeProblem)
    is_populated(EP) && return(solutions(base_fibre(EP)))
    return(nothing)
end

function solutions(EP::EnumerativeProblem)
    return(base_solutions(EP))
end

function base_parameters(EP::EnumerativeProblem)
    is_populated(EP) && return(parameters(base_fibre(EP)))
    return(nothing)
end


function ambient_dimension(EP::EnumerativeProblem)
    length(variables(system(EP)))
end

function n_polynomials(EP::EnumerativeProblem)
    length(expressions(system(EP)))
end

function n_parameters(EP::EnumerativeProblem)
    length(parameters(system(EP)))
end

function variables(EP::EnumerativeProblem)
    variables(system(EP))
end

function parameters(EP::EnumerativeProblem)
    parameters(system(EP))
end


function degree(EP::EnumerativeProblem)
    if !is_populated(EP)
        println("Warning: this function will alter the inputted enumerative problem by computing a base fibre")
        populate!(EP)
    end
    return(length(base_solutions(EP)))
end

