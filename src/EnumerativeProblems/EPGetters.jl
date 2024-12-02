################################################################
################################################################
##############       Fundamental Getters         ###############
################################################################
################################################################


#If one changes their mind about what to call EnumerativeProblem.base_fibre, or EnumerativeProblem.system
#   one should only need to change the fundamental getters and fundamental populators (see EPPopulators.jl)
#   for enumerative problems. Everything esle should work fine. 

function base_fibre(EP::EnumerativeProblem) :: Fibre
    if is_populated(EP)
        EP.base_fibre
    else
        println("Base fibre has not been compute. Call populate!(::EnumerativeProblem) or solve!(::EnumerativeProblem)")
    end
end

function system(EP::EnumerativeProblem) :: System
    EP.system
end

function data(EP::EnumerativeProblem) :: Dict{Symbol,Any}
    EP.Data
end
function homotopy_continuation_options(EP::EnumerativeProblem) 
    EP.HomotopyContinuationOptions
end

function tracker_options(EP::EnumerativeProblem) :: TrackerOptions
    homotopy_continuation_options(EP)[:tracker_options]
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

