################################################################
################################################################
##############       Fundamental Populators      ###############
################################################################
################################################################

#If one changes their mind about what to call EnumerativeProblem.base_fibre, or EnumerativeProblem.system
#   one should only need to change the fundamental getters and fundamental populators (see EPPopulators.jl)
#   for enumerative problems. Everything esle should work fine. 


function update_base_fibre!(EP::EnumerativeProblem,fibre::Fibre)
    EP.base_fibre = fibre
    return(EP)
end

################################################################
################################################################
##############       Simple  Populators          ###############
################################################################
################################################################


function populate!(EP::EnumerativeProblem)
    new_param = randn(ComplexF64,n_parameters(EP))
    S = EP(new_param)
    update_base_fibre!(EP,(S,new_param))
    return(EP)
end


function solve!(EP::EnumerativeProblem)
    return(populate!(EP))
end


function monodromy_solve!(EP::EnumerativeProblem)
    println("Warning: monodromy solve finds a fibre within only one component")
    F = system(EP)
    MS = monodromy_solve(F)
    update_base_fibre!(EP,(solutions(MS),parameters(MS)))
    return(EP)
end
