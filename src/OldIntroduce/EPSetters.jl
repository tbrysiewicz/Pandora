################################################################
################################################################
##############       Fundamental Setters         ###############
################################################################
################################################################

#If one changes their mind about what to call EnumerativeProblem.base_fibre, or EnumerativeProblem.system
#   one should only need to change the fundamental getters and fundamental populators (see EPPopulators.jl)
#   for enumerative problems. Everything esle should work fine. 

function update_base_fibre!(EP::EnumerativeProblem,fibre::Fibre)
    EP.base_fibre = fibre
    data(EP)[:degree]=length(solutions(fibre))
    return(EP)
end

################################################################
################################################################
##############       Simple  Populators          ###############
################################################################
################################################################





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
