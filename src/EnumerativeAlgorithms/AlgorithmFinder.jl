function find_algorithm(EProp::EnumerativeProperty,EP::EnumerativeProblem)
    potential_algorithms = filter(EA->output_property(EA)==EProp,MAIN_ALGORITHMS)
    println("There is a total of ",length(potential_algorithms), " algorithm(s) which are coded in Pandora.jl to compute ",name(EProp),":")
    counter=0
    for p in potential_algorithms
        counter=counter+1
        println("      ",counter,") ",name(p))
    end
    if length(potential_algorithms)==0
        return(nothing)
    else
        return(potential_algorithms[1])
    end
end