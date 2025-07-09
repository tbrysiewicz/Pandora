
function monodromy_summary(EP::EnumerativeProblem)
    summary = raw"\section{Monodromy Data}"
    
    if knows(EP, MONODROMY_GROUP)
        summary *= raw"The monodromy group of the enumerative problem "
        if is_natural_symmetric_group(monodromy_group(EP))
            summary *= raw"is the full symmetric group of degree "*string(degree(EP))*"."
        elseif is_natural_alternating_group(monodromy_group(EP))
            summary *= raw"is the full alternating group of degree "*string(degree(EP))*"."
        else
            summary *= raw"is a subgroup of the symmetric group of degree "*string(degree(EP))*"."
        
            summary *=" given by the following generators: "
            summary *= raw"\begin{align*}"*"\n "
            for (i, gen) in enumerate(gens(monodromy_group(EP)))
                summary *= raw"""g_{"""*string(i)*raw"""} &= """*string(gen)*raw"""\\\\"""*"\n "
            end
            summary *= raw"\end{align*}"*"\n "
            summary *= raw"The order of the monodromy group is $"*string(order(monodromy_group(EP)))*
                    raw"$, and it is "
            if is_transitive(monodromy_group(EP))
                summary *= raw"transitive."
                if is_primitive(monodromy_group(EP))
                    summary *= raw" It is also primitive."
                else
                    summary *= raw" It is not primitive."

                    summary *= raw"The minimal block representatives of the monodromy group are given by the following partitions: "
                    summary *= raw"\begin{align*}"*"\n "
                    for (i, part) in enumerate(minimal_block_reps(monodromy_group(EP)))
                        summary *= raw"""B_{"""*string(i)*raw"""} &= """*string(part)*"\n "
                    end
                    summary *= raw"\end{align*}"*"\n "
                end
            else
                summary *= raw"intransitive."
            end
        end
    else
        summary *= raw"The monodromy group of the enumerative problem is not known."
    end  
    return(summary)
end