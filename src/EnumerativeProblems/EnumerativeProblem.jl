################################################################
################################################################
############## Constructor and Show Function####################
################################################################
################################################################

mutable struct EnumerativeProblem <: AbstractEnumerativeProblem
    system :: System
    knowledge :: Vector{Knowledge}
    hc_options :: Dict{Any,Any}

    function EnumerativeProblem(F::System; populate = true)
        EP = new()
        EP.system = F
        EP.knowledge = Vector{Knowledge}()
        EP.hc_options = Dict{Any,Any}()
        EP.hc_options[:tracker_options]=TrackerOptions()

        know!(EP,system,F)

        if populate
            learn!(EP,base_fibre; algorithm = polyhedral)
            learn!(EP,enumerative_degree; algorithm = degree_from_base_fibre)
        end
        return(EP)
    end
end




function Base.show(io::IO, EP::EnumerativeProblem)
    tenspaces="          "
    print(io,"\n\n")
    print(io,tenspaces," X := V(")
    if n_polynomials(EP)==1
        print(io,"f_1")
    else
        print(io,"f_1..f_",n_polynomials(EP),"")
    end
    print(io,") ⊆ C^",ambient_dimension(EP)," x C^",n_parameters(EP),"\n");
    println(io,tenspaces," |")
    println(io,tenspaces," |")
    print(io,tenspaces," | π ")
    if knows(EP,enumerative_degree)
        println(io,"  ",degree(EP),"-to-1")
    else
        println(io,"   ???-to-1")
    end
    println(io,tenspaces," |")
    println(io,tenspaces," V")
    println(io,tenspaces,"C^",n_parameters(EP),"\n")
    println(io,"An enumerative problem in ",ambient_dimension(EP)," variable(s) cut out by ", 
                n_polynomials(EP)," condition(s) over ", n_parameters(EP)," parameter(s).")
    println(io,"The following information is known about this problem:")
    for k in knowledge(EP)
        println(io,"-",property(k))
    end
end