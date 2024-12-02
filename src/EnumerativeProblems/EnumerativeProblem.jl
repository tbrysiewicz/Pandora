################################################################
################################################################
############## Constructor and Show Function####################
################################################################
################################################################

mutable struct EnumerativeProblem <: AbstractEnumerativeProblem
    system::System       #EPs consist of a system and
    base_fibre::Fibre    #a fibre of solutions  and parameters

    function EnumerativeProblem(F::System)
        EP = new()
        EP.system=F
        EP.Data = Dict{Symbol,Any}()
        EP.Reasoning = Dict{Symbol,Symbol}()
        EP.Interpreter = Dict{Symbol,Any}()
        EP.HomotopyContinuationOptions = Dict{Any,Any}()
        EP.HomotopyContinuationOptions[:tracker_options]=TrackerOptions()
        return(EP)
    end

    function EnumerativeProblem(F::System,fibre::Fibre)
        EP = new()
        EP.system = F
        EP.base_fibre = fibre
        EP.Data = Dict{Symbol,Any}()
        EP.Reasoning = Dict{Symbol,Symbol}()
        EP.Interpreter = Dict{Symbol,Any}()
        EP.HomotopyContinuationOptions = Dict{Any,Any}()
        EP.HomotopyContinuationOptions[:tracker_options]=TrackerOptions()
        return(EP)
    end

    Data::Dict{Symbol,Any}          #Cache information about the problem
    Reasoning::Dict{Symbol,Symbol}  #Justify the key value pairs in Data via algorithm names
    Interpreter::Dict{Symbol,Any}   #Customize how to interpret solutions
    HomotopyContinuationOptions::Dict{Any,Any}

end



function Base.show(io::IO, E::EnumerativeProblem)
    tenspaces="          "
    print(io,"\n\n")
    print(io,tenspaces," X := V(")
    if n_polynomials(E)==1
        print(io,"f_1")
    else
        print(io,"f_1..f_",n_polynomials(E),"")
    end
    print(io,") ⊆ C^",ambient_dimension(E)," x C^",n_parameters(E),"\n");
    println(io,tenspaces," |")
    println(io,tenspaces," |")
    print(io,tenspaces," | π ")
    if haskey(data(E),:degree)
        println(io,"  ",degree(E),"-to-1")
    else
        println(io,"   ???-to-1")
    end
    println(io,tenspaces," |")
    println(io,tenspaces," V")
    println(io,tenspaces,"C^",n_parameters(E),"\n")
    println(io,"An enumerative problem in ",ambient_dimension(E)," variable(s) cut out by ", 
                n_polynomials(E)," condition(s) over ", n_parameters(E)," parameter(s).")
    println(io,"The following information is known about this problem:")
    for k in keys(data(E))
        println(io,"-",k)
    end
end


