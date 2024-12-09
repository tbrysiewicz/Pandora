function solutions(fibre::Fibre) :: Vector{Vector{ComplexF64}}
    return(fibre[1])
end

function parameters(fibre::Fibre) :: Vector{ComplexF64}
    return(fibre[2])
end


function accepts_signature_and_return_type(f, ::Type{T}) where T
    methods_matching = filter(method -> method.sig.types[2] == T, methods(f))
    if isempty(methods_matching)
        return false, nothing  
    else
        return_type = Core.Compiler.return_type(f, Tuple{T})
        return true, return_type
    end
end

#Any output of the function should be comparable
function is_score(f::Function)
    (acceptable, output_sig) = accepts_signature_and_return_type(f,Fibre)
    if acceptable
        if hasmethod(isless,Tuple{output_sig,output_sig})
            return(true)
        else
            println("Cannot compare outputs of given function")
            return(false)
        end
    end
    (racceptable, routput_sig) = accepts_signature_and_return_type(f,RealFibre)
    if racceptable && hasmethod(isless,Tuple{routput_sig,routput_sig})
        println("Currently, score input must be fibres in terms of complex solution sets and complex parameters.")
        return(false)
    end

    (Csol_acceptable,Csol_output_sig) = accepts_signature_and_return_type(f,Vector{Vector{ComplexF64}})
    (Rsol_acceptable,Rsol_output_sig) = accepts_signature_and_return_type(f,Vector{Vector{Float64}})
    if Csol_acceptable==true
        println("Your function f takes solutions as input. Turn it into a score function using one of the following \n -max_score(f) \n -min_score(f) \n -sum_score(f) \n -etc")
        return(false)
    end
    if Rsol_acceptable==true
        println("Your function (1) takes solutions (but not fibres) as input, and (2) does so only with real solutions.")
    end
    println("Your function, as written, is not recognized as a score function.  
    Consider writing it as to enforce the input type.")
    return(false)
end
