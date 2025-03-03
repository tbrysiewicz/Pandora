
function specialize(F::System; P = nothing)
    if P==nothing   
        P = randn(Float64,length(parameters(F)))
    end
    return(System(evaluate(expressions(F), parameters(F)=>P)))
end