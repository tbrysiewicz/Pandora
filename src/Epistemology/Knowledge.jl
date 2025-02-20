mutable struct Knowledge{T}
    property :: EnumerativeProperty{T}
    value :: T
    input_knowledge :: Vector{Knowledge}
    input_parameters :: Dict{Symbol,Any}
    algorithm :: EnumerativeAlgorithm
end

function property(K::Knowledge)
    K.property
end
function value(K::Knowledge)
    K.value
end
function input_knowledge(K::Knowledge)
    K.input_knowledge
end
function input_parameters(K::Knowledge)
    K.input_parameters
end
function algorithm(K::Knowledge)
    K.algorithm
end

function Base.show(io::IO, K::Knowledge)
    p = property(K)
    println(io,"Knowledge node for [",property(K),"]:")
    println(io,value(K))
    println(io,"----------------------------")
    println(io,"Algorithm used: [",name(algorithm(K)),"]")
    if length(input_knowledge(K))>0
        println(io,"Ran on knowledge nodes for ")
    end
    for i in input_knowledge(K)
        println(io,"     ",name(property(i)))
    end
    if length(keys(input_parameters(K)))>0
        println(io," with input parameters")
        for ip in keys(input_parameters(K))
            println(io," ",ip,": ",(input_parameters(K))[ip])
        end
    end
end
