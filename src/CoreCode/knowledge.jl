export
    KnowledgeNode,
    attribute,
    value,
    input_knowledge,
    input_kwargs,
    algorithm

"""
`KnowledgeNode` stores one piece of knowledge for an enumerative problem.

It stores the attribute, its value, the input knowledge used to compute it, the
keyword arguments used, and the algorithm that produced it.
"""
mutable struct KnowledgeNode{T}
    attribute::Union{EnumerativeProperty{T}, EnumerativeData{T}}
    value::T
    input_knowledge::Vector{KnowledgeNode}
    input_kwargs::Dict{Symbol, Any}
    algorithm::Function
end

const Knowledge = Vector{KnowledgeNode}

attribute(K::KnowledgeNode) = K.attribute
value(K::KnowledgeNode) = K.value
input_knowledge(K::KnowledgeNode) = K.input_knowledge
input_kwargs(K::KnowledgeNode) = K.input_kwargs
algorithm(K::KnowledgeNode) = K.algorithm

known_properties(K::Knowledge) = unique([attribute(k) for k in K])
known_data(K::Knowledge) = unique([attribute(k) for k in K])

function Base.show(io::IO, K::Knowledge)
    for (i, k) in enumerate(K)
        print(io, i, ") ", string(k), "\n")
    end
end

function Base.show(io::IO, ::MIME"text/plain", K::Knowledge)
    for (i, k) in enumerate(K)
        print(io, i, ") ", string(k), "\n")
    end
end

function Base.string(K::KnowledgeNode{T}) where T
    s = "[" * string(attribute(K))
    s *= "] as computed by (" * string(name(algorithm(K))) * ") applied to "
    if length(input_knowledge(K)) == 0
        s *= "(nothing)."
    else
        s *= "[" * string(attribute(input_knowledge(K)[1]))
        for i in input_knowledge(K)[2:end]
            s *= ", " * string(attribute(i))
        end
        s *= "]"
    end
    return s
end

function Base.show(io::IO, K::KnowledgeNode)
    show_knowledge_tree(io, K; isroot = true)
end

function Base.show(io::IO, ::MIME"text/plain", K::KnowledgeNode)
    show_knowledge_tree(io, K; isroot = true)
end

function Base.display(K::KnowledgeNode)
    show_knowledge_tree(stdout, K; isroot = true)
end

function show_knowledge_tree(io::IO, K::KnowledgeNode; prefix = "", islast = true, isroot = false)
    connector = islast ? "└──" : "├──"
    line = prefix * connector * " [" * string(attribute(K)) * "]"
    if isroot
        line *= " as computed by (" * string(name(algorithm(K))) * ")"
    end
    println(io, line)

    new_prefix = prefix * (islast ? "    " : "│   ")
    children = input_knowledge(K)
    n = length(children)
    for (i, child) in enumerate(children)
        show_knowledge_tree(io, child; prefix = new_prefix, islast = i == n)
    end
end

function knowledge_tree(K::KnowledgeNode; prefix = "", islast = true)
    show_knowledge_tree(stdout, K; prefix = prefix, islast = islast)
end

function is_certified(K::KnowledgeNode)
    return reliability(K) == [:user_given]
end

function reliability_consensus(reliability_bucket::Vector{Symbol})
    reliability_bucket = unique(filter(r -> r != :certified, reliability_bucket))
end

function reliability(K::KnowledgeNode)
    reliability_bucket = [ALGORITHM_DATA[algorithm(K)].reliability]
    for i in input_knowledge(K)
        for r in reliability(i)
            push!(reliability_bucket, r)
        end
    end
    return reliability_consensus(reliability_bucket)
end

function combined_input_knowledge(K::Knowledge)
    vcat([input_knowledge(k) for k in K]...)
end

function check_same_attribute(K::Knowledge)
    isempty(K) && error("Cannot combine an empty knowledge list")
    EA = attribute(first(K))
    all(k -> attribute(k) == EA, K) || error("Cannot combine knowledge nodes with different attributes")
    return EA
end

function combine_properties(K::Knowledge, EA::EnumerativeProperty)
    first_value = value(first(K))
    all(k -> value(k) == first_value, K) || error("Cannot combine property knowledge nodes with conflicting values")

    T = get_type(EA)
    return KnowledgeNode{T}(EA, first_value, combined_input_knowledge(K), Dict{Symbol, Any}(), conjunction)
end

function combine_data(K::Knowledge, EA::EnumerativeData)
    T = get_type(EA)
    combined_value = reduce(vcat, [value(k) for k in K])
    combined_value isa T || error("Cannot combine data values into the declared attribute type")

    return KnowledgeNode{T}(EA, combined_value, combined_input_knowledge(K), Dict{Symbol, Any}(), conjunction)
end

function combine_knowledge(K::Knowledge)
    EA = check_same_attribute(K)
    return is_enumerative_data(EA) ? combine_data(K, EA) : combine_properties(K, EA)
end
