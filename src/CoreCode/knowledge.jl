export
    KnowledgeNode,
    property,
    value,
    input_knowledge,
    input_kwargs,
    algorithm

"""
`KnowledgeNode` records a property or data value for an enumerative problem.

It stores the attribute, its value, the input knowledge used to compute it, the
keyword arguments used, and the algorithm that produced it.
"""
mutable struct KnowledgeNode{T}
    property::Union{EnumerativeProperty{T}, EnumerativeData{T}}
    value::T
    input_knowledge::Vector{KnowledgeNode}
    input_kwargs::Dict{Symbol, Any}
    algorithm::Function
end

const Knowledge = Vector{KnowledgeNode}

property(K::KnowledgeNode) = K.property
value(K::KnowledgeNode) = K.value
input_knowledge(K::KnowledgeNode) = K.input_knowledge
input_kwargs(K::KnowledgeNode) = K.input_kwargs
algorithm(K::KnowledgeNode) = K.algorithm

known_properties(K::Knowledge) = unique([property(k) for k in K])
known_data(K::Knowledge) = unique([property(k) for k in K])

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
    s = "[" * string(property(K))
    s *= "] as computed by (" * string(name(algorithm(K))) * ") applied to "
    if length(input_knowledge(K)) == 0
        s *= "(nothing)."
    else
        s *= "[" * string(property(input_knowledge(K)[1]))
        for i in input_knowledge(K)[2:end]
            s *= ", " * string(property(i))
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
    line = prefix * connector * " [" * string(property(K)) * "]"
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

function combine_knowledge(K::Knowledge)
    new_type = type((map(k -> get_type(property(k)), K)))
    new_name = join([name(property(k)) for k in K], " & ")
    combined_EP = EnumerativeProperty{new_type}(new_name)
    combined_value = (map(k -> value(k), K))
    combined_input_knowledge = vcat([input_knowledge(k) for k in K]...)
    combined_input_kwargs = Dict{Symbol, Any}()
    new_knowledge = KnowledgeNode{new_type}(combined_EP, combined_value, combined_input_knowledge, combined_input_kwargs, conjunction)
    return new_knowledge
end
