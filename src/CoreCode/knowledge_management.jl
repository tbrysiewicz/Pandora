function knowledge_agrees_with_kwargs(k::KnowledgeNode; kwargs...)
    isempty(kwargs) && return true
    ikwargs = input_kwargs(k)
    isempty(ikwargs) && return true

    d = Dict(kwargs)
    for i in keys(ikwargs)
        if haskey(d, i)
            if d[i] != ikwargs[i]
                return false
            end
        end
    end
    return true
end

function attribute_bucket(EP::EnumerativeProblem, EA::EnumerativeAttribute)
    is_enumerative_data(EA) ? data(EP) : properties(EP)
end

function knows(EP::EnumerativeProblem, EA::EnumerativeAttribute; kwargs...)
    # TODO: must check that the knowledge node has input kwargs which agree with kwargs.
    candidate_knowledge = filter(k -> attribute(k) == EA, attribute_bucket(EP, EA))
    kwarg_agreement = filter(k -> knowledge_agrees_with_kwargs(k; kwargs...), candidate_knowledge)
    if length(kwarg_agreement) > 0
        return true
    else
        return false
    end
end

function find_algorithm(EA::EnumerativeAttribute, EP::EnumerativeProblem)
    potential_algorithms = algorithms_which_return(EA)
    @vprintln("Pandora.jl is automatically finding an algorithm to compute ", EA, ". To specify an algorithm, call again with algorithm=>[nameofalgorithm]")
    @vprintln("There is a total of ", length(potential_algorithms), " algorithm(s) in Pandora.jl which compute(s) ", name(EA), ":")
    counter = 0
    length(potential_algorithms) == 0 && return nothing
    algorithm_to_use = potential_algorithms[1]
    for p in potential_algorithms
        counter += 1
        @vprint("      ", counter, ") ")
        if p == algorithm_to_use
            @vprintln("[USING] ", name(p))
        else
            @vprintln(name(p))
        end
    end
    return algorithm_to_use
end

function learn!(
    EP::EnumerativeProblem,
    EA::Union{EnumerativeProperty{T}, EnumerativeData{T}};
    algorithm = nothing,
    kwargs...
) where T
    # If the algorithm is not given, find an algorithm known to Pandora which
    # computes the given enumerative attribute.
    if algorithm === nothing
        algorithm = find_algorithm(EA, EP)
    end

    f = algorithm
    @assert(ALGORITHM_DATA[f].output_attribute == EA)
    input_knowledge = [get_knowledge!(i, EP; kwargs...) for i in input_attributes(f)]
    input_knowledge_values = [value(kn) for kn in input_knowledge]
    kwargs_to_pass = copy(default_kwargs(f))
    if !isempty(kwargs)
        for kv in kwargs
            if haskey(kwargs_to_pass, kv[1])
                kwargs_to_pass[kv[1]] = kv[2]
            end
        end
    end
    o = f(input_knowledge_values...; kwargs_to_pass...)
    new_knowledge = KnowledgeNode{T}(EA, o, input_knowledge, kwargs_to_pass, f)
    return know!(EP, new_knowledge)
end

"""
    compute(EA::EnumerativeAttribute, EP::EnumerativeProblem; algorithm = nothing, recompute_depth = 0, kwargs...)

Compute the value of an enumerative attribute `EA` for an enumerative problem `EP`.
If `recompute_depth` is 1, it will recompute the attribute even if it is already known,
using the known input knowledge for that algorithm. If `recompute_depth` is 2, it will
recompute the attribute and recompute the input knowledge using recompute depth 1.
In general, `recompute_depth` is the recursive limit of recomputing knowledge, but
user-given information is always returned.
"""
function compute(
    EA::Union{EnumerativeProperty{T}, EnumerativeData{T}},
    EP::EnumerativeProblem;
    algorithm = nothing,
    recompute_depth = 0,
    kwargs...
) where T
    K = get_knowledge(EA, EP; kwargs...)
    if K !== nothing
        if length(K.input_knowledge) == 0
            return value(K)
            # @vprintln("Returning value of ", EA, " from knowledge: ", value(K), ".")
        end
    end

    # println("Recompute depth is currently: ", recompute_depth, " on computation of ", EA, ".")
    if recompute_depth == 0
        # print("Checking if the knowledge of ", EA, " is already known in EP: ", EP, "...")
        g = get_knowledge_value(EA, EP; kwargs...)
        if g !== nothing
            # print("Yes, it is known: ", g, "\n")
            return g
        end
        # println("No, it is not known. Computing it now...")
    end

    recompute_depth = max(recompute_depth - 1, 0)
    # If the algorithm is not given, find an algorithm known to Pandora which
    # computes the given enumerative attribute.
    if algorithm === nothing
        algorithm = find_algorithm(EA, EP)
    end

    f = algorithm
    @assert(ALGORITHM_DATA[f].output_attribute == EA)
    input_values = [compute(i, EP; recompute_depth = recompute_depth) for i in input_attributes(f)]
    kwargs_to_pass = copy(default_kwargs(f))
    if !isempty(kwargs)
        for kv in kwargs
            if haskey(kwargs_to_pass, kv[1])
                kwargs_to_pass[kv[1]] = kv[2]
            end
        end
    end
    o = f(input_values...; kwargs_to_pass...)
    return o
end

function know!(EP::EnumerativeProblem, K::KnowledgeNode)
    if is_enumerative_data(attribute(K))
        push!(EP.data, K)
    else
        push!(EP.properties, K)
    end
    return K
end

function know!(EP::EnumerativeProblem, EA::Union{EnumerativeProperty{T}, EnumerativeData{T}}, value::T) where T
    K = KnowledgeNode(EA, value, Vector{KnowledgeNode}(), Dict{Symbol, Any}(), user_given)
    know!(EP, K)
end

function get_knowledge(EA::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    if !knows(EP, EA; kwargs...)
        return nothing
    end

    candidate_knowledge = filter(k -> attribute(k) == EA, attribute_bucket(EP, EA))
    kwarg_agreement = filter(k -> knowledge_agrees_with_kwargs(k; kwargs...), candidate_knowledge)
    if length(kwarg_agreement) == 1
        return kwarg_agreement[1]
    else
        if is_enumerative_data(EA)
            @vprintln("The data [", EA, "] with given keyword arguments has more than one node. Returning the most recent one.")
            return kwarg_agreement[end]
        else
            @vprintln("Warning: the property [", EA, "] with given keyword arguments has more than one knowledge node.")
            return kwarg_agreement[1]
        end
    end
end

function get_knowledge_value(EA::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    K = get_knowledge(EA, EP; kwargs...)
    if K !== nothing
        return value(K)
    else
        return nothing
    end
end

function get_knowledge!(EA::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EA, EP; kwargs...)
    if gk !== nothing
        return gk
    else
        learn!(EP, EA; kwargs...)
        return get_knowledge!(EA, EP; kwargs...)
    end
end

function remove_knowledge!(EP::EnumerativeProblem, EA::EnumerativeAttribute)
    filter!(k -> attribute(k) != EA, attribute_bucket(EP, EA))
    return EP
end

get_knowledge_value!(EA::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...) = value(get_knowledge!(EA, EP; kwargs...))

function (EProp::EnumerativeProperty)(EP::EnumerativeProblem; learn = true, recompute_depth = 1, kwargs...)
    if learn == true
        return get_knowledge_value!(EProp, EP; kwargs...)
    else
        return compute(EProp, EP; recompute_depth = recompute_depth, kwargs...)
    end
end

function (EData::EnumerativeData)(EP::EnumerativeProblem; learn = true, recompute_depth = 1, kwargs...)
    if learn == true
        return get_knowledge_value!(EData, EP; kwargs...)
    else
        return compute(EData, EP; recompute_depth = recompute_depth, kwargs...)
    end
end

function algorithms_which_return(EA::EnumerativeAttribute)
    filter(A -> output_attribute(ALGORITHM_DATA[A]) == EA, collect(keys(ALGORITHM_DATA)))
end
