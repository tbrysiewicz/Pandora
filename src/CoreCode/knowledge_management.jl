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

function attribute_bucket(EP::EnumerativeProblem, EAttr::EnumerativeAttribute)
    is_enumerative_data(EAttr) ? data(EP) : properties(EP)
end

function knows(EP::EnumerativeProblem, EProp::EnumerativeAttribute; kwargs...)
    # TODO: must check that the knowledge node has input kwargs which agree with kwargs.
    candidate_knowledge = filter(k -> property(k) == EProp, attribute_bucket(EP, EProp))
    kwarg_agreement = filter(k -> knowledge_agrees_with_kwargs(k; kwargs...), candidate_knowledge)
    if length(kwarg_agreement) > 0
        return true
    else
        return false
    end
end

function find_algorithm(EProp::EnumerativeAttribute, EP::EnumerativeProblem)
    potential_algorithms = algorithms_which_return(EProp)
    @vprintln("Pandora.jl is automatically finding an algorithm to compute ", EProp, ". To specify an algorithm, call again with algorithm=>[nameofalgorithm]")
    @vprintln("There is a total of ", length(potential_algorithms), " algorithm(s) in Pandora.jl which compute(s) ", name(EProp), ":")
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
    EProp::Union{EnumerativeProperty{T}, EnumerativeData{T}};
    algorithm = nothing,
    kwargs...
) where T
    # If the algorithm is not given, find an algorithm known to Pandora which
    # computes the given enumerative property.
    if algorithm === nothing
        algorithm = find_algorithm(EProp, EP)
    end

    f = algorithm
    @assert(ALGORITHM_DATA[f].output_property == EProp)
    input_knowledge = [get_knowledge!(i, EP; kwargs...) for i in input_properties(f)]
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
    new_knowledge = KnowledgeNode{T}(EProp, o, input_knowledge, kwargs_to_pass, f)
    return know!(EP, new_knowledge)
end

"""
    compute(EProp::EnumerativeProperty{T}, EP::EnumerativeProblem; algorithm = nothing, recompute_depth = 0, kwargs...) where T

Compute the value of an enumerative property `EProp` for an enumerative problem `EP`.
If `recompute_depth` is 1, it will recompute the property even if it is already known,
using the known input knowledge for that algorithm. If `recompute_depth` is 2, it will
recompute the property and recompute the input knowledge using recompute depth 1.
In general, `recompute_depth` is the recursive limit of recomputing knowledge, but
user-given information is always returned.
"""
function compute(
    EProp::Union{EnumerativeProperty{T}, EnumerativeData{T}},
    EP::EnumerativeProblem;
    algorithm = nothing,
    recompute_depth = 0,
    kwargs...
) where T
    K = get_knowledge(EProp, EP; kwargs...)
    if K !== nothing
        if length(K.input_knowledge) == 0
            return value(K)
            # @vprintln("Returning value of ", EProp, " from knowledge: ", value(K), ".")
        end
    end

    # println("Recompute depth is currently: ", recompute_depth, " on computation of ", EProp, ".")
    if recompute_depth == 0
        # print("Checking if the knowledge of ", EProp, " is already known in EP: ", EP, "...")
        g = get_knowledge_value(EProp, EP; kwargs...)
        if g !== nothing
            # print("Yes, it is known: ", g, "\n")
            return g
        end
        # println("No, it is not known. Computing it now...")
    end

    recompute_depth = max(recompute_depth - 1, 0)
    # If the algorithm is not given, find an algorithm known to Pandora which
    # computes the given enumerative property.
    if algorithm === nothing
        algorithm = find_algorithm(EProp, EP)
    end

    f = algorithm
    @assert(ALGORITHM_DATA[f].output_property == EProp)
    input_values = [compute(i, EP; recompute_depth = recompute_depth) for i in input_properties(f)]
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
    if is_enumerative_data(property(K))
        push!(EP.data, K)
    else
        push!(EP.properties, K)
    end
    return K
end

function know!(EP::EnumerativeProblem, EProp::Union{EnumerativeProperty{T}, EnumerativeData{T}}, value::T) where T
    K = KnowledgeNode(EProp, value, Vector{KnowledgeNode}(), Dict{Symbol, Any}(), user_given)
    know!(EP, K)
end

function get_knowledge(EProp::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    if !knows(EP, EProp; kwargs...)
        return nothing
    end

    candidate_knowledge = filter(k -> property(k) == EProp, attribute_bucket(EP, EProp))
    kwarg_agreement = filter(k -> knowledge_agrees_with_kwargs(k; kwargs...), candidate_knowledge)
    if length(kwarg_agreement) == 1
        return kwarg_agreement[1]
    else
        if is_enumerative_data(EProp)
            @vprintln("The data [", EProp, "] with given keyword arguments has more than one node. Returning the most recent one.")
            return kwarg_agreement[end]
        else
            @vprintln("Warning: the property [", EProp, "] with given keyword arguments has more than one knowledge node.")
            return kwarg_agreement[1]
        end
    end
end

function get_knowledge_value(EProp::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    K = get_knowledge(EProp, EP; kwargs...)
    if K !== nothing
        return value(K)
    else
        return nothing
    end
end

function get_knowledge!(EProp::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...)
    gk = get_knowledge(EProp, EP; kwargs...)
    if gk !== nothing
        return gk
    else
        learn!(EP, EProp; kwargs...)
        return get_knowledge!(EProp, EP; kwargs...)
    end
end

function remove_knowledge!(EP::EnumerativeProblem, EProp::EnumerativeAttribute)
    filter!(k -> property(k) != EProp, attribute_bucket(EP, EProp))
    return EP
end

get_knowledge_value!(EProp::EnumerativeAttribute, EP::EnumerativeProblem; kwargs...) = value(get_knowledge!(EProp, EP; kwargs...))

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

function algorithms_which_return(EProp::EnumerativeAttribute)
    filter(A -> output_property(ALGORITHM_DATA[A]) == EProp, collect(keys(ALGORITHM_DATA)))
end
