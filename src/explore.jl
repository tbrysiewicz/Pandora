export
    explore,
    reality_exploration,
    reality_exploration!,
    fibre_datum,
    fibre_data,
    explore_reality,
    tally


function tally(S)
    # Count the number of occurrences of each unique solution
    counts = Dict{Any, Int}()
    for s in S
        if haskey(counts, s)
            counts[s] += 1
        else
            counts[s] = 1
        end
    end
    return counts
end



"""
    explore(EP::EnumerativeProblem, f::Function;  sampler = UniformSampler{Float64}(n_parameters(EP)), real_parameters=true, n_samples=1000, as_fibres = false, kwargs...)
    explore(EP::EnumerativeProblem, F::AbstractVector{<:Function}; sampler = UniformSampler{Float64}(n_parameters(EP)), real_parameters=true, n_samples=1000, as_fibres = false, kwargs...)

Explores the parameter space of an EnumerativeProblem `EP` by sampling parameters, solving `EP` over each parameter, and applying `f` or each function in `F`.
It returns a vector of function values for each function in `F` or a vector of `FibreDatum` objects if `as_fibres` is set to `true`. 
The function accepts the following keyword arguments:
- `n_samples`: The number of samples to draw from the parameter space (default is 1000).
- `sampler`: A sampler function to generate samples from the parameter space (default is `UniformSampler{Float64}(n_parameters(EP))`).
- `real_parameters`: If `true`, the sampler will generate real parameters; if `false`, it will generate complex parameters.
- `as_fibres`: If `true`, the function returns a vector of `FibreDatum` objects, each containing the function values for the fibre solutions; if `false`, it returns a vector of vectors of function values for each function in `F`.

Note: The function filters out fibres that are not generic according to the `valid_fibre` function. 
"""
function explore(EP::EnumerativeProblem, f::Function;  sampler = UniformSampler{Float64}(n_parameters(EP)), real_parameters=true, n_samples=1000,  as_fibres = false,kwargs...)::Vector{FibreDatum}
    explore(EP, [f]; n_samples = n_samples, sampler = sampler, real_parameters=real_parameters,  as_fibres = as_fibres, kwargs...) 
end

function explore(EP::EnumerativeProblem, F::AbstractVector{<:Function}; sampler = UniformSampler{Float64}(n_parameters(EP)), real_parameters=true, n_samples=1000,  as_fibres = false, kwargs...)
    #set the sampler up
    if sampler == nothing
        if real_parameters
            sampler = UniformSampler{Float64}(n_parameters(EP))
        else
            sampler = UniformSampler{ComplexF64}(n_parameters(EP))
        end
    end

    P = sampler(n_samples)
    fibres = map(Fibre,collect(zip(EP(P),P)))
    filter!(fibre -> valid_fibre(EP, fibre), fibres)
    @vprintln("Number of valid fibres:$(length(fibres))")
    if as_fibres == false
        return([[f(fib[1]) for fib in fibres] for f in F])
    end

    fibres = map(FibreDatum, fibres)
    for fib in fibres
        for f in F 
            fib.function_values[f] = f(fib.solutions)
        end
    end
    return(fibres)
end


function tally_and_witness_fibres(Fibres,f; certify_function = nothing)
    # Count the number of occurrences of each unique solution and return a witness
    counts = Dict{Any, Int}()
    witnesses = Dict{Any, Any}()
    for i in eachindex(Fibres)
        s = Fibres[i][1]
        v = f(s)
        if haskey(counts, v)
            counts[v] += 1
            if !haskey(witnesses, v)
                @vprint("Trying again to certify that function value $v is found in fibre number $(i)")
                if certify_function(Fibres[i],v)
                    witnesses[v] = Fibres[i][2]  # Store the first occurrence as a witness
                else

                end
            end
        else
            counts[v] = 1
            if certify_function !== nothing 
                # If a certify function is provided, use it to certify the solution
                @vprint("Certifying that function value $v is found in fibre number $(i)")
                if certify_function(Fibres[i],v)
                    witnesses[v] = Fibres[i][2]  # Store the first occurrence as a witness
                else

                end
            else
                witnesses[v] = Fibres[i][2]  # Store the first occurrence as a witness
            end
        end
    end
    return counts, witnesses
end



function explore_reality(EP::EnumerativeProblem; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), certify = true, as_fibres = false, kwargs...)
    certification_function = (fibre,n) -> true  # Default certification function that always returns true
    if certify
        certification_function = (fibre,n) -> certify_n_real(EP, fibre) == n
    end

    F = explore(EP, n_real_solutions; sampler = sampler, n_samples = n_samples, real_parameters = true, as_fibres = true)
    real_values = sort(unique([f.function_values[n_real_solutions] for f in F]))
    @vprintln("Number of real solutions found: $real_values")
    certified_fibres = Vector{FibreDatum}()
    unknown_real_values = filter(u -> all(fb -> get(fb.function_values,n_real_solutions,nothing) != u, fibre_data(EP)), real_values)
    length(unknown_real_values)>0 && @vprintln("Certifying a witness for each real solution count not already known in the knowledge base: $unknown_real_values")
    for u in unknown_real_values
        fibres_with_u = filter(fibre -> fibre.function_values[n_real_solutions] == u, F)
        certified = false
        certified_fibre = nothing
        i=1
        while certified == false && i <= length(fibres_with_u)
            #@vprintln("Trying to certify that function value $u is found in instance $(i) of $(length(fibres_with_u))")
            fib = fibres_with_u[i]
            certify!(fib, EP)
            if n_real_solutions_certified(fib) == u
                #@vprintln("Certified that function value $u is found in fibre number $(i)")
                certified = true
                certified_fibre = fib
                push!(certified_fibres, fib)
            else
                #@vprintln("Failed to certify that function value $u is found in fibre number $(i)")
            end
            i += 1
        end
        if certified == false
            @vprintln("Failed to certify a witness for function value $u")
        else
            #@vprintln("Successfully certified a witness for function value $u")
            KN = KnowledgeNode(
                FIBRE_DATUM,
                certified_fibre,
                [get_knowledge(SYSTEM,EP), get_knowledge(BASE_FIBRE,EP)],
                Dict{Symbol, Any}(),
                explore_reality
                )
            know!(EP, KN)
        end
    end
    if as_fibres
        return F
    else
        return([fib.function_values[n_real_solutions] for fib in F])
    end
    
end

const FIBRE_DATA = EnumerativeProperty{Vector{FibreDatum}}("fibre_data")


ALGORITHM_DATA[explore_reality] = AlgorithmDatum(
    name = "explore_reality",
    description = "Explores the reality of solutions in an EnumerativeProblem by sampling parameters and solving for real solutions. It attempts to certify a witness for each unique count of real solutions found.",
    input_properties = [SYSTEM, BASE_FIBRE],
    output_property = FIBRE_DATA,
    reliability = :certified,
    automated = false
)



function find_all_real_counts(EP::EnumerativeProblem)
    E = explore(EP, n_real_solutions; sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples = 100, real_parameters = true, as_fibres = true)
    real_counts = sort(unique([f.function_values[n_real_solutions] for f in E]))
    @vprintln("Found the following real solution counts: $real_counts")
    possible_real_counts = degree(EP) % 2:2:degree(EP)
    @vprintln("Possible real solution counts: $possible_real_counts")
    missing_counts = setdiff(possible_real_counts, real_counts)
    found_fibres = []
    if length(missing_counts) > 0
        @vprintln("Missing real solution counts: $missing_counts")
        for k in missing_counts
            @vprint("Trying to find a witness for real solution count $k")
            O = find_k_real_solutions(EP, k; n_samples = 20, max_steps=10)
            if n_real_solutions(record_fibre(O)[1]) == k
                @vprintln("Success: Found a real solution count $k")
                push!(found_fibres, record_fibre(O))
                push!(E,FibreDatum(record_fibre(O)))
            else
                @vprintln("Failed to find a real solution count $k")
            end
        end
    end
    found_counts = [n_real_solutions(f) for f in found_fibres]
    still_missing_counts = setdiff(missing_counts,found_counts)
    all_found_counts = unique(vcat(found_counts,real_counts))
    certified_fibres = []
    for c in all_found_counts
        cfib = E[findfirst(x->n_real_solutions(x)==c, E)]
        certify!(cfib,EP)
        push!(certified_fibres,cfib)
        cfib.function_values[n_real_solutions] = n_real_solutions(cfib)
    end
    sort!(certified_fibres, by = x -> n_real_solutions(x))
    return(certified_fibres)
end