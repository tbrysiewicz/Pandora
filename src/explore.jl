export
    explore,
    reality_exploration


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
                    print("    check")
                else

                    print("    failed")
                end
            end
        else
            counts[v] = 1
            if certify_function !== nothing 
                # If a certify function is provided, use it to certify the solution
                @vprint("Certifying that function value $v is found in fibre number $(i)")
                if certify_function(Fibres[i],v)
                    witnesses[v] = Fibres[i][2]  # Store the first occurrence as a witness
                    print("    check")
                else

                    print("    failed")
                end
            else
                witnesses[v] = Fibres[i][2]  # Store the first occurrence as a witness
            end
            print("\n")
        end
    end
    return counts, witnesses
end

"""
    explore(EP::EnumerativeProblem, f::Function; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), kwargs...)

Explores the parameter space of an EnumerativeProblem `EP` by sampling parameters, solving `EP` over each parameter, and applying `f`. 
It returns a count of the unique values returned by `f` and a witness for each unique value.
The function accepts the following keyword arguments:
- `n_samples`: The number of samples to draw from the parameter space (default is 1000).
- `sampler`: A sampler function to generate samples from the parameter space (default is `UniformSampler{Float64}(n_parameters(EP))`).  

Note: The function filters out fibres that are not generic according to the `valid_fibre` function. 
"""
function explore(EP::EnumerativeProblem, f::Function; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), kwargs...)
    # Sample points in the parameter space
    P = sampler(n_samples)
    # Solve EP for each sampled point
    fibres = map(Fibre,collect(zip([EP(p) for p in P],P)))
    filter!(fibre -> valid_fibre(EP, fibre), fibres)
    return tally_and_witness_fibres(fibres, f)
end

function explore(EP::EnumerativeProblem, F::Vector{Function}; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), kwargs...)
    # Sample points in the parameter space
    P = sampler(n_samples)
    # Solve EP for each sampled point
    fibres = map(Fibre,collect(zip([EP(p) for p in P],P)))
    filter!(fibre -> valid_fibre(EP, fibre), fibres)
    return [tally_and_witness_fibres(fibres, f) for f in F]
end


function reality_exploration(EP::EnumerativeProblem; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), certify = true, kwargs...)
    certification_function = (fibre,n) -> true  # Default certification function that always returns true
    if certify
        certification_function = (fibre,n) -> certify_n_real(EP, fibre) == n
    end
    # Sample points in the parameter space
    P = sampler(n_samples)
    # Solve EP for each sampled point
    fibres = map(Fibre,collect(zip([EP(p) for p in P],P)))
    filter!(fibre -> valid_fibre(EP, fibre), fibres)
    (T,W) = tally_and_witness_fibres(fibres, n_real_solutions; certify_function = certification_function)
end

const FIBRE_DATUM = EnumerativeProperty{Dict{Symbol, Any}}("fibre_datum")

function reality_exploration!(EP::EnumerativeProblem; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), certify = true, kwargs...)
    (T,W) = reality_exploration(EP; n_samples = n_samples, sampler = sampler, certify = certify, kwargs...)

    ParameterKnowledge = Vector{KnowledgeNode{Dict{Symbol, Any}}}([])
    for (v,P) in W

        KN = KnowledgeNode{Dict{Symbol, Any}}(
            FIBRE_DATUM,
            Dict{Symbol, Any}(
                :parameter => P,
                :n_real_solutions => v
            ),
            Vector{KnowledgeNode}(),
            Dict{Symbol, Any}(),
            reality_exploration
        )
        push!(EP.knowledge, KN)
    end
end

