export
    explore


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
    explore(EP::EnumerativeProblem, f::Function; n_samples = 1000)
"""
function explore(EP::EnumerativeProblem, f::Function; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), kwargs...)
    # Sample points in the parameter space
    P = sampler(n_samples)
    # Solve EP for each sampled point
    sols = [EP(p) for p in P]
    # Apply the function to each solution
    results = [f(s) for s in sols]
    tally(results)
end

function explore(EP::EnumerativeProblem, F::Vector{Function}; n_samples = 1000, sampler = UniformSampler{Float64}(n_parameters(EP)), kwargs...)
    # Sample points in the parameter space
    P = sampler(n_samples)
    # Solve EP for each sampled point
    sols = [EP(p) for p in P]
    # Apply each function to each solution
    results = [tuple(f(s) for f in F) for s in sols]
    return([tally([f(s) for s in sols]) for f in F])
end