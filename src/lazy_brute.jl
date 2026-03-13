export lazy_brute

#Lazy_brute takes an enumerative problem and tries to find as many real solutions as possible to it.
# it does this by choosing random parameters P=[p1,...,pN] and instead of solving for the entire
# fibre, it only tracks a random start solutin to a target solution using EP(pi,rand(1:degree(EP))
# If the solution is real, then that parameter is saved in a list of promising parameters
# Then, for the list of promising parameters, another random solution is tracked. If that is real, it 
# is added to a list of even MORE promising parameters (perhaps it makes sense to have a dictionary
#  like D[k] = list of parameters which have found k real solutions?
# The function should work in parallel. 
function lazy_brute(EP::EnumerativeProblem, N::Int; n_solves=100000, break_at = 5, record_check = 2, local_search = false, boolean_function = Pandora.is_real)
    function count_function(S)
        count = 0
        for s in S
            if boolean_function(s)
                count += 1
            end
        end
        return count
    end
    record_reals = 0
    record_param = Vector{Float64}(randn(Float64, n_parameters(EP)))
    d = Pandora.degree(EP)
    by_k = Dict{Int, Vector{Vector{Float64}}}(0 => [local_search ? randn(Float64, n_parameters(EP)) + record_param : randn(Float64, n_parameters(EP)) for _ in 1:N])
    for i in 1:break_at
        by_k[i] = Vector{Vector{Float64}}()
    end

    active_keys = filter(k -> !isempty(get(by_k, k, [])), keys(by_k))
    while n_solves > 0 && !isempty(active_keys) && maximum(active_keys) < break_at


        for k in 0:break_at
            println("k = $k: $(length(get(by_k, k, []))) parameters")
        end

        # Get the top key (highest k value) for which there are parameters to test - i.e. nonempty values
        active_keys = filter(k -> !isempty(get(by_k, k, [])), keys(by_k))
        if isempty(active_keys)
            break
        end
        top_k = maximum(active_keys)


        kparameters = by_k[top_k]
        #Reset this bucket; parameters will be re-added as needed
        by_k[top_k] = Vector{Vector{Float64}}()
        
        if top_k >= record_check
            all_solutions = solutions.(first.(EP(kparameters)))
            m = max([count_function(s) for s in all_solutions]...)
            first_best = findfirst(s -> count_function(s) == m, all_solutions)
            first_best_param = kparameters[first_best]
            if m > record_reals
                record_reals = m
                record_param = first_best_param
                println("New record of $record_reals real solutions found for parameters: $record_param")
            else
                println("Checked parameters with $top_k real solutions, but no new record found. Current record is $record_reals real solutions for parameters: $record_param")
            end
        end

        println("Solving for $(length(kparameters)) parameters with $top_k real solutions...")
        S = EP(kparameters,rand(1:d))
        for (i,s) in enumerate(S)
            if boolean_function(s)
                push!(by_k[top_k + 1], kparameters[i])
            end
        end
        
        n_solves -= length(S)
        # Replenish if dictionary is under N parameters
        total_params = sum(length(by_k[k]) for k in keys(by_k))
        if total_params < N
            by_k[0] = vcat(get(by_k, 0, []),  [local_search ? 0.1*randn(Float64, n_parameters(EP)) + record_param : randn(Float64, n_parameters(EP)) for _ in 1:N])
        end
        active_keys = filter(k -> !isempty(get(by_k, k, [])), keys(by_k))
        println("$n_solves solves remaining. Total parameters in by_k: $total_params")
        println("Record: $record_reals")
    end
    #Either we have exhausted our n_solves or we have found 5 real solutions for some parameters.
    #If the latter, we solve for all parameters in that list to find the real solutions.
    active_keys = filter(k -> !isempty(get(by_k, k, [])), keys(by_k))
    if !isempty(active_keys) && maximum(active_keys) == break_at
        promising_params = by_k[break_at]
        println("Found $(length(promising_params)) parameters with $break_at real solutions. Solving for all of them...")
        all_solutions = solutions.(first.(EP(promising_params)))
        #find the one with the most real solutions and return the Fibre 
        m = max([count_function(s) for s in all_solutions]...)
        first_best = findfirst(s -> count_function(s) == m, all_solutions)
        first_best_param = promising_params[first_best]
        first_best_solutions = all_solutions[first_best]
        if m == d 
            println("Found parameters with all solutions real! Returning one of those...")
        else
            println("Most real solutions we found was $m out of $d. Returning one of those...")
        end

        return(Fibre(first_best_solutions, first_best_param))
    else
        println("Exhausted n_solves without finding parameters with $break_at real solutions. Returning the best we found...")
        #still take the bucket of highest k and return the one with the most like in earlier lines
        if isempty(active_keys)
            promising_params = get(by_k, 0, Vector{Vector{Float64}}())
            println("No parameters with real solutions found. Returning a random parameter...")
        else
            promising_params = by_k[maximum(active_keys)]
        end
        if isempty(promising_params)
            error("No promising parameters remain in by_k.")
        end
    end
    return(record_param)
end