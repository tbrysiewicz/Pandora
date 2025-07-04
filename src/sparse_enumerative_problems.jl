
export 
    is_lacunary,
    is_triangular

function is_lacunary(EP::EnumerativeProblem)
    poly_supports = support(EP)
    shifted_supports = [A .- A[:, 1] for A in poly_supports]
    #horizontally concatenate the shifted supports
    concatenated_supports = hcat(shifted_supports...)
    if prod(elementary_divisors(concatenated_supports)) == 1
        return false
    else
        return true
    end
end


function is_triangular(EP::EnumerativeProblem)
    poly_supports = newton_polytopes(EP)
    for c in combinations(poly_supports)
        if dim(sum(c)) <= length(c) && length(c) != n_variables(EP)
            return true
        end
    end
    return(false)
end