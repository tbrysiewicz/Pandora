
export 
    is_lacunary

function is_lacunary(EP::EnumerativeProblem)
    """
    Checks if the natural sparse embedding of the enumerative problem is into a lacunary sparse system. 
    """
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

