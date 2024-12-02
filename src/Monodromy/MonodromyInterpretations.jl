function is_decomposable(EP::EnumerativeProblem)
    !is_primitive(galois_group(EP))
end

function is_irreducible(EP::EnumerativeProblem)
    is_transitive(galois_group(EP))
end