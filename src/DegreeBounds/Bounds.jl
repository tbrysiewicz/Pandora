
function compute_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=true)
end

function compute_affine_bkk(EP::EnumerativeProblem)
  paths_to_track(specialize(system(EP));only_torus=false)
end

function compute_bezout(EP::EnumerativeProblem)
    return(prod(degree_sequence(EP)))
end

function degree_sequence(EP::EnumerativeProblem)
    Ms = support_coefficients(system(EP))[1]
    deg_seq = [max(map(x->sum(x),eachcol(M))...) for M in Ms]
    return(deg_seq)
end