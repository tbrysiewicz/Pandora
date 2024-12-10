

function initialize_real_optimizer(EP::EnumerativeProblem)
    S = initialize_real_sampler(n_parameters(EP),2*n_parameters(EP))
    O = Optimizer(EP,S,dietmaier_pair)
    O.taboo = x->-n_real_solutions(x)
    O.error_checker = x->!(valid_fibre(EP,x)&&valid_real_fibre(EP,x))
    return(O)
end