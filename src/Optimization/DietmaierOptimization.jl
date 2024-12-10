

function initialize_real_optimizer(EP::EnumerativeProblem)
    S = initialize_real_sampler(n_parameters(EP),2*n_parameters(EP))
    O = Optimizer(EP,S,dietmaier_pair)
end