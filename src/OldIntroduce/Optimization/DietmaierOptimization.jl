
function dietmaier_scheme(EP::EnumerativeProblem) :: ScoringScheme
    objective = dietmaier_pair
    taboo = x -> -n_real_solutions(x)
    error_checker = x->!(valid_fibre(EP,x)&&valid_real_fibre(EP,x))

    d = degree(EP)

    function totally_real(optimizer::Optimizer)
        optimizer.record_objective[1]==d
    end 

    
    barrier = x->(n_real_solutions(x),0.0) ##replace with 1/(min dist between real solutions)
    barrier_weight = 0.001

    SS = ScoringScheme(objective; barrier=barrier, barrier_weight=barrier_weight, 
                            taboo=taboo, error_checker=error_checker, 
                            goal = totally_real, name = "Dietmaier")
    return(SS)
end

function initialize_real_optimizer(EP::EnumerativeProblem)
    SS = dietmaier_scheme(EP)
    S = initialize_real_sampler(n_parameters(EP),2*n_parameters(EP))
    O = Optimizer(EP,S,SS)
    return(O)
end