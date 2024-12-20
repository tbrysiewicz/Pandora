function optimizer_data(optimizer::Optimizer)
    optimizer.optimizer_data
end

function sampler(optimizer::Optimizer)
    optimizer.sampler
end

function sample(optimizer::Optimizer)
    sample(optimizer.sampler)
end

function enumerative_problem(optimizer::Optimizer)
    optimizer.EP
end

function objective_function(optimizer::Optimizer)
    optimizer.scoring_scheme.objective
end

function record_objective(optimizer::Optimizer)
    optimizer.record_objective
end

function error_checker(optimizer::Optimizer)
    optimizer.scoring_scheme.error_checker
end

function taboo(optimizer::Optimizer)
    optimizer.scoring_scheme.taboo
end

function record_fibre(optimizer::Optimizer)
    optimizer.record_fibre
end

function weighted_objective(optimizer::Optimizer)
    weighted_objective(optimizer.scoring_scheme)
end