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
    optimizer.objective
end