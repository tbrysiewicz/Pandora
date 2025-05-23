###The mutable struct OptimizerData 
@kwdef mutable struct OptimizerData
    step :: Int = 0
    parameters_solved :: Int = 0
    steps_no_progress :: Int = 0
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end

### The mutable struct ScoringScheme
mutable struct ScoringScheme
    objective                       #function we really want to optimize
    barrier                         #barrier penalty
    barrier_weight :: Float64       #weight of barrier penalty
    taboo
    error_checker
    goal
    name :: String
end



### The mutable struct Optimizer:
@kwdef mutable struct Optimizer
    EP::EnumerativeProblem

    solver_fibre :: Fibre

    record_fibre :: Fibre
    record_objective :: Any
    path :: Vector{Vector{ComplexF64}}

    sampler :: Sampler                       #used to sample parameter space (frequently updated)
    scoring_scheme :: ScoringScheme
    optimizer_data :: OptimizerData


    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, SS :: ScoringScheme)
        optimizer = new()

        optimizer.scoring_scheme = SS

        optimizer.EP = EP
        optimizer.solver_fibre = base_fibre(EP)
        optimizer.sampler=sampler
        p = base_parameters(EP)
        
        if is_real(sampler)
            p = Vector{ComplexF64}(real(p))
        end


        optimizer.record_fibre = (EP(p),p)
        optimizer.record_objective = SS.objective(optimizer.record_fibre)
        optimizer.path=[p]

        optimizer.optimizer_data = OptimizerData()

        return(optimizer)
    end
    
    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, objective)
        SS = ScoringScheme(objective)
        Optimizer(EP,sampler,SS)
    end

end



function ScoringScheme(objective; 
    barrier = zero_function, barrier_weight = 0.0, 
    taboo = zero_function, error_checker = false_function,
    goal = nothing, name = "")
    function more_than_one_hundred_steps(O::Optimizer)
        step(optimizer_data(O)) >100
    end
    if goal == nothing
        ScoringScheme(objective,barrier,barrier_weight,taboo,error_checker, more_than_one_hundred_steps,name)
    else
        ScoringScheme(objective,barrier,barrier_weight,taboo,error_checker, goal,name)
    end
end














