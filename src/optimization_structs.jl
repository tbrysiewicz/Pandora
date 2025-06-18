
#The struct OptimizerData : keeps track of the Data involved while using the optimizer.

@kwdef mutable struct OptimizerData
    step :: Int = 0                             #Step count 
    parameters_solved :: Int = 0                
    steps_no_progress :: Int = 0                
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end

##Getters for the OptimizerData

step(OD::OptimizerData) = OD.step
parameters_solved(OD::OptimizerData) = OD.parameters_solved
steps_no_progress(OD::OptimizerData) = OD.steps_no_progress
steps_no_major_progress(OD::OptimizerData) = OD.steps_no_major_progress
error_proportion(OD::OptimizerData) = OD.error_proportion
taboo_proportion(OD::OptimizerData) = OD.taboo_proportion
improvement_proportion(OD::OptimizerData) = OD.improvement_proportion

##Setters for OptimizerData

function set_step!(OD::OptimizerData, new_step:: Int)
    OD.step = new_step
    return nothing
end

function set_parameters_solved!(OD::OptimizerData, new_params_solved:: Int)
    OD.parameters_solved = new_params_solved
    return nothing
end

function set_steps_no_progress!(OD::OptimizerData, new_steps_no_progress :: Int)  
    OD.steps_no_progress = new_steps_no_progress
    return nothing
end              
    
function set_steps_no_major_progress!(OD::OptimizerData, new_steps_no_major_progress :: Int)
    OD.steps_no_major_progress = new_steps_no_major_progress
    return nothing
end
   
function set_error_proportion!(OD::OptimizerData, new_error_proportion:: Float64)
    OD.error_proportion = new_error_proportion
    return nothing
end

function set_taboo_proportion!(OD::OptimizerData, new_taboo_proportion::Float64)
    OD.taboo_proportion = new_taboo_proportion
    return nothing
end

function set_improvement_proportion!(OD::OptimizerData, new_improvement_proportion :: Float64)
    OD.improvement_proportion = new_improvement_proportion
    return nothing
end



#The mutable Struct ScoringScheme

mutable struct ScoringScheme
    objective                               # The function that we really want to be optimized.
    barrier
    barrier_weight :: Float64               # weight of barrier penalty
    taboo
    error_checker                  
    goal                                    # goal (:: Union{Bool, Nothing}) that stops 
                                            # the optimizer, when it achieves the value "true".
    name :: String                          # name of the ScoringScheme.
end 


# Base.show for ScoringScheme

function Base.show(io::IO, SS::ScoringScheme)
    if SS.name==""
        print(io,"A scoring scheme with barrier weight ",SS.barrier_weight,".")
    else
        print(io,"The ",SS.name," scoring scheme with barrier weight ",SS.barrier_weight, ".")
    end
end


# Getters for ScoringScheme

objective(SS:: ScoringScheme) = SS.objective
barrier(SS::ScoringScheme) = SS.barrier
barrier_weight(SS::ScoringScheme) = SS.barrier_weight
taboo(SS::ScoringScheme) = SS.taboo
error_checker(SS::ScoringScheme) = SS.error_checker
goal(SS::ScoringScheme) = SS.goal
name(SS::ScoringScheme) = SS.name

### Function weighted objective:

function weighted_objective(SS::ScoringScheme)
    SS_objective = objective(SS)
    SS_barrier_weight = barrier_weight(SS)
    SS_barrier = barrier(SS)
    x->(1-SS_barrier_weight).*SS_objective(x).+(SS_barrier_weight).*SS_barrier(x)
end

## Setters for ScoringScheme to be introduced.






# The mutable struct Optimizer:

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
        optimizer.record_objective = SS.objective(optimizer.record_fibre[1])
        optimizer.path=[p]

        optimizer.optimizer_data = OptimizerData()

        return(optimizer)
    end
    
    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, objective; goal = nothing)
        SS = ScoringScheme(objective, goal = goal)
        Optimizer(EP,sampler,SS)
    end

end

# Base.show for Optimizer

function Base.show(io::IO, optimizer::Optimizer)
    tenspaces="          "
    println(io,"Optimizer setup for an enumerative problem of degree ",degree(optimizer.EP),".")
    print(io,"\n\n")
    println(io,tenspaces,tenspaces,"Progress Data",tenspaces)
    println("---------------------------------------------------------")
    println(io,"  Current score: ",optimizer.record_objective)
    od = optimizer.optimizer_data
    println("                         Step #:  ",od.step)
    println("        Total parameters solved:  ",od.parameters_solved)
    println("Steps since last minor progress:  ",od.steps_no_progress)
    println("      since last major progress:  ",od.steps_no_major_progress)
    println("          Last error proportion:  ",od.error_proportion)
    println("          Last taboo proportion:  ",od.taboo_proportion)
    println("    Last improvement proportion:  ",od.improvement_proportion)
    println("                         Radius:  ", maximum(abs, vec(optimizer.sampler.transform_matrix)))
    #=println("The parameters of the record fibre is:")
    display(optimizer.record_fibre[2])
    =#
    #info regarding radius or last improvement size 
    #info regarding conditioning of the ellipse sampler
    #interpretive information 
    #soft-goal info
end


##The getters for the struct Optimizer.

enumerative_problem(optimizer::Optimizer) = optimizer.EP
solver_fibre(optimizer::Optimizer) = optimizer.solver_fibre
record_fibre(optimizer::Optimizer) = optimizer.record_fibre
record_objective(optimizer::Optimizer) = optimizer.record_objective
path(optimizer::Optimizer) = optimizer.path
sampler(optimizer::Optimizer) = optimizer.sampler
scoring_scheme(optimizer::Optimizer) = optimizer.scoring_scheme
optimizer_data(optimizer::Optimizer) = optimizer.optimizer_data

### For sampler:
sample(optimizer::Optimizer) = sample(sampler(optimizer))

### For scoring_scheme:
objective_function(optimizer::Optimizer) = objective(scoring_scheme(optimizer))
barrier(optimizer::Optimizer) = barrier(scoring_scheme(optimizer))
barrier_weight(optimizer::Optimizer) = barrier_weight(scoring_scheme(optimizer))
taboo(optimizer::Optimizer) = taboo(scoring_scheme(optimizer))
error_checker(optimizer::Optimizer) = error_checker(scoring_scheme(optimizer))
goal(optimizer::Optimizer) = goal(scoring_scheme(optimizer))
name(optimizer::Optimizer) = name(scoring_scheme(optimizer))


weighted_objective(optimizer::Optimizer) = weighted_objective(scoring_scheme(optimizer))


## Constructor for ScoringScheme
function ScoringScheme(objective; 
    barrier = zero_function,
    barrier_weight = 0.0, 
    taboo = zero_function, 
    error_checker = false_function,
    goal = nothing, 
    name = "")

    function more_than_one_hundred_steps(O::Optimizer)
        step(optimizer_data(O)) > 100
    end

    if goal == nothing
        ScoringScheme(objective,
                      barrier,
                      barrier_weight,
                      taboo,
                      error_checker, 
                      more_than_one_hundred_steps,
                      name)
    else
        ScoringScheme(objective,
                      barrier,
                      barrier_weight,
                      taboo,
                      error_checker, 
                      goal,
                      name)
    end
end

