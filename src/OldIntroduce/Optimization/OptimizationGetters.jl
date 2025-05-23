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

### The mutable struct SCoringScheme

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

### Getters for OptimizerData
#=
@kwdef mutable struct OptimizerData
    step :: Int = 0
    parameters_solved :: Int = 0
    steps_no_progress :: Int = 0
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end
=#

function step(OD::OptimizerData)
    OD.step
end

function parameters_solved(OD::OptimizerData)
    OD.parameters_solved
end

function steps_no_progress(OD::OptimizerData)
    OD.steps_no_progress
end

function steps_no_major_progress(OD::OptimizerData)
    OD.steps_no_major_progress
end

function error_proportion(OD::OptimizerData)
    OD.error_proportion
end

function taboo_proportion(OD::OptimizerData)
    OD.taboo_proportion
end

function improvement_proportion(OD::OptimizerData)
    OD.improvement_proportion
end

### Getters for the struct Scoring scheme:
#=
mutable struct ScoringScheme
    objective                       #function we really want to optimize
    barrier                         #barrier penalty
    barrier_weight :: Float64       #weight of barrier penalty
    taboo
    error_checker
    goal
    name :: String
end
=#

function objective(SS:: ScoringScheme)
    SS.objective
end

function barrier(SS::ScoringScheme)
    SS.barrier
end

function barrier_weight(SS::ScoringScheme)
    SS.barrier_weight
end

function taboo(SS::ScoringScheme)
    SS.taboo
end

function error_checker(SS::ScoringScheme)
    SS.error_checker
end

function goal(SS::ScoringScheme)
    SS.goal
end

function name(SS::ScoringScheme)
    SS.name
end


### Function weighted objective:

function weighted_objective(SS::ScoringScheme)
    SS_objective = objective(SS)
    SS_barrier_weight = barrier_weight(SS)
    SS_barrier = barrier(SS)
    x->(1-SS_barrier_weight).*SS_objective(x).+(SS_barrier_weight).*SS_barrier(x)
end


### Getters for the struct Optimizer:

#=
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
=#

function enumerative_problem(optimizer::Optimizer)
    optimizer.EP
end

function solver_fibre(optimizer::Optimizer)
    optimizer.solver_fibre
end

function record_fibre(optimizer::Optimizer)
    optimizer.record_fibre
end

function record_objective(optimizer::Optimizer)
    optimizer.record_objective
end

function path(optimizer::Optimizer)
    optimizer.path
end

function sampler(optimizer::Optimizer)
    optimizer.sampler
end

## For sampler:
function sample(optimizer::Optimizer)
    sample(sampler(optimizer))
end
##


function scoring_scheme(optimizer::Optimizer)
    optimizer.scoring_scheme
end

##For scoring_scheme:
function objective_function(optimizer::Optimizer)
    objective(scoring_scheme(optimizer))
end

function barrier(optimizer::Optimizer)
    barrier(scoring_scheme(optimizer))
end

function barrier_weight(optimizer::Optimizer)
    barrier_weight(scoring_scheme(optimizer))
end

function taboo(optimizer::Optimizer)
    taboo(scoring_scheme(optimizer))
end

function error_checker(optimizer::Optimizer)
    error_checker(scoring_scheme(optimizer))
end

function goal(optimizer::Optimizer)
    goal(scoring_scheme(optimizer))
end

function name(optimizer::Optimizer)
    name(scoring_scheme(optimizer))
end

#
function weighted_objective(optimizer::Optimizer)
    weighted_objective(scoring_scheme(optimizer))
end
##


function optimizer_data(optimizer::Optimizer)
    optimizer.optimizer_data
end






### Base.show for ScoringScheme

function Base.show(io::IO, SS::ScoringScheme)
    if SS.name==""
        print(io,"A scoring scheme with barrier weight ",SS.barrier_weight,".")
    else
        print(io,"The ",SS.name," scoring scheme with barrier weight ",SS.barrier_weight, ".")
    end
end

### Base.show for Optimizer

function Base.show(io::IO, optimizer::Optimizer)
    tenspaces="          "
    print(io,"Optimizer setup for an enumerative problem of degree ",degree(optimizer.EP),".")
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
    #info regarding radius or last improvement size 
    #info regarding conditioning of the ellipse sampler
    #interpretive information 
    #soft-goal info
end
