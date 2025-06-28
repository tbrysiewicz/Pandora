export 
    ScoringScheme




"""
    The struct ScoringScheme collects the information required to run optimization code in Pandora.
    Its fields are:
    - `objective`: A function from `Fibre` to some type `T` where elements of `T` can be compared. This function represents what the user truly wants to optimize locally.
    - `barrier`: A function from `Fibre` to `T`, which is used to penalize the objective function.
    - `barrier_weight`: A `Float64` between 0 and 1 which determines the weight of the barrier penalty
    - `taboo`: A function from `Fibre` to `S`, where elements of `S` can be compared.
    - `goal`: A function from `Optimizer` to `Bool` 
    - `name`: A `String` that names the scoring scheme.
    The function `improve(EP::EnumerativeProblem, F::Fibre, SS:ScoringScheme; Sampler, n_samples)` 
        will sample `n_samples` new parameters, using `Sampler`, solve for the fibres `FF` of the EnumerativeProblem `EP` for these parameters,
        and return the fibre of `FF` which is "best" according to the `ScoringScheme` `SS`.
    A fibre `F'` is better than a fibre `F` according to the `ScoringScheme` `SS` if:
    - `error_checker(SS)(F')` is `false` and `error_checker(SS)(F)` is `true`, or if both are `false`, then 
    - `taboo(SS)(F')` is less than `taboo(SS)(F)`, or if both are equal, then
    - `weighted_objective(SS)(F')`=`(1-barrier_weight)*objective(SS)(F')+barrier_weight*barrier(SS)(F')` is greater than `weighted_objective` of `F`.
"""
mutable struct ScoringScheme
    objective                               # The function that we really want to be locally optimized.
    barrier
    barrier_weight :: Float64               # weight of barrier penalty
    taboo            
    name :: String                          # name of the ScoringScheme.
end 


# Base.show for ScoringScheme
function Base.show(io::IO, SS::ScoringScheme)
    print(io,"A scoring scheme.")
end


#=

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

## Setters for ScoringScheme.
set_objective!(SS::ScoringScheme, new_objective) = (SS.objective = new_objective; nothing)
set_barrier!(SS::ScoringScheme, new_barrier) = (SS.barrier = new_barrier; nothing)
set_barrier_weight!(SS::ScoringScheme, new_barrier_weight::Float64) = (SS.barrier_weight = new_barrier_weight; nothing)
set_taboo!(SS::ScoringScheme, new_taboo) = (SS.taboo = new_taboo; nothing)
set_error_checker!(SS::ScoringScheme, new_error_checker) = (SS.error_checker = new_error_checker; nothing)
set_goal!(SS::ScoringScheme, new_goal::Union{Bool, Nothing}) = (SS.goal = new_goal; nothing)
set_name!(SS::ScoringScheme, new_name::String) = (SS.name = new_name; nothing)




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
set_step!(OD::OptimizerData, new_step::Int) = (OD.step = new_step; nothing)
set_parameters_solved!(OD::OptimizerData, new_params_solved::Int) = (OD.parameters_solved = new_params_solved; nothing)
set_steps_no_progress!(OD::OptimizerData, new_steps_no_progress::Int) = (OD.steps_no_progress = new_steps_no_progress; nothing)
set_steps_no_major_progress!(OD::OptimizerData, new_steps_no_major_progress::Int) = (OD.steps_no_major_progress = new_steps_no_major_progress; nothing)
set_error_proportion!(OD::OptimizerData, new_error_proportion::Float64) = (OD.error_proportion = new_error_proportion; nothing)
set_taboo_proportion!(OD::OptimizerData, new_taboo_proportion::Float64) = (OD.taboo_proportion = new_taboo_proportion; nothing)
set_improvement_proportion!(OD::OptimizerData, new_improvement_proportion::Float64) = (OD.improvement_proportion = new_improvement_proportion; nothing)


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
        optimizer.record_objective = objective(SS)(EP(p))
        optimizer.path=[p]

        optimizer.optimizer_data = OptimizerData()

        return(optimizer)
    end
    
    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, objective; goal = O-> step(optimizer_data(O)) > 100)
        SS = ScoringScheme(objective, goal = goal)
        Optimizer(EP,sampler,SS)
    end

end


function more_than_one_hundred_steps(O::Optimizer)
    step(optimizer_data(O)) > 100
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

#Setters
set_enumerative_problem!(optimizer::Optimizer, new_enumerative_problem::EnumerativeProblem) = (optimizer.EP = new_enumerative_problem; nothing)
set_solver_fibre!(optimizer::Optimizer, new_solver_fibre::Fibre) = (optimizer.solver_fibre = new_solver_fibre; nothing)
set_record_fibre!(optimizer::Optimizer, new_record_fibre::Fibre) = (optimizer.record_fibre = new_record_fibre; nothing)
set_record_objective!(optimizer::Optimizer, new_record_objective) = (optimizer.record_objective = new_record_objective; nothing)
set_path!(optimizer::Optimizer, new_path::Vector{Vector{ComplexF64}}) = (optimizer.path = new_path; nothing)
set_sampler!(optimizer::Optimizer, new_sampler::Sampler) = (optimizer.sampler = new_sampler; nothing)
set_scoring_scheme!(optimizer::Optimizer, new_scoring_scheme::ScoringScheme) = (optimizer.scoring_scheme = new_scoring_scheme; nothing)
set_optimizer_data!(optimizer::Optimizer, new_optimizer_data::OptimizerData) = (optimizer.optimizer_data = new_optimizer_data; nothing)


##Sampler:


## Constructor for ScoringScheme
function ScoringScheme(objective; 
    barrier = zero_function,
    barrier_weight = 0.0, 
    taboo = zero_function, 
    error_checker = false_function,
    goal = nothing, 
    name = "")


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

export dietmaier_scheme, initialize_real_optimizer

function dietmaier_scheme(EP::EnumerativeProblem) :: ScoringScheme
    objective = dietmaier_pair
    taboo = x -> -n_real_solutions(x)
    error_checker = x->!(valid_fibre(EP,x)&&valid_real_fibre(EP,x))

    d = degree(EP)

    function totally_real(optimizer::Optimizer)
        record_objective(optimizer)[1]==d
    end 

    
    barrier = x->(n_real_solutions(x),0.0) ##replace with 1/(min dist between real solutions)
    barrier_weight = 0.001

    SS = ScoringScheme(objective; barrier=barrier, barrier_weight=barrier_weight, 
                            taboo=taboo, error_checker=error_checker, 
                            goal = totally_real, name = "Dietmaier")
    return(SS)
end

# The function find_error_fibres uses the error_checker function from the optimizer to check for 
# errors and return the fibres that has an error.

function find_error_fibres(optimizer::Optimizer, 
                           new_fibres, 
                           N)

    error_fibres = findall(i->error_checker(optimizer)(new_fibres[i]), 1:N)
    return (error_fibres)
end

#The function find_taboo_fibres finds the fibres that are taboo, from the set of all the new_fibres.

function find_taboo_fibres(optimizer::Optimizer, 
                           new_fibres::Vector{Fibre},  
                           error_fibres, 
                           N)

     #Gets the taboo function.
     taboo_function = taboo(optimizer)

     #Gets the current taboo_score.
     current_taboo_score = taboo_function(record_fibre(optimizer)) 

     #We evaluate the taboo scores for all the non error_fibres, and store it as a list.
     taboo_score_list = [in(i,error_fibres) ? nothing : taboo_function(new_fibres[i]) for i in 1:N]
     #Finds the indices for which the taboo_function is more than the allowed current_taboo_score, and store it as a list.
     taboo_fibres = findall(i->(in(i,error_fibres)==false && taboo_score_list[i]>current_taboo_score), 1:N)
     
     return (current_taboo_score, 
             taboo_score_list, 
             taboo_fibres)
end

# The function find_weighted_score_list finds the 

function find_weighted_score_list(optimizer::Optimizer,
                                  new_fibres, 
                                  error_fibres, 
                                  taboo_fibres, 
                                  N)

    #Gets the weighted objective.
    w_obj = weighted_objective(optimizer)

    #Evaluates the weighted objective of the record_fibre.
    wo_record = w_obj(record_fibre(optimizer)[1])

    #Evaluates the weighted_objective of the new fibres, as long as a fibre is not in taboo_fibres or error_fibres. 
    wo_score_list = [in(i,union(taboo_fibres,error_fibres)) ? nothing : w_obj(new_fibres[i][1]) for i in 1:N]
    
    return (wo_record, wo_score_list)
end

#The function find_improvement_fibres finds the fibres that have "improved", from the set of all the new_fibres.

function find_improvement_fibres(error_fibres, 
                                 current_taboo_score,
                                 taboo_score_list,  
                                 taboo_fibres, 
                                 wo_score_list, 
                                 wo_record, 
                                 N)

    #Finds the non-error, non-taboo fibres for which the w_obj have improved.
    improvement_fibres = findall(i->in(i,union(error_fibres, taboo_fibres))==false && wo_score_list[i]>wo_record,1:N)
   
    #Among the improvement_fibres we seperate the ones whose taboo scores have decreased and stayed the same into
    # major_improvement_fibres and minor_improvement_fibres respectively.
    major_improvement_fibres = findall(i->in(i,improvement_fibres) && taboo_score_list[i]<current_taboo_score,1:N)
    minor_improvement_fibres = findall(i->in(i,improvement_fibres) && taboo_score_list[i]==current_taboo_score, 1:N)
    
    #And the non-error non-taboo fibres that showed no improvement is collected in no_improvement_fibres.
    no_improvement_fibres = findall(i->in(i,union(taboo_fibres,error_fibres,improvement_fibres))==false, 1:N)
    
    return (improvement_fibres, 
            major_improvement_fibres, 
            minor_improvement_fibres, 
            no_improvement_fibres)
end

function find_best_score_index(wo_score_list,  
                               improvement_fibres)
    best_score_index = []

    #Making sure that improvement_fibres are not empty
    if length(improvement_fibres)>0

        #filtering out the scores of the non-error, non-taboo fibres and their indices.
        true_scores = filter(x->x!=nothing, wo_score_list)
        true_score_indices = findall(x->x!=nothing, wo_score_list)

        #@assert throws an error if it is empty
        @assert length(true_score_indices)>0
        #Find the index (a.k.a arg) of the maximum in true_scores.
        b = argmax(true_scores)                     # Maybe need to change argmax to something else to include the ToSet functionality. Check this.
        
        #using that to get the index of the maximum in wo_score_list.
        best_score_index = [true_score_indices[b]]
    end
    return (best_score_index)
end

# This is the function which increments all relevant step datum in an optimizer data.
# It is used in update_optimizer! which is used in improve!.
function increment!(OD::OptimizerData, N::Int) 
    set_step!(OD, step(OD)+1)                                            #Setters to be added.
    set_steps_no_major_progress!(OD, steps_no_major_progress(OD)+1)
    set_steps_no_progress!(OD, steps_no_progress(OD)+1)
    set_parameters_solved!(OD, parameters_solved(OD) + N)
end

# The function partition_samples partitions the collection of samples into collections of 
# error_fibres (which satisfies the error condition and thus causes errors),
# taboo_fibres (which are non-error fibres that has taboo beyond the accepted limit),
# improvement_fibres (which are non-error, non-taboo fibres that does give us an improvement) and
# no_improvement_fibres (which are non-error, non-taboo fibres that does NOT give us an improvement).

function partition_samples(optimizer:: Optimizer, new_fibres :: Vector{Fibre}) 

    N = length(new_fibres)

    # Finds the error_fibres.
    (error_fibres) = find_error_fibres(optimizer,
                                     new_fibres,
                                     N)
    
    # Finds the taboo_fibres.
    (current_taboo_score,
     taboo_score_list,
     taboo_fibres) = find_taboo_fibres(optimizer, 
                                       new_fibres, 
                                       error_fibres, 
                                       N)

    # Finds the weighted objective scores of the fibres as wo_score_list.
    (wo_record,
     wo_score_list) = find_weighted_score_list(optimizer, 
                                               new_fibres, 
                                               error_fibres, 
                                               taboo_fibres, 
                                               N)
    
     # Finds the improvement fibres; and also its partition into 
     # major_improvement_fibres and
     # minor_improvement_fibres.
    (improvement_fibres,
     major_improvement_fibres,
     minor_improvement_fibres,
     no_improvement_fibres) = find_improvement_fibres(error_fibres, 
                                                      current_taboo_score,
                                                      taboo_score_list, 
                                                      taboo_fibres, 
                                                      wo_score_list, 
                                                      wo_record, 
                                                      N)
    

    (best_score_index) = find_best_score_index(wo_score_list, improvement_fibres)

    return (error_fibres,
            taboo_fibres,
            minor_improvement_fibres,
            major_improvement_fibres,
            no_improvement_fibres,
            best_score_index)
end



import LinearAlgebra: diagm, svd

# function update_optimizer! updates the fields of the optimizer::Optimizer according
# to the new fibres; especially the fields of the optimizer_data::OptimizerData.

function update_optimizer!(optimizer :: Optimizer, new_fibres :: Vector{Fibre})
    OD = optimizer_data(optimizer)
    N = length(new_fibres)

    #Make a default increment update to optimizer data (steps += 1, params solved +=N)
    increment!(OD,N)

    #Get the values of the folllowing using the function partition_samples
    (error_fibs,
     taboo_fibs,                          #Check where this is used again.
     min_improv, 
     maj_improv,
     non_improv, 
     best_score_index) = partition_samples(optimizer,new_fibres)
    
    #Updating OD, to include the proportion of error fibres, taboo fibres, and
    #improved fibres to the total number of fibres N.
    set_error_proportion!(OD, length(error_fibs)/N)                           #Setters to be added.
    set_taboo_proportion!(OD, length(taboo_fibs)/N)
    set_improvement_proportion!(OD, (length(min_improv)+length(maj_improv))/N)

    #Update record fibre/solver fibre
    if improvement_proportion(OD)>0.0
        set_steps_no_progress!(OD,0)                                           #Setters to be added.
        println("Minor progress")
        if length(maj_improv)>0
            set_steps_no_major_progress!(OD, 0)                               #Setters to be added.
            println("Major progress")
        end
        optimizer.sampler.translation = 2*new_fibres[best_score_index[1]][2] - record_fibre(optimizer)[2]
        optimizer.record_fibre = new_fibres[best_score_index[1]]         #Setters to be added.
        optimizer.record_objective = weighted_objective(optimizer)(record_fibre(optimizer)[1])
    end

    #shaping the samples according to non-error non-taboo samples using SVD.
    #non_error_taboo_fibre_values = new_fibres[union(min_improv,
    #                                                maj_improv, 
    #                                                non_improv)]
    #A = reduce(hcat,([non_error_taboo_fibre_values[i][2] for i in 1:length(non_error_taboo_fibre_values)]))
    #u,s,v = svd(A)
    #optimizer.sampler.transform_matrix = u*diagm(s)
    
    
    
    #Updating the path
    optimizer.path = push!(path(optimizer), record_fibre(optimizer)[2])
    
    return optimizer
end



# update_sampler! - the function that updates the sampler. Specifically, this
# functions deals with updating the transform_matrix field of the sampler,
# 1) to make sure the new samples are biased towards the direction of the record,
# 2) to make sure that the proportions (such as taboo and errors) stay within the bounds we define.


function update_sampler!(optimizer,new_fibres)
    #Gets the optimizer_data and the sampler.
    OD = optimizer_data(optimizer)
    S = sampler(optimizer)

    #Proportion updates
    if error_proportion(OD)>0.9
        println("Lot's of errors - code adjustment of numerical tolerances")
    end
    if taboo_proportion(OD)>0.8
        println("Lot's of taboo fibres - reducing radius of sampler")
        set_transform_matrix!(S, transform_matrix(S)*0.9)
    end
    if taboo_proportion(OD)<0.2 && improvement_proportion(OD)>0.0
        println("Not enough taboo fibres - we can be bolder in our sampling -increase radius")
        set_transform_matrix!(S, transform_matrix(S)*1.1)
    end
    if improvement_proportion(OD)>0.2
        println("Lot's of improvement - increase radius")
        set_transform_matrix!(S, transform_matrix(S)*1.1)
    end
    if improvement_proportion(OD)<0.01 && steps_no_progress(OD)>5
        println("Almost no improvement - decrease radius - maybe we are at a local max")
        set_transform_matrix!(S, transform_matrix(S)*0.5)
    end
end


# The function improve! takes in optimizer::Optimizer, and updates it for 
# a single step of a new sample of parameters. The sampling is done by 
# getting the sampler from optimizer, and the OptimizerData and sampler
# values are respectively updated by the update_optimizer! and update_sampler!
# functions respectively.

function improve!(optimizer::Optimizer) 

    #Get relevant fields
    EP = enumerative_problem(optimizer)
    S = sampler(optimizer)

    #Sample new parameters
    new_samples = sample(S)
    N = n_samples(S)

    #Solve for the samples
    new_solutions = EP(new_samples)
    new_fibres = [(new_solutions[i],new_samples[i]) for i in 1:N]
    
    #Update optimizer data values and record fibre
    update_optimizer!(optimizer,new_fibres)
    #Update sampler and solver fibre
    update_sampler!(optimizer,new_fibres)

    return(optimizer)
end


# The function optimize! updates the input O::Optimizer by applying the function 
# improve! to it until the goal is satisified.

 function optimize!(O::Optimizer; visualize_optimizer = false)
    G = goal(O)
    if visualize_optimizer && (length(parameters(enumerative_problem(O)))==2)
        #   visualize(enumerative_problem(O))[2]
        while (G(O)==false) 
            improve!(O)
            println(O)
        end
        _path = path(O)
        scatter!([real(_path[i][1]) for i in 1:length(_path)],
                 [real(_path[i][2]) for i in 1:length(_path)], 
                 zcolor = [i for i in 1:length(_path)], 
                 colormap = :viridis)
    #return(O)
    elseif visualize_optimizer  && (length(parameters(enumerative_problem(O))) !=2)
        println("Parameter space is not of dimension 2. Cannot be visualised")
        while (G(O)==false) 
            improve!(O)
            println(O)
        end
    return(O)
    else
        while (G(O)==false) 
            improve!(O)
            println(O)
        end
    return(O)
    end
end

function optimize(EP::EnumerativeProblem, SS::ScoringScheme; sampler = nothing, n_samples = 100,visualize_optimizer = false)#add kwarg for is_less
    if isnothing(sampler)
        sampler = initialize_real_sampler(n_parameters(EP), n_samples)
    end
    
    O = Optimizer(EP, sampler, SS::ScoringScheme)
    optimize!(O,visualize_optimizer = visualize_optimizer)
 end 


function optimize(EP::EnumerativeProblem, objective::Function; sampler = nothing, n_samples =100,visualize_optimizer = false)#add kwarg for is_less
    if isnothing(sampler)
        sampler = initialize_real_sampler(n_parameters(EP), n_samples)
    end
    
    O = Optimizer(EP, sampler, objective::Function)
    optimize!(O,visualize_optimizer = visualize_optimizer)
end

function optimize_n_real_solutions(EP::EnumerativeProblem)
    D = dietmaier_scheme(EP)
    optimize(EP, D)
end


#=

function optimize!(O::Optimizer)
    G = goal(O)
    while (G(O)==false)
            improve!(O)
            println(O)
    end
    return(O)
 end

=#















=#