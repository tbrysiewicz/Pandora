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


# The function updates_from_proportions! updates the sampler to help make sure 
# the (error,taboo) proportions lay within the defined limits. It is called 
# in the function update_sampler!.

function updates_from_proportions!(optimizer::Optimizer)
    #Gets the optimizer_data and the sampler.
    OD = optimizer_data(optimizer)
    S = sampler(optimizer)

    #Proportion updates
    if error_proportion(OD)>0.9
        println("Lot's of errors - code adjustment of numerical tolerances")
    end
    if taboo_proportion(OD)>0.8
        println("Lot's of taboo fibres - reducing radius of sampler")
        S.transform_matrix = S.transform_matrix*0.9
    end
    if taboo_proportion(OD)<0.2 && improvement_proportion(OD)>0.0
        println("Not enough taboo fibres - we can be bolder in our sampling -increase radius")
        S.transform_matrix = S.transform_matrix*1.1
    end
    if improvement_proportion(OD)>0.2
        println("Lot's of improvement - increase radius")
        S.transform_matrix = S.transform_matrix*1.1
    end
    if improvement_proportion(OD)<0.01 && steps_no_progress(OD)>5
        println("Almost no improvement - decrease radius - maybe we are at a local max")
        S.transform_matrix = S.transform_matrix*0.5
    end
    if improvement_proportion(OD)>0.0
        println("Since there was an improvement, we need to adjust the translator of sampler")
        S.translation = parameters(record_fibre(optimizer))
    end
    return nothing
end


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
    OD.error_proportion = length(error_fibs)/N                           #Setters to be added.
    OD.taboo_proportion = length(taboo_fibs)/N
    OD.improvement_proportion = (length(min_improv)+length(maj_improv))/N

    #Update record fibre/solver fibre

    if OD.improvement_proportion>0.0
        OD.steps_no_progress=0                                           #Setters to be added.
        println("Minor progress")
        if length(maj_improv)>0
            OD.steps_no_major_progress = 0                               #Setters to be added.
            println("Major progress")
        end
        optimizer.record_fibre = new_fibres[best_score_index[1]]         #Setters to be added.
        optimizer.record_objective = weighted_objective(optimizer)(record_fibre(optimizer)[1])
    end
    optimizer.path = push!(path(optimizer), record_fibre(optimizer)[2])
    return optimizer
end



# update_sampler! - the function that updates the sampler. Specifically, this
# functions deals with updating the transform_matrix field of the sampler,
# 1) to make sure the new samples are biased towards the direction of the record,
# 2) to make sure that the proportions (such as taboo and errors) stay within the bounds we define.


import LinearAlgebra: diagm

function update_sampler!(optimizer,new_fibres)
    #getting the data
    vec = record_fibre(optimizer)[2]
    S = sampler(optimizer)
    
    #Fixing the singular values for scaling
    σ1 = 2.0
    σ0 = 1.0
    U = hcat(vec,nullspace(transpose(vec)))
    Σ = diagm(vcat(σ1,[σ0 for i in 2:length(vec)]))
    B = U*Σ*transpose(U)

    #SVD computation
    #=
    param_of_fibres = [new_fibres[i][2] for i in 1:length(new_fibres)]
    mean_param = sum(params for params in param_of_fibres)/length(param_of_fibres)
    A = stack([w_obj(params).*(params-mean_param) for params in param_of_fibres])
    u, s, v = svd(A, full = true)
    Σ = diagm(s)                                    # To obtain the diagonal matrix with the
    n = length(s)                                                # diagonal entries given by s.
    for i in (n+1):N
        Σ = hcat(Σ, [0 for k in 1:n])
    end
    =#

    # S.transform_matrix updated according to the SVD.
    S.transform_matrix = B
    # S.translation needs to be updated, probably to the parameter of the record_fibre.

    #Warnings and modifications to optimizer due to (error, taboo etc.) proportions.
    #updates_from_proportions!(optimizer)
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
    #update_sampler!(optimizer,new_fibres)

    return(optimizer)
end


# The function optimize! updates the input O::Optimizer by applying the function 
# improve! to it until the goal is satisified.

function optimize!(O::Optimizer; max_steps = 100)
    G = goal(O)
    while ((G(O)==false) && (steps_no_major_progress(optimizer_data(O))<max_steps))
            improve!(O)
            println(O)
    end
    return(O)
 end

 #=
function optimize!(O::Optimizer; max_steps = 100, visualize_optimizer = false)
    G = goal(O)
    if visualize_optimizer && (length(parameters(enumerative_problem(O)))==2)
        visualize(enumerative_problem(O))[2]
        while ((G(O)==false) && (steps_no_major_progress(optimizer_data(O))<max_steps))
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
        while ((G(O)==false) && (steps_no_major_progress(optimizer_data(O))<max_steps))
            improve!(O)
            println(O)
        end
    return(O)
    else
        while ((G(O)==false) && (steps_no_major_progress(optimizer_data(O))<max_steps))
            improve!(O)
            println(O)
        end
    return(O)
    end
end
=#


 # optimizer! function without visualize:
 


##

#=

@kwdef mutable struct OptimizerData
    step :: Int = 0                             #Step count 
    parameters_solved :: Int = 0                
    steps_no_progress :: Int = 0
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end


### The mutable struct ScoringScheme
mutable struct ScoringScheme
    objective                       #function that we really want to optimize
    barrier                         #barrier penalty
    barrier_weight :: Float64       #weight of barrier penalty
    taboo
    error_checker                   #function that checks for errors
    goal                            #goal (:: Union{Bool, Nothing}) that stops 
                                        the optimizer, when it achieves the value "true".
    name :: String                  #name of the ScoringScheme
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

=#












