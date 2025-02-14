
@kwdef mutable struct OptimizerData
    step :: Int = 0
    parameters_solved :: Int = 0
    steps_no_progress :: Int = 0
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end


mutable struct ScoringScheme
    objective #function we really want to optimize
    barrier #barrier penalty
    barrier_weight :: Float64 #weight of barrier penalty
    taboo
    error_checker
    goal
    name :: String
end


function weighted_objective(SS::ScoringScheme)
    x->(1-SS.barrier_weight).*SS.objective(x).+(SS.barrier_weight).*SS.barrier(x)
end



function Base.show(io::IO, SS::ScoringScheme)
    if SS.name==""
        print(io,"A scoring scheme with barrier weight ",SS.barrier_weight,".")
    else
        print(io,"The ",SS.name," scoring scheme with barrier weight ",SS.barrier_weight, ".")
    end
end


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


function ScoringScheme(objective; 
    barrier = zero_function, barrier_weight = 0.0, 
    taboo = zero_function, error_checker = false_function,
    goal = nothing, name = "")
    function more_than_one_hundred_steps(O::Optimizer)
        O.optimizer_data.step >100
    end
    if goal == nothing
        ScoringScheme(objective,barrier,barrier_weight,taboo,error_checker, more_than_one_hundred_steps,name)
    else
        ScoringScheme(objective,barrier,barrier_weight,taboo,error_checker, goal,name)
    end
end

#This is the function which increments all relevant step datum in an optimizer data
#  prior to running an improvement step
function increment!(OD::OptimizerData, N::Int) 
    OD.step = OD.step+1
    OD.steps_no_major_progress = OD.steps_no_major_progress+1
    OD.steps_no_progress = OD.steps_no_progress+1

    OD.parameters_solved = OD.parameters_solved + N
end



function improve!(optimizer::Optimizer) 

    #Extract relevant names
    OD = optimizer_data(optimizer)
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
    update_settings!(optimizer,new_fibres)

    return(optimizer)
end


function optimize!(O::Optimizer; max_steps = 100)
    G = goal(O)
    while G(O)==false && O.optimizer_data.steps_no_major_progress<max_steps
        improve!(O)
        println(O)
    end
    return(O)
end


function update_settings!(optimizer,new_fibres)
    OD = optimizer_data(optimizer)
    S = sampler(optimizer)
    if OD.error_proportion>0.9
        println("Lot's of errors - code adjustment of numerical tolerances")
    end
    if OD.taboo_proportion>0.8
        println("Lot's of taboo fibres - reducing radius of sampler")
        S.transform_matrix = S.transform_matrix*0.9
    end
    if OD.taboo_proportion<0.2 && OD.improvement_proportion>0.0
        println("Not enough taboo fibres - we can be bolder in our sampling -increase radius")
        S.transform_matrix = S.transform_matrix*1.1
    end
    if OD.improvement_proportion>0.2
        println("Lot's of improvement - increase radius")
        S.transform_matrix = S.transform_matrix*1.1
    end
    if OD.improvement_proportion<0.01 && OD.steps_no_progress>5
        println("Almost no improvement - decrease radius - maybe we are at a local max")
        S.transform_matrix = S.transform_matrix*0.5
    end
    if OD.improvement_proportion>0.0
        println("Since there was an improvement, we need to adjust the translator of sampler")
        S.translation = parameters(optimizer.record_fibre)
    end
end

function partition_optimizer_samples(optimizer:: Optimizer, new_fibres :: Vector{Fibre}) 
    n = length(new_fibres)

    error_fibres = findall(i->error_checker(optimizer)(new_fibres[i]), 1:n)
    
    taboo_function = taboo(optimizer)
    cts = taboo_function(record_fibre(optimizer)) ##Current taboo score
    taboo_score_list = [in(i,error_fibres) ? nothing : taboo_function(new_fibres[i]) for i in 1:n]
    
    taboo_fibres = findall(i->in(i,error_fibres)==false && taboo_score_list[i]>cts, 1:n)

    w_obj = weighted_objective(optimizer)
    wo_record = w_obj(record_fibre(optimizer))
    wo_score_list = [in(i,union(taboo_fibres,error_fibres)) ? nothing : w_obj(new_fibres[i]) for i in 1:n]

    improvement_fibres = findall(i->in(i,union(taboo_fibres,error_fibres))==false && 
                                    wo_score_list[i]>wo_record,1:n)
    major_improvement_fibres = findall(i->in(i,improvement_fibres) && taboo_score_list[i]<cts,1:n)
    minor_improvement_fibres = findall(i->in(i,improvement_fibres) && taboo_score_list[i]==cts, 1:n)
    no_improvement_fibres = findall(i->in(i,union(taboo_fibres,
                                                  error_fibres,
                                                  improvement_fibres))==false, 1:n)
    
    best_score_index = []
    if length(improvement_fibres)>0
        true_scores = filter(x->x!=nothing,wo_score_list)
        true_score_indices = findall(x->x!=nothing,wo_score_list)
        @assert length(true_score_indices)>0
        b = argmax(true_scores)
        best_score_index = [true_score_indices[b]]
    end
    result = (error_fibres,taboo_fibres,no_improvement_fibres,minor_improvement_fibres,major_improvement_fibres,best_score_index)


    return(result)

end

function update_optimizer!(optimizer :: Optimizer, new_fibres :: Vector{Fibre})
    OD = optimizer_data(optimizer)
    N = length(new_fibres)

    #Make a default increment update to optimizer data (steps += 1, params solved +=N)
    increment!(OD,N)

    (error_fibs, taboo_fibs, non_improv, min_improv, maj_improv, best_score_index) = partition_optimizer_samples(optimizer,new_fibres)
    
    #Updating OD
    OD.error_proportion = length(error_fibs)/N
    OD.taboo_proportion = length(taboo_fibs)/N
    OD.improvement_proportion = (length(min_improv)+length(maj_improv))/N




    #Update record fibre/solver fibre

    if OD.improvement_proportion>0.0
        OD.steps_no_progress=0
        println("Minor progress")
        if length(maj_improv)>0
            OD.steps_no_major_progress = 0
            println("Major progress")
        end
        optimizer.record_fibre = new_fibres[best_score_index[1]]
        optimizer.record_objective = weighted_objective(optimizer)(optimizer.record_fibre)
    end


end