
@kwdef mutable struct OptimizerData
    step :: Int = 0
    parameters_solved :: Int = 0
    steps_no_progress :: Int = 0
    steps_no_major_progress :: Int = 0
    error_proportion :: Float64 = 0.0
    taboo_proportion :: Float64 = 0.0
    improvement_proportion :: Float64 = 1.0
end

@kwdef mutable struct Optimizer
    EP::EnumerativeProblem
    solver_fibre :: Fibre
    record_fibre :: Fibre
    record_objective :: Any
    path :: Vector{Vector{ComplexF64}}

    sampler :: Sampler                       #used to sample parameter space (frequently updated)

    error_checker  = always_false        #Run any computed fibres through this boolean-valued function to determine if they are errors
    taboo  = x->0.0                          #Eliminate any fibres for which taboo score decreases - these are not considered valid candidates
                                             #Steps where taboo increases are considered major improvements

    barrier = x->0.0                         #Among remaining candidates, evaluate objective(fibre)-barrier_weight*barrier(fibre)
    objective = x->0.0
    barrier_weight ::Float64  = 0.0    

    goal  = OD -> OD.step>100

    optimizer_data :: OptimizerData


    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, objective)
        optimizer = new()
        optimizer.objective = objective
        optimizer.EP = EP
        optimizer.solver_fibre = base_fibre(EP)
        optimizer.sampler=sampler
        p = base_parameters(EP)
        
        if is_real(sampler)
            p = Vector{ComplexF64}(real(p))
        end

        optimizer.error_checker = always_false
        optimizer.taboo = zero
        optimizer.barrier_weight = 0.0

        function more_than_one_hundred_steps(O::Optimizer)
            optimizer_data(O).steps >100
        end

        optimizer.goal = more_than_one_hundred_steps

        optimizer.record_fibre = (EP(p),p)
        optimizer.record_objective = objective(optimizer.record_fibre)
        optimizer.path=[p]

        optimizer.optimizer_data = OptimizerData()

        return(optimizer)
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
    #info regarding radius or last improvement size 
    #info regarding conditioning of the ellipse sampler
    #interpretive information 
    #soft-goal info
end


function improve!(optimizer::Optimizer) 

    OD = optimizer_data(optimizer)
    OD.step = OD.step+1
    OD.steps_no_major_progress = OD.steps_no_major_progress+1
    OD.steps_no_progress = OD.steps_no_progress+1

    #First sample some new points
    S = sampler(optimizer)
    new_samples = sample(S)
    N = n_samples(S)

    #solve for the samples
    OD.parameters_solved = OD.parameters_solved + N

    EP = enumerative_problem(optimizer)
    new_solutions = EP(new_samples)
    new_fibres = [(new_solutions[i],new_samples[i]) for i in 1:N]

    update_record!(optimizer,new_fibres)
    update_settings!(optimizer,new_fibres)

    return(optimizer)
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
    if OD.improvement_proportion<0.01
        println("Almost no improvement - decrease radius - maybe we are at a local max")
        S.transform_matrix = S.transform_matrix*0.5
    end
    if OD.improvement_proportion>0.0
        println("Since there was an improvement, we need to adjust the translator of sampler")
        S.translation = parameters(optimizer.record_fibre)
    end
end

function update_record!(optimizer,new_fibres)
    OD = optimizer_data(optimizer)
    N = length(new_fibres)

    non_error_fibres = filter(x->error_checker(optimizer)(x)==false,new_fibres)

    if isempty(non_error_fibres)
        OD.error_proportion = 1.0
        OD.taboo_proportion = 1.0
        OD.improvement_proportion = 0.0
        return()
    end

    OD.error_proportion = 1-(length(non_error_fibres)/N)

    current_taboo_score = taboo(optimizer)(record_fibre(optimizer))
    new_taboo_scores = map(x->taboo(optimizer)(x[1]),non_error_fibres)

    m = min(new_taboo_scores...)

    non_taboo_fibres = filter(x->taboo(optimizer)(x[1])>=current_taboo_score,non_error_fibres)

    if isempty(non_taboo_fibres)
        OD.taboo_proportion = 1.0
        OD.improvement_proportion = 0.0
        return()
    end

    OD.taboo_proportion = 1-(length(non_taboo_fibres)/N)
    if m<current_taboo_score #This indicates a major improvement - non_taboo_fibres changes, but #non_taboo doesn't
        current_taboo_score = m
        OD.steps_no_major_progress = 0
        println("Major Progress")
    end

    non_taboo_fibres = filter(x->taboo(optimizer)(x[1])>=current_taboo_score,non_error_fibres)

    new_scores = map(objective_function(optimizer),non_taboo_fibres)

    improvement_indices = findall(x->x>record_objective(optimizer),new_scores)
    improvements = non_taboo_fibres[improvement_indices]


    if isempty(improvements)
        OD.improvement_proportion = 0.0
        return()
    end
    println("Minor Progress")

    best_score_index = argmax(new_scores[improvement_indices])
    best_score = new_scores[improvement_indices][best_score_index]
    best_fibre = new_fibres[improvement_indices][best_score_index]
    optimizer.record_fibre = best_fibre
    optimizer.record_objective = best_score
    push!(optimizer.path,parameters(best_fibre))
    OD.steps_no_progress = 0

    OD.improvement_proportion = length(improvements)/N
end