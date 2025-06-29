export 
    ScoringScheme,
    dietmaier_scheme,
    dietmaier_optimizer,
    improve,
    improve!,
    update_optimizer_data!,
    update_optimizer_parameters!,
    optimizer_run,
    Optimizer,
    OptimizerData,
    optimize!,
    optimize_n_real_solutions




"""
    The struct ScoringScheme collects the information required to run optimization code in Pandora.
    Its fields are:
    - `objective`: A function from `Fibre` to some type `T` where elements of `T` can be compared. This function represents what the user truly wants to optimize locally.
    - `taboo`: A function from `Fibre` to `S`, where elements of `S` can be compared.
    - `name`: A `String` that names the scoring scheme.
    The function `improve(EP::EnumerativeProblem, F::Fibre, SS:ScoringScheme; Sampler, n_samples)` 
        will sample `n_samples` new parameters, using `Sampler`, solve for the fibres `FF` of the EnumerativeProblem `EP` for these parameters,
        and return the fibre of `FF` which is "best" according to the `ScoringScheme` `SS`.
    A fibre `F'` is better than a fibre `F` according to the `ScoringScheme` `SS` if:
    - `taboo(SS)(F')` is less than `taboo(SS)(F)`, or if both are equal, then
    - `objective(SS)(F')` is larger than `objective(SS)(F)`.
"""
@kwdef mutable struct ScoringScheme
    objective                               # The function that we really want to be locally optimized.
    taboo            
    name :: String                          # name of the ScoringScheme.
end 



objective(SS:: ScoringScheme) = SS.objective
taboo(SS::ScoringScheme) = SS.taboo
name(SS::ScoringScheme) = SS.name
# Base.show for ScoringScheme
function Base.show(io::IO, SS::ScoringScheme)
    if SS.name == ""
        print(io,"A scoring scheme")
    else
        print(io, name(SS)," scoring scheme")
    end
end


function dietmaier_scheme() :: ScoringScheme
    return(ScoringScheme(objective = S -> -dietmaier(S), taboo = S->-n_real_solutions(S), name = "Dietmaier"))
end

function (SS::ScoringScheme)(S::Vector{Vector{ComplexF64}}) 
    return (taboo(SS)(S), objective(SS)(S))
end


"""
    optimizer_run(EP::EnumerativeProblem, SS::ScoringScheme; sampler::Sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples::Int = 100, solver_fibre::Fibre = base_fibre(EP))

This function will sample `n_samples` new parameters, using `sampler`, solve for the fibres `FF` of the EnumerativeProblem `EP` for these parameters,
    and return a triple containing:
- `fibre_data`: A dictionary mapping each fibre `F` to its score according to
- `best_fibre`: The fibre `F` which is "best" according to the `ScoringScheme` `SS`.
- `best_score`: The score of the `best_fibre` according to the `ScoringScheme` `SS`.
"""
function optimizer_run(EP::EnumerativeProblem,  SS::ScoringScheme; sampler::Sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples::Int = 100, solver_fibre::Fibre = base_fibre(EP))
    # This function will sample `n_samples` new parameters, using `sampler`, solve for the fibres `FF` of the EnumerativeProblem `EP` for these parameters,
    # and return the fibre of `FF` which is "best" according to the `ScoringScheme` `SS`.
    
    new_parameters = sampler(n_samples)
    new_solutions = EP(solver_fibre, new_parameters)
    fibre_data = Dict{Fibre, Any}()
    
    for (sol, p) in zip(new_solutions, new_parameters)
        F = Fibre((sol, p))
        if valid_fibre(EP, F)
            fibre_data[F] = SS(F[1])
        else
            fibre_data[F] = nothing
        end
    end
    
    # Find the best fibre according to the scoring scheme.
    best_fibre = nothing
    best_score = nothing
    for (F, score) in fibre_data
        if score !== nothing
            if best_fibre === nothing || (score[1] < fibre_data[best_fibre][1] || (score[1] == fibre_data[best_fibre][1] && score[2] > fibre_data[best_fibre][2])) #
                best_fibre = F
                best_score = score
            end
        end
    end

    @error "No valid fibres found." if best_fibre === nothing
    return(fibre_data, best_fibre, best_score)
end


#The struct OptimizerData : keeps track of the Data involved while using the optimizer.
@kwdef mutable struct OptimizerData
    record_fibre :: Fibre = nothing 
    record_score :: Any = nothing 

    step :: Int = 0                             
    parameters_solved :: Int = 0                
    steps_no_objective_progress :: Int = 0                
    steps_no_taboo_progress :: Int = 0
    
    last_error_proportion :: Float64 = 0.0
    last_taboo_proportion :: Float64 = 0.0
    last_improvement_proportion :: Float64 = 1.0
    
    path :: Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}() 
end



@kwdef mutable struct Optimizer
    EP::EnumerativeProblem
    solver_fibre :: Fibre

    sampler :: TransformedSampler                       #used to sample parameter space (frequently updated)
    scoring_scheme :: ScoringScheme
    optimizer_data :: OptimizerData
    goal :: Function
end


function optimizer_run(O::Optimizer)
    return optimizer_run(O.EP, O.scoring_scheme; sampler = O.sampler, n_samples = 100, solver_fibre = O.solver_fibre)
end

# Base.show for Optimizer

function radius(O::Optimizer)
    det(O.sampler.affine_transformation.transform_matrix)^(1/n_parameters(O.EP))
end

function Base.show(io::IO, optimizer::Optimizer)
    tenspaces="          "
    println(io,"Optimizer for an enumerative problem of degree ",degree(optimizer.EP))
    println("--------------------------------------------------------------------------")
    od = optimizer.optimizer_data
    println("Progress:                       Current score:  ",od.record_score)
    println("                                       Step #:  ",od.step)
    println("                      Total parameters solved:  ",od.parameters_solved)
    println("          Steps since last objective progress:  ",od.steps_no_objective_progress)
    println("                    since last taboo progress:  ",od.steps_no_taboo_progress)
    println("                        Last error proportion:  ",od.last_error_proportion)
    println("                        Last taboo proportion:  ",od.last_taboo_proportion)
    println("                  Last improvement proportion:  ",od.last_improvement_proportion)
    println("                                       Radius: ", radius(optimizer))

end

function dietmaier_optimizer(EP::EnumerativeProblem)
    SS = dietmaier_scheme()
    sampler = UniformSampler{Float64}(n_parameters(EP))
    p = sampler(1)[1]
    start_fibre = Fibre((EP(p), p))
    record_score = SS(start_fibre[1])
    OD = OptimizerData(record_fibre = start_fibre, record_score = record_score, path = [start_fibre[2]])
    Optimizer(EP = EP, solver_fibre = base_fibre(EP), sampler = sampler, scoring_scheme = SS, optimizer_data = OD, goal = OptD -> OptD.record_score[1] == -degree(EP))
end

function success(O::Optimizer)
    # Check if the optimizer has reached its goal.
    return O.goal(O.optimizer_data)
end

function improve!(O::Optimizer; n_samples::Int = 100)

    (fibre_data, best_fibre, best_score) = optimizer_run(O.EP, O.scoring_scheme; sampler = O.sampler, n_samples = n_samples, solver_fibre = O.solver_fibre)

    update_optimizer_data!(O, fibre_data, best_fibre, best_score)
    update_optimizer_parameters!(O)
    
    return(O)
end

function optimize!(O::Optimizer; n_samples::Int = 100, max_steps::Int = 1000)
    # This function will run the optimizer until it reaches its goal or until it has taken `max_steps` steps.
    for step in 1:max_steps
        if success(O)
            println("Optimizer has reached its goal after ", step-1, " steps.")
            return O
        end
        improve!(O; n_samples = n_samples)
        update_optimizer_parameters!(O; optimistic = true)
    end
    println("Optimizer did not reach its goal after ", max_steps, " steps.")
    return O
end

function update_optimizer_parameters!(O::Optimizer; optimistic = true)
    sampler = O.sampler
    OD = O.optimizer_data
    if optimistic
        # If we are being optimistic, the step in the parameter space we just took, we will double
        last_parameter = OD.path[end-1]
        current_parameter = OD.path[end]
        t = current_parameter+(current_parameter - last_parameter)
        if is_real(t) 
            t = real(t)
        end
        set_translation!(sampler, t)
    else
        # If we are being pessimistic, we will just take a step in the direction of the last step.
        last_parameter = OD.path[end]
        set_translation!(sampler, last_parameter)
    end

    #Now we update the transformation matrix either by scaling or PCA 
    #TODO: Implement PCA
    # For now, we will just scale
    if OD.last_error_proportion > 0.9 
        # Change solver fibre
        O.solver_fibre = OD.EP(OD.record_fibre[2]+ 0.1*randn(ComplexF64, n_parameters(O.EP)))
        println("Changing solver fibre since so many errors occurred.")
    end
    if OD.last_taboo_proportion > 0.8
        # Scale down the sampler by 0.9
        scale!(sampler, 0.9)
        println("Scaling down sampler since so many taboo fibres were found: ", radius(O))
    end
    if OD.last_taboo_proportion < 0.2 && OD.last_improvement_proportion > 0.0 
        # Scale up the sampler by 1.1
        scale!(sampler, 1.1)
        println("Scaling up sampler since taboo fibres were not found and some improvement was made: ", radius(O))
    end 
    if OD.last_improvement_proportion > 0.4
        # If we are making good progress, we will scale up the sampler by 1.2
        scale!(sampler, 1.2)
        println("Scaling up sampler since good improvement was made: ", radius(O))
    end
    if OD.last_improvement_proportion == 0.0 && OD.steps_no_objective_progress > 10
        # If we are not making any progress, we will scale down the sampler by 0.8
        scale!(sampler, 0.5)
        println("Scaling down sampler since no improvement was made for a long time. Probably at a local optimum: ", radius(O))
    end
    print(".")
end

function update_optimizer_data!(O::Optimizer, fibre_data::Dict{Fibre, Any}, best_fibre::Fibre, best_score::Any)
    OD = O.optimizer_data
    # Update the optimizer data with the new best fibre and score.
    old_record_score = OD.record_score
    if OD.record_score[1] > best_score[1] # If the new best score has improved (decreased) taboo. 
        OD.steps_no_taboo_progress = 0
        OD.record_score = best_score
        OD.record_fibre = best_fibre
        println("Taboo improvement: ", best_score)
    else
        OD.steps_no_taboo_progress += 1
        if OD.record_score[1] == best_score[1] && OD.record_score[2] < best_score[2] # If the new best score has improved (increased) objective.
            OD.steps_no_objective_progress = 0
            OD.record_score = best_score
            OD.record_fibre = best_fibre
            println("Objective improvement: ", best_score)
        else
            OD.steps_no_objective_progress += 1
        end
    end

    # Update the step count and parameters solved.
    OD.step += 1
    N = length(fibre_data)
    OD.parameters_solved += N

    # Update the path with the new best fibre's parameters.
    push!(OD.path, OD.record_fibre[2])

    # Calculate error and taboo proportions.
    error_fibre_count = count(x -> valid_fibre(O.EP,x)==false, keys(fibre_data))
    taboo_fibre_count = count(x -> valid_fibre(O.EP,x) && taboo(O.scoring_scheme)(x[1]) > old_record_score[1], keys(fibre_data))
    improvement_fibre_count = count(x -> valid_fibre(O.EP,x) && taboo(O.scoring_scheme)(x[1]) < old_record_score[1] && objective(O.scoring_scheme)(x[1])>old_record_score[2], keys(fibre_data))

    OD.last_error_proportion = error_fibre_count / N
    OD.last_taboo_proportion = taboo_fibre_count / N
    OD.last_improvement_proportion = improvement_fibre_count / N
end

function optimize_n_real_solutions(EP::EnumerativeProblem; n_samples::Int = 2*n_parameters(EP), max_steps::Int = 100)
    # This function will run the optimizer until it reaches its goal of finding a fibre with n_real_solutions equal to degree(EP).
    O = dietmaier_optimizer(EP)
    return optimize!(O; n_samples = n_samples, max_steps = max_steps)
end