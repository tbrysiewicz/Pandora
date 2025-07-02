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
    maximize_n_real_solutions,
    anti_dietmaier_scheme,
    anti_dietmaier_optimizer,
    change_scoring_scheme!,
    record_fibre,
    record_score,
    obtain_k_reals_scheme,
    obtain_k_optimizer,
    record_parameters,
    record_solutions




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



"""
    objective(SS::ScoringScheme)

This function returns the objective function of the ScoringScheme `SS`.
"""
objective(SS:: ScoringScheme) = SS.objective
"""
    taboo(SS::ScoringScheme)
This function returns the taboo function of the ScoringScheme `SS`.
"""
taboo(SS::ScoringScheme) = SS.taboo
"""
    name(SS::ScoringScheme)

This function returns the name of the ScoringScheme `SS`.
"""
name(SS::ScoringScheme) = SS.name
# Base.show for ScoringScheme
function Base.show(io::IO, SS::ScoringScheme)
    if SS.name == ""
        print(io,"A scoring scheme")
    else
        print(io, name(SS)," scoring scheme")
    end
end

"""
    dietmaier_scheme(d)

This function returns a `ScoringScheme` that is used to minimize 
    the smallest imaginary part of a non-real solution. Doing so
    is a heuristic for "pushing" non-real solutions together to 
    produce more real solutions.
The taboo function of this scheme is the difference between
    the number of real solutions and `d`. By default, `k` is the degree
    of the EnumerativeProblem, but it can be set to any integer.
"""
function dietmaier_scheme(d) :: ScoringScheme
    return(ScoringScheme(objective = S -> -dietmaier(S), taboo = S->abs(d-n_real_solutions(S)), name = "Dietmaier"))
end

"""
    anti_dietmaier_scheme(d::Int)

This function returns a `ScoringScheme` that is used to minimize the distance between real solutions.
    Doing so is a heuristic for "pushing" real solutions together to produce 
    more non-real solutions.
The taboo function of this scheme is the the number of real solutions.
"""
function anti_dietmaier_scheme(d::Int) :: ScoringScheme
    return(ScoringScheme(objective = S -> -anti_dietmaier(S), taboo = S->abs(n_real_solutions(S)), name = "Anti Dietmaier"))
end

function obtain_k_reals_scheme(d::Int,k::Int) :: ScoringScheme
    function obtain_k_objective(S::Vector{Vector{ComplexF64}})
        nr = n_real_solutions(S)
        if nr < k 
            return -dietmaier(S)
        else
            return -anti_dietmaier(S)
        end
    end
    function obtain_k_taboo(S::Vector{Vector{ComplexF64}})
        nr = n_real_solutions(S)
        return abs(nr - k)
    end
    return ScoringScheme(objective = obtain_k_objective, 
                         taboo = obtain_k_taboo, 
                         name = "Obtain $k real solutions")
end


function obtain_k_optimizer(EP::EnumerativeProblem, k::Int)
    SS = obtain_k_reals_scheme(degree(EP), k)
    sampler = UniformSampler{Float64}(n_parameters(EP))
    p = sampler(1)[1]
    start_fibre = Fibre((EP(p), p))
    record_score = SS(start_fibre[1])
    OD = OptimizerData(record_fibre = start_fibre, 
                       record_score = record_score, 
                       path = [start_fibre[2]])
    Optimizer(EP = EP, 
              solver_fibre = base_fibre(EP), 
              sampler = sampler, 
              scoring_scheme = SS, 
              optimizer_data = OD, 
              goal = OptD -> OptD.record_score[1] == 0.0)
end

"""
    (SS::ScoringScheme)(S::Vector{Vector{ComplexF64}})

    Evaluates the scoring scheme `SS` on a solution set `S`.
    Such an evaluation returns a pair containing the taboo score
    and the objective score. 
"""
function (SS::ScoringScheme)(S::Vector{Vector{ComplexF64}}) 
    return (taboo(SS)(S), objective(SS)(S))
end


"""
    optimizer_run(EP::EnumerativeProblem, SS::ScoringScheme; sampler::Sampler = UniformSampler{Float64}(n_parameters(EP)), n_samples::Int = 100, solver_fibre::Fibre = base_fibre(EP))
    optimizer_run(O::Optimizer)

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
    if best_fibre === nothing
        @error "No valid fibres found in the sampled data."
    end
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


"""
    Optimizer

    This is the main optimization struct in Pandora.jl for EnumerativeProblems.
    It contains the following fields:
    - `EP`: The EnumerativeProblem to be solved in order to evaluate the objective function.
    - `solver_fibre`: The fibre used to solve the EnumerativeProblem - this will be adjusted
        automatically to speed up the optimization process.
    - `sampler`: A `TransformedSampler` used to sample the parameter space.
    - `scoring_scheme`: A `ScoringScheme` that defines how to score the solution sets in the fibres.
    - `optimizer_data`: An `OptimizerData` struct that keeps track of the optimization process.
    - `goal`: A function from `OptimizerData` to `Bool` that defines the goal of the optimization process. 
"""
@kwdef mutable struct Optimizer
    EP::EnumerativeProblem
    solver_fibre :: Fibre

    sampler :: TransformedSampler                       #used to sample parameter space (frequently updated)
    scoring_scheme :: ScoringScheme
    optimizer_data :: OptimizerData
    goal :: Function
end

record_fibre(O::Optimizer) = O.optimizer_data.record_fibre
record_score(O::Optimizer) = O.optimizer_data.record_score
record_parameters(O::Optimizer) = O.optimizer_data.record_fibre[2]
record_solutions(O::Optimizer) = O.optimizer_data.record_fibre[1]

function optimizer_run(O::Optimizer)
    return optimizer_run(O.EP, O.scoring_scheme; sampler = O.sampler, n_samples = 100, solver_fibre = O.solver_fibre)
end


"""
    radius(O::Optimizer)

    This function calculates a number associated to a sampler which we call
        the `radius`. It is defined as the determinant of the affine transformation matrix
        of the sampler raised to the power of `1/n_parameters(O.EP)`.
        When the affine transformation matrix is a scalar matrix c*I, the radius is c. 
    
    The radius is printed to the user throughout the optimization process to give
       an idea of how the optimization is progressing. Small radius means that the sampler
       was forced to shrink so that fewer samples were taboo. 
"""
function radius(O::Optimizer)
    det(O.sampler.affine_transformation.transform_matrix)^(1/n_parameters(O.EP))
end

function Base.show(io::IO, optimizer::Optimizer)
    tenspaces="          "
    println(io,"Optimizer for an enumerative problem of degree ",degree(optimizer.EP))
    if VERBOSE[]==true
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
end

function dietmaier_optimizer(EP::EnumerativeProblem)
    SS = dietmaier_scheme(degree(EP))
    sampler = UniformSampler{Float64}(n_parameters(EP))
    p = sampler(1)[1]
    start_fibre = Fibre((EP(p), p))
    record_score = SS(start_fibre[1])
    OD = OptimizerData(record_fibre = start_fibre, 
                       record_score = record_score, 
                       path = [start_fibre[2]])
    Optimizer(EP = EP, 
              solver_fibre = base_fibre(EP), 
              sampler = sampler, 
              scoring_scheme = SS, 
              optimizer_data = OD, 
              goal = OptD -> OptD.record_score[1] == 0)
end

function anti_dietmaier_optimizer(EP::EnumerativeProblem)
    SS = anti_dietmaier_scheme(degree(EP))
    sampler = UniformSampler{Float64}(n_parameters(EP))
    p = sampler(1)[1]
    start_fibre = Fibre((EP(p), p))
    record_score = SS(start_fibre[1])
    OD = OptimizerData(
        record_fibre = start_fibre,
        record_score = record_score,
        path = [start_fibre[2]]
    )
    Optimizer(
        EP = EP,
        solver_fibre = base_fibre(EP),
        sampler = sampler,
        scoring_scheme = SS,
        optimizer_data = OD,
        goal = OptD -> OptD.record_score[1] == 0.0
    )
end

"""
    success(O::Optimizer)

This function checks if the optimizer has reached its goal.
"""
function success(O::Optimizer)
    # Check if the optimizer has reached its goal.
    return O.goal(O.optimizer_data)
end

"""
    improve!(O::Optimizer; n_samples::Int = 100)

    This is the main function for the optimization process.
    It calls `optimizer_run` to sample new parameters, solve for the fibres, and
       evaluate their scores. 
    It then updates the `OptimizerData` with the relevant data of the `optimizer_run`,
       For example, the proportion of errors, taboo fibres, and improvements, the new record fibres and scores, etc.
    Finally, it updates the parameters of the optimizer based on the results of the optimization step.
       For example, it may scale the sampler or change the solver fibre based on the results.
"""
function improve!(O::Optimizer; n_samples::Int = 100)

    (fibre_data, best_fibre, best_score) = optimizer_run(O.EP, O.scoring_scheme; sampler = O.sampler, n_samples = n_samples, solver_fibre = O.solver_fibre)

    update_optimizer_data!(O, fibre_data, best_fibre, best_score)
    update_optimizer_parameters!(O)
    
    return(O)
end

"""
    optimize!(O::Optimizer; n_samples::Int = 100, max_steps::Int = 100)

This function will run the optimizer until it reaches its goal or until it has taken `max_steps` steps.
"""
function optimize!(O::Optimizer; n_samples::Int = 100, max_steps::Int = 100)
    # This function will run the optimizer until it reaches its goal or until it has taken `max_steps` steps.
    for step in 1:max_steps
        if success(O)
            println("Optimizer has reached its goal after ", step-1, " steps.")
            return O
        end
        improve!(O; n_samples = n_samples)
    end
    println("Optimizer did not reach its goal after ", max_steps, " steps.")
    return O
end

#TODO: Refactor this function to run through a collection of functions 
#      which interpret `OptimizerData`. For example, function is_at_local_optimum(OD::OptimizerData) which checks if the last improvement proportion is 0.0 and the steps no objective progress is greater than 10.
#      would return true if it is believed the optimizer is at a local optimum.
#      Such functions could be combined into a struct which comes to form 
#      an attitude about the optimizer's current state. These attitudes could 
#      then be used to update the optimizer parameters. 
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
    if OD.last_error_proportion > 0.2
        # Change solver fibre
        new_parameter = OD.record_fibre[2] + 0.1*randn(ComplexF64, n_parameters(O.EP))
        O.solver_fibre = (O.EP(new_parameter),new_parameter)
        @vprintln("Changing solver fibre since so many errors occurred.")
    end
    if OD.last_taboo_proportion > 0.8
        # Scale down the sampler by 0.9
        scale!(sampler, 0.9)
        @vprintln("Scaling down sampler since so many taboo fibres were found: ", radius(O))
    end
    if OD.last_taboo_proportion < 0.2 && OD.last_improvement_proportion > 0.0 
        # Scale up the sampler by 1.1
        scale!(sampler, 1.1)
        @vprintln("Scaling up sampler since taboo fibres were not found and some improvement was made: ", radius(O))
    end 
    if OD.last_improvement_proportion > 0.4
        # If we are making good progress, we will scale up the sampler by 1.2
        scale!(sampler, 1.2)
        @vprintln("Scaling up sampler since good improvement was made: ", radius(O))
    end
    if OD.last_improvement_proportion == 0.0 && OD.steps_no_objective_progress > 10
        # If we are not making any progress, we will scale down the sampler by 0.8
        scale!(sampler, 0.5)
        @vprintln("Scaling down sampler since no improvement was made for a long time. Probably at a local optimum: ", radius(O))
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
        @vprintln("Taboo improvement")
    else
        OD.steps_no_taboo_progress += 1
        if OD.record_score[1] == best_score[1] && OD.record_score[2] < best_score[2] # If the new best score has improved (increased) objective.
            OD.steps_no_objective_progress = 0
            OD.record_score = best_score
            OD.record_fibre = best_fibre
            @vprintln("Objective improvement")
        else
            OD.steps_no_objective_progress += 1
        end
    end

    @vprintln("Current score (Taboo, Objective): ", record_score(O))

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

function maximize_n_real_solutions(EP::EnumerativeProblem; n_samples::Int = 2*n_parameters(EP), max_steps::Int = 100)
    # This function will run the optimizer until it reaches its goal of finding a fibre with n_real_solutions equal to degree(EP).
    O = dietmaier_optimizer(EP)
    return optimize!(O; n_samples = n_samples, max_steps = max_steps)
end


function find_k_real_solutions(EP::EnumerativeProblem, k::Int; n_samples::Int = 2*n_parameters(EP), max_steps::Int = 100)
    # This function will run the optimizer until it reaches its goal of finding a fibre with n_real_solutions equal to k.
    O = obtain_k_optimizer(EP, k)
    return optimize!(O; n_samples = n_samples, max_steps = max_steps)
end


"""
    change_scoring_scheme!(O::Optimizer, SS::ScoringScheme)

"""
function change_scoring_scheme!(O::Optimizer, SS::ScoringScheme)
    # This function changes the scoring scheme of the optimizer.
    O.scoring_scheme = SS
    O.optimizer_data.record_score = SS(O.optimizer_data.record_fibre[1])
    @vprintln("Changed scoring scheme to ", name(SS))
end


