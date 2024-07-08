export
    optimize_real,
    improve!,
    optimize,
    initialize_optimizer

const Fibre = Union{Tuple{Result,Vector{ComplexF64}},Tuple{Result,Vector{Float64}}}
const RealFibre = Tuple{Result,Vector{Float64}}

function is_score(f::Function)
    RT = Base.return_types(f)
    for r in RT
        if hasmethod(isless,Tuple{r,r})==false 
            return(false)
        end
    end
    return(true)
end


function non_real_min_norm(R::Result)
    record = nothing                #Here, record ans record_fibre are just variable names, not parts of the struct OptimizerData.
    for r in R
        if HomotopyContinuation.is_real(r)==false
            if record!= nothing
                nm = norm(imag(r.solution))
                if nm<record
                    record=nm
                end
            else
                record = norm(imag(r.solution))
            end
        end
    end
    if record==nothing
        return(0.0)
    else
        return(record)
    end
end


function real_min_dist(R::Result)
    S = HomotopyContinuation.real_solutions(R)
    if length(S)==0
      return(0.0)
    end
    record=Inf
    n = length(S)
    for i in 1:n
      for j in 1:n
        if i!=j
           s1=S[i]
           s2=S[j]
           nm = norm(s1-s2)
           if nm<record
              record=nm
           end
        end
       end
    end
    return(record)
    end
    
    



function cts_real(RP::Fibre)
    return((nreal(RP[1]),-non_real_min_norm(RP[1])))
end

function real_taboo(RP::Fibre)
    return(-length(HomotopyContinuation.real_solutions(RP[1])))
end

function real_barrier(RP::Fibre)
    return(0,1/real_min_dist(RP[1]))
end


function max_score(fibres::Vector{RealFibre},objective_score::Function,barrier_score::Function,barrier_weight::Float64)
    #initialize record holder
    record_fibre = fibres[1]
    f = x->objective_score(x).-barrier_weight.*barrier_score(x)
    record = f(record_fibre)
    for i in 2:length(fibres)
        newscore = f(fibres[i])
        if newscore>record
            record_fibre = fibres[i]
            record = newscore
        end
    end
    return((record,record_fibre))
end

#The Optimizer struct holds all of the information gathered during an optimization session

mutable struct Optimizer
    EP::EnumerativeProblem

    #If history is enabled, this tracks the history of record fibres and their record scores
    history::Vector{Tuple{Fibre,Any}}

    #current_fibre is the current fibre
    current_fibre::Fibre
    current_score::Any

    #record_fibre holds the fibre which has the best score - it may not be the last fibre in the history
    record_fibre::Fibre
    record_score::Any

    #solver_fibre is the fibre the optimizer is using as a start system
        #It cannot be "current fibre" since current fibre is often over real parameters
        # so instead it is best to use a complex parameter `nearby`
    solver_fibre::Fibre

    #The following scores, and the barrier_weight govern the local search
    objective_score::Function  #This is the score trying to be optimized
                            # any fibre whose score is higher than the current record is a candidate for the next fibre
    taboo_score::Function      #This is a score function, where any fibre's score which is smaller than the current fibre's score is considered taboo
                            # one class of taboo scores is  1/0-valued score where 0 is considered taboo
    barrier_score::Function    #This is a barrier score which penalizes fibres which evaluate to large values
    barrier_weight::Float64 #This is a weight w so that the true objective function is (1-w)objective-(w)barrier

    #sampling_ellipsoid is a matrix which represents a linear transformation
    #  the sampling is the image of the normal distribution N(0,1) under this transformation
    sampling_ellipsoid::Matrix{Float64}

    goal::Union{Function,Nothing} #This is the goal function. Optimization will stop if goal is reached
end

function Base.show(io::IO, O::Optimizer)
    println("---------------------------------------------------------------------------")
    println("Optimizer for an enumerative problem with ",degree(O.EP), " many solutions.")
    println("  Current barrier-weighted objective function: ",O.current_score)
    println("   Record barrier-weighted objective function: ", O.record_score)
    goal_passed = O.goal(O)
    if typeof(O.goal)==Nothing
        println("  Goal:  n/a")
    else
        if  goal_passed==true
            println("  Goal: Reached")
        else
            println("  Goal: Not Reached")
        end
    end
    println("---------------------------------------------------------------------------")
end


#(history, current_fibre, current_score, record_fibre, record_score, solver_fibre, objective_score, taboo_score, barrier_score, barrier_weight, sampling_ellipsoid)
function initialize_optimizer(EP::EnumerativeProblem,SC::Function;
                                taboo_score=x->0,
                                barrier_score=x->0,
                                barrier_weight=0.1,
                                goal = nothing)
    history = Vector{Tuple{Fibre,Any}}([])
    !(is_populated(EP)) && populate_base_fibre(EP)
    solver_fibre = base_fibre(EP)
    objective_score = SC
    S = solve_over_params(EP,[randn(Float64,n_parameters(EP)) for i in 1:10])
    if length(S)==0 #if solving over real parameters failed, return nothing
        return(nothing)
    end
    current_fibre = S[1]
    record_fibre = current_fibre
    current_score = objective_score(current_fibre)
    record_score = current_score
    sampling_ellipsoid = Matrix{Float64}(LinearAlgebra.I,n_parameters(EP),n_parameters(EP))
    push!(history,(record_fibre,record_score))
    data = (EP, history, current_fibre, current_score, record_fibre, record_score, solver_fibre, SC, taboo_score, barrier_score, barrier_weight, sampling_ellipsoid,goal)
    Optimizer(data...)
end

function non_taboo_fibres(O::Optimizer,F::Vector{RealFibre})
    taboo_level = current_taboo_level(O)
    filter(x->taboo_level>=O.taboo_score(x),F)
end

function current_taboo_level(O::Optimizer)
    O.taboo_score(O.current_fibre)
end

#improve will
#  1) sample from T(N(0,1)) where T is given by the sampling ellipsoid
#  2) solve for fibres over those samples
#  3) discard any taboo samples
#  4) evaluate the weighted objective function (penalized by barrier weight) on remaining fibres
#  5) replace current fibre with the record-holder of step (4)
#  6) update history
#  7) return the data from steps 2,3,4,5 for the meta_strategy to adjust the optimizer settings

function improve!(O::Optimizer; n_samples=nothing)
    improvement_info = Dict{Any,Any}()
    improvement_info["n_samples"] = n_samples

    status = ""
    #1)
    if n_samples==nothing
        n_samples = floor(3*n_parameters(O.EP))
    end
    search_vectors = map(x->O.sampling_ellipsoid*x,[randn(Float64,n_parameters(O.EP)) for i in 1:ceil(n_samples/2)])
    search_vectors = vcat(search_vectors,map(x->-x,search_vectors))
    samples = [b+O.current_fibre[2] for b in search_vectors]
    #2)
    sols = solve_over_params(O.EP,samples; start_fibre = O.solver_fibre)
    improvement_info["fibres"] = sols
    improvement_info["n_fibres"] = length(sols)
    println("Number of samples: ",n_samples)
    println("Number of successful tracks: ",length(sols))

    taboo_level = current_taboo_level(O)
    #3)
    non_taboo = non_taboo_fibres(O,sols)
    improvement_info["n_non_taboo"] = length(non_taboo)
    improvement_info["non_taboo_fibres"] = non_taboo
    if length(non_taboo)>0
        #4)
        (record,record_fibre) = max_score(non_taboo,O.objective_score,O.barrier_score,O.barrier_weight)
        #5)
        if record>O.current_score
            status = "Improved Current Score"
            O.current_score = record
            O.current_fibre = record_fibre
            #6)
            push!(O.history,(record_fibre,record))
            if record>O.record_score
                status = "Improved Record Score"
                O.record_score = record
                O.record_fibre = record_fibre
            end
        else
            status = "No Improvement"
        end
    else
        status = "All Are Taboo"
    end
    println("Current Score: ",O.current_score)
    return(improvement_info)
end

function optimize!(O::Optimizer; n_trials = 10)
    trials = 0
    while trials<n_trials && O.goal(O)==false
        trials = trials+1
        println("-----------------------------------------------Trial: ",trials)
        information = improve!(O)
        println(information["n_fibres"])
        println(information["n_non_taboo"])
        update_sampler_radius!(O,information)
        if trials%10 == 0
            update_solver_fibre!(O)
        end
    end
    return(O)
end

function all_real_goal(O::Optimizer)
    degree(O.EP) == length(HomotopyContinuation.real_solutions(O.record_fibre[1]))
end

function optimize_real(EP::EnumerativeProblem; n_trials = Inf)
    O = initialize_optimizer(EP,cts_real;taboo_score = real_taboo, barrier_score = real_barrier, goal = all_real_goal);
    trials = 0
    O = optimize!(O; n_trials = n_trials)
    #=
    while O.record_score[1]<degree(EP)
        trials = trials + 1
        (s,nt,stat) = improve!(O)
        if stat == "Improved Record Score"
            trials = 0
        end
        if nt==0
            println("Shrinking Sampler")
            O.sampling_ellipsoid=O.sampling_ellipsoid/2
        elseif nt==s
            println("Growing sampler")
            O.sampling_ellipsoid=O.sampling_ellipsoid*2
        end
        println("Steps:", trials)
        if trials>10
            O.sampling_ellipsoid=O.sampling_ellipsoid/100
        end
        println(O)
    end
    =#
end