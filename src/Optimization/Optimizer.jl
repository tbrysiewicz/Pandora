
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

    error_checker = x->false           #Run any computed fibres through this boolean-valued function to determine if they are errors
    taboo  = x->0.0                          #Eliminate any fibres for which taboo score decreases - these are not considered valid candidates
                                             #Steps where taboo increases are considered major improvements

    barrier = x->0.0                         #Among remaining candidates, evaluate objective(fibre)-barrier_weight*barrier(fibre)
    objective = x->0.0
    barrier_weight ::Float64  = 0.0    

    goal  = OD -> OD.step>100

    optimizer_data :: OptimizerData


    function Optimizer(EP::EnumerativeProblem, sampler :: Sampler, objective)
        optimizer = new()
        optimizer.EP = EP
        optimizer.solver_fibre = base_fibre(EP)
        optimizer.sampler=sampler
        p = base_parameters(EP)
        
        if is_real(sampler)
            p = Vector{ComplexF64}(real(p))
        end

        optimizer.record_fibre = (EP(p),p)
        optimizer.record_objective = objective(optimizer.record_fibre)

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

