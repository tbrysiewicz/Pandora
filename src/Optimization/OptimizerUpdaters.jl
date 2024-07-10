function scale_sampler_radius!(O::Optimizer,k)
    O.sampling_ellipsoid = k*(O.sampling_ellipsoid)
end

#This function will scale the sampler based on how many samples were taboo
function update_sampler_radius!(O::Optimizer,information::Dict{Any,Any}; desired_interval =[0.7,0.8])
    if O.taboo_score == TrivialScore
        return()
    end
    n_fibres = get(information,"n_fibres","na")
    n_non_taboo = get(information,"n_non_taboo","na")
    if n_fibres!="na" && n_non_taboo!="na"
        non_taboo_proportion = n_non_taboo/n_fibres
        println("Non taboo proportion:",non_taboo_proportion)
        if desired_interval[1]<non_taboo_proportion
            #there are sufficiently many non-taboo
            if non_taboo_proportion<desired_interval[2]
                #non_taboo proportion is within the desired interval - do nothing
            else
                #there are too many non-taboo
                scale_sampler_radius!(O,2+rand(Float64))
                println("                                                                 radius up: ",eigmax(O.sampling_ellipsoid))
            end
        else
            #there are too few non-taboo
            scale_sampler_radius!(O,rand(Float64))
            println("                                                                 radius down: ",eigmax(O.sampling_ellipsoid))
        end
    else
        println("!!!!!!!!!!!!!",n_fibres,"  ",n_non_taboo)
    end
end

function update_solver_fibre!(O::Optimizer)
    println("Updating solver fibre")
    new_parameters = [O.current_fibre[2]+im*randn(Float64,n_parameters(O.EP)) for i in 1:3]
    sols  = solve_over_params(O.EP,new_parameters)
    #println(sols)
    if length(sols)==0
        println("Failed to update solver fibre")
    else
        println("Solver fibre updated")
        O.solver_fibre = sols[1]
    end
end