function initialize_Subdivision(EP::EnumerativeProblem, xlims::Vector, ylims::Vector, fibre_function = x->HomotopyContinuation.nreal(x[1]), resolution = 1000)
	xlength = xlims[2] - xlims[1]
	ylength = ylims[2] - ylims[1]

	if xlength == ylength
		num_x_divs = Int(floor(sqrt(resolution)))
		num_y_divs = Int(floor(sqrt(resolution)))
	else
		xlength_proportion = (xlength)/(xlength + ylength)
		num_x_divs = sqrt((xlength_proportion*resolution)/(1 - xlength_proportion))
		num_y_divs = resolution/num_x_divs
		num_x_divs = Int(floor(num_x_divs))
		num_y_divs = Int(floor(num_y_divs))
	end
			
	x_values = range(xlims[1], xlims[2], num_x_divs)
	y_values = range(ylims[1], ylims[2], num_y_divs)
	shift_amount = (x_values[2] - x_values[1])/2
	parameters_1 = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]
	
	#=
	start_parameters = randn(ComplexF64, length(parameters(system(EP))))
	start_solutions = solutions(HomotopyContinuation.solve(system(EP), target_parameters= start_parameters))
	subdivision_data = HomotopyContinuation.solve(system(EP), start_solutions, start_parameters = start_parameters, target_parameters = parameters_1)
	=#
    function_values = []
	
	for s in subdivision_data
		push!(function_values, (s[2], fibre_function(s)))
	end

	polygons = []
	
	parameter_array = reshape(parameters_1, length(x_values), length(y_values))
	for i in 1:length(x_values)-1
		for j in 1:length(y_values)-1
			tl = parameter_array[i,j]
			tr = parameter_array[i+1,j]
			bl = parameter_array[i, j+1]
			br = parameter_array[i+1, j+1]
			tl_index = findfirst(x->x[1] == tl, function_values)
			tr_index = findfirst(x->x[1] == tr, function_values)
			bl_index = findfirst(x->x[1] == bl, function_values)
			br_index = findfirst(x->x[1] == br, function_values)

		if isodd(i)
			push!(polygons, [tl_index, tr_index, bl_index])
			push!(polygons, [bl_index, tr_index, br_index])
		else
			push!(polygons, [tl_index, tr_index, br_index])
			push!(polygons, [bl_index, tl_index, br_index])
		end
		end
	end
	
	initial_graphmesh = GraphMesh(function_values)
	initial_subdivision = Subdivision(initial_graphmesh,polygons)
	
	return initial_graphmesh, initial_subdivision
end


function visualize(GM::GraphMesh)
    #Create plot, display plot, return plot
end


function visualize(EP::EnumerativeProblem; fibre_function = n_real_points,kwargs)
    #Create graph mesh
    #Refine graph mesh to user specs
    #Simplify unnecessary mesh points from the point of view of the plot 
    #Call visualize on graph mesh and return that output, all while cache-ing the graphmesh
end

