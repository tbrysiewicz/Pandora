import Base: getindex

using Plots: scatter

export
    initialize_valued_subdivision,
    visualize,
    input_points,
    output_values,
    graph_mesh,
    polygons


#Never re-order the function_cache
mutable struct GraphMesh
    function_cache :: Vector{Tuple{Vector{Float64},Float64}}
end

input_points(GM::GraphMesh) = getindex.(GM.function_cache, 1)
output_values(GM::GraphMesh) = getindex.(GM.function_cache, 2)

#Introduces functionality of a call like GM[4] to get the 4-th function cache pair in GM
Base.getindex(GM::GraphMesh,GM_index) = GM.function_cache[GM_index]


mutable struct ValuedSubdivision
    GM :: GraphMesh
    Polygons :: Vector{Vector{Int64}} #Indices 
end

graph_mesh(SD::ValuedSubdivision) = SD.GM
polygons(SD::ValuedSubdivision) = SD.Polygons





function initialize_valued_subdivision(EP::EnumerativeProblem; xlims::Vector = [-1,1], ylims::Vector = [-1,1], 
								fibre_function = x->n_real_solutions(x), resolution = 1000)
	xlength = xlims[2] - xlims[1]
	ylength = ylims[2] - ylims[1]

    #Decide how many x-vals vs y-vals to produce approximately 'resolution' many data pts
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

    #Shift x-values by half a step every-other-row to make hexagonal lattice
    shift_amount = (x_values[2] - x_values[1])/2 	
    parameters_1 = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]

    #Solve for all pts in the initial hexagonal lattice
    subdivision_solutions = solve(EP, parameters_1)
    function_cache = []

    length(subdivision_solutions) != length(parameters_1) && error("Did not solve for each parameter")


    for i in eachindex(subdivision_solutions)
        push!(function_cache, (parameters_1[i], fibre_function(subdivision_solutions[i])))
    end


    #Create the polygons in the hexagonal lattice
	polygons = []	
	parameter_array = reshape(parameters_1, length(x_values), length(y_values))
	for i in 1:length(x_values)-1
		for j in 1:length(y_values)-1
			tl = parameter_array[i,j]
			tr = parameter_array[i+1,j]
			bl = parameter_array[i, j+1]
			br = parameter_array[i+1, j+1]
			tl_index = findfirst(x->x[1] == tl, function_cache)
			tr_index = findfirst(x->x[1] == tr, function_cache)
			bl_index = findfirst(x->x[1] == bl, function_cache)
			br_index = findfirst(x->x[1] == br, function_cache)

			if isodd(i)
				push!(polygons, [tl_index, tr_index, bl_index])
				push!(polygons, [bl_index, tr_index, br_index])
			else
				push!(polygons, [tl_index, tr_index, br_index])
				push!(polygons, [bl_index, tl_index, br_index])
			end
		end
	end
	
	initial_graphmesh = GraphMesh(function_cache)
	initial_subdivision = ValuedSubdivision(initial_graphmesh,polygons)
	
	return initial_subdivision
end






#Visualization
function visualize(GM::GraphMesh)
    scatter(first.(input_points(GM)), last.(input_points(GM)), 
    zcolor = output_values(GM), legend = false, colorbar = true)
end

#TODO: Write this cleanly
function visualize(VS::ValuedSubdivision)
	GM = graph_mesh(VS)
	#scatter([p[1] for p in input_points(GM)],[p[2] for p in input_points(GM)],zcolor = output_values(GM), legend = false, colorbar=true)
end


