#I pulled the following three functions from main branch, only changed getters
function gram_schmidt(basis_vectors::Vector)
	M = hcat(basis_vectors...)
	Q = Matrix(qr(M).Q)
	return(collect(eachcol(Q)))
end

#Given P = [p_1...p_k] parameters, this function considers the affine span of P 
#   as p_1+span(p_i-p_1) where span(p_i-p_1) = span(b_1...b_k-1) where the bi's are
#   an orthonormal basis for the span. 

#TODO: Which cache values can be inherited by the restriction?

function restrict(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
	n = length(P)
    @var t[1:n-1]
    basis_vectors = gram_schmidt([(P[i]-P[1]) for i in 2:n])
    affine_span = P[1] + sum([t[i].*basis_vectors[i] for i in 1:n-1])
    new_expressions = [subs(f,parameters(EP)=>affine_span) for f in expressions(EP)]
    return(EnumerativeProblem(System(new_expressions,variables=variables(EP),parameters=t)))
end

#If the user doesn't care WHICH plane to restrict to, do one through the current base
# parameters so that EP inherits the solutions

function planar_restriction(EP::EnumerativeProblem)
	P = [randn(Float64,n_parameters(EP)) for i in 1:3]
	return(restrict(EP,P))
end

function initialize_subdivision(EP::EnumerativeProblem, xlims::Vector, ylims::Vector, 
								fibre_function = x->n_real_solutions(x), resolution = 1000)
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
	subdivision_solutions = solve(EP, parameters_1)
    function_cache = []
	
	for i in subdivision_solutions
		push!(function_cache, (parameters_1[findfirst(x->x==i, subdivision_solutions)], fibre_function(i)))
	end

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
	initial_subdivision = Subdivision(initial_graphmesh,polygons)
	
	return initial_subdivision
end
#=
function draw_triangle(triangle, color_value; label = false, label_text = "real solutions")
	triangle = Shape([(t[1], t[2]) for t in triangle])
	c = cgrad(:thermal, rev = false)[color_value]
	if label == true
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, labels = label_text)
	else
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, labels = false)
	end
end


function draw_subdivision(sd::Subdivision, xlims::Vector, ylims::Vector)
	plotting_values = []
	for i in subdivision_1.GM.function_cache
		if (i[2] in plotting_values) == false
			push!(plotting_values, i[2])
		end
	end
	plotting_values = sort(plotting_values)
	number_of_plotting_values = length(plotting_values)
	my_plot = plot(xlims = xlims, ylims = ylims, legend = true)
	values_included_in_legend = []
	for T in sd.Polygons
		if is_complete(T, sd.GM)[1] == true
			current_plotting_value = value_from_index(T[1], sd.GM)[2]
			if (current_plotting_value in values_included_in_legend) == false
				draw_triangle(T, findfirst(x->x == current_plotting_value, plotting_values)/number_of_plotting_values, label = true, label_text = "$(current_plotting_value) real solutions")
				push!(values_included_in_legend, current_plotting_value)
			else
				draw_triangle(T, findfirst(x->x == current_plotting_value, plotting_values)/number_of_plotting_values)
			end
		else
			continue
		end
	end
	return my_plot
end

function visualize(EP::EnumerativeProblem; xlims = [-2,2], ylims = [-2,2], fibre_function = x->n_real_solutions(x), initial_resolution = 1000, refinement_resolution = 4*1000)
	if n_parameters(EP) > 2
		EP = planar_restriction(EP)
	end
	sd = initialize_subdivision(EP, xlims, ylims, fibre_function, initial_resolution)
	while refinement_resolution > 0
		resolution_used = refine!(sd, EP, refinement_resolution, fibre_function = fibre_function)
		refinement_resolution -= resolution_used
		if refinement_resolution <= 3
			break
		end
	end
	my_plot = draw_subdivision(sd, xlims, ylims)
	return my_plot
end
=#
function visualize(GM::GraphMesh)
	scatter([p[1] for p in input_points(GM)],[p[2] for p in input_points(GM)],zcolor = output_values(GM), legend = false, colorbar=true)
end

#TODO: Write this cleanly
function visualize(SD::Subdivision)
	GM = graph_mesh(SD)
	#scatter([p[1] for p in input_points(GM)],[p[2] for p in input_points(GM)],zcolor = output_values(GM), legend = false, colorbar=true)

end




function visualize(EP::EnumerativeProblem; fibre_function = n_real_points,kwargs)
    #Create graph mesh
    #Refine graph mesh to user specs
    #Simplify unnecessary mesh points from the point of view of the plot 
    #Call visualize on graph mesh and return that output, all while cache-ing the graphmesh
end

