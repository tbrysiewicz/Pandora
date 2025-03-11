import Base: getindex, iterate

using Plots: scatter, plot, plot!, Shape, cgrad
using DelaunayTriangulation: triangulate, each_solid_triangle, triangle_vertices, get_point
using LinearAlgebra: qr



export
    initialize_valued_subdivision,
    visualize,
    input_points,
    output_values,
    graph_mesh,
    polygons,
	function_cache

#GRAPHMESH
#Never re-order the function_cache
mutable struct GraphMesh
    function_cache :: Vector{Tuple{Vector{Float64},Float64}}
end

input_points(GM::GraphMesh) = getindex.(GM.function_cache, 1)
output_values(GM::GraphMesh) = getindex.(GM.function_cache, 2)

#Introduces functionality of a call like GM[4] to get the 4-th function cache pair in GM
Base.getindex(GM::GraphMesh,GM_index) = GM.function_cache[GM_index]

#Allows you to iterate over a GM (iterate over pairs in the function_cache)
function Base.iterate(GM::GraphMesh, state = 1)
	if state <= length(GM.function_cache)
		return GM.function_cache[state], state + 1
	else
		return nothing
	end
end

#VALUEDSUBDIVISION
mutable struct ValuedSubdivision
    GM :: GraphMesh
    Polygons :: Vector{Vector{Int64}} #Indices 
end

#GETTERS
graph_mesh(SD::ValuedSubdivision) = SD.GM
polygons(SD::ValuedSubdivision) = SD.Polygons
function_cache(GM::GraphMesh) = GM.function_cache

#SETTERS
function set_polygons!(SD::ValuedSubdivision, P::Vector{Vector{Int64}})
	SD.Polygons = P
end

function push_to_graph_mesh!(GM::GraphMesh, V::Tuple{Vector{Float64},Float64})
	push!(GM.function_cache, V)
end


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

#REFINEMENT
function is_complete(p::Vector{Int}, GM::GraphMesh) #checks if all polygon vertices share the same fibre function value. 
	v1_value = GM[p[1]][2]
	for i in p[2:end]
		vval = GM[i][2]
		if vval!=v1_value
			return(false)
		end
	end
	return(true)
end

function point_insertion(p::Vector{Int}, GM::GraphMesh) #returns a random point in the convex hull of a triangle from the graphmesh
	polygon_vertices = [GM[x][1] for x in p]
	c = randn(Float64, 3)
	c = map(abs, c)
	c = c./sum(c)
	random_convex_combination = [polygon_vertices[i].*c[i] for i in 1:3]
	random_convex_combination = sum(random_convex_combination)
	return [random_convex_combination]
end

function midpoint(p,q)
    return((p+q)./2)
end

function triforce_point_insertion(p::Vector{Int}, GM::GraphMesh) #returns "triforce" points
	polygon_vertices = [GM[x][1] for x in p]
	midpoint_1 = midpoint(polygon_vertices[1], polygon_vertices[2])
	midpoint_2 = midpoint(polygon_vertices[1], polygon_vertices[3])
	midpoint_3 = midpoint(polygon_vertices[2], polygon_vertices[3])
	return [midpoint_1, midpoint_2, midpoint_3]
end

function delaunay_triangulation!(SD::ValuedSubdivision)
    vertices = []
    for v in graph_mesh(SD)
        push!(vertices, v[1])
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	triangles::Vector{Vector{Int64}} = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, function_cache(graph_mesh(SD)))
        index_2 = findfirst(x->x[1] == vertex_2, function_cache(graph_mesh(SD)))
        index_3 = findfirst(x->x[1] == vertex_3, function_cache(graph_mesh(SD)))
		push!(triangles, [index_1, index_2, index_3])
	end
    set_polygons!(SD, triangles)
	return nothing
end

function refine!(SD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; fibre_function = x->n_real_solutions(x))
	new_parameters::Vector{Vector{Float64}} = []
	resolution_used = 0
	for T in polygons(SD)
		resolution_used >= resolution && break #checking to see if there is any resolution "left" to insert another point
		if is_complete(T, graph_mesh(SD)) == false
			parameters_to_insert_in_polygon = point_insertion(T, graph_mesh(SD))
			for i in parameters_to_insert_in_polygon 
				push!(new_parameters, i)
			end
			resolution_used +=length(parameters_to_insert_in_polygon)
		end
	end
	inserted_point_solutions = solve(EP, new_parameters)
	length(inserted_point_solutions) != length(new_parameters) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		push_to_graph_mesh!(graph_mesh(SD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end
	delaunay_triangulation!(SD)
	return resolution_used
end

#PLOTTING
function draw_triangle(triangle, color_value, GM::GraphMesh; label = false, label_text = "real solutions")
	triangle = [input_points(GM)[i] for i in triangle]
	triangle = Shape([(t[1], t[2]) for t in triangle])
	c = cgrad(:thermal, rev = false)[color_value]
	if label == true
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, label = label_text)
	else
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, label = false)
	end
end

function draw_valued_subdivision(SD::ValuedSubdivision; xlims = [-1,1], ylims = [-1,1]) #currently this function can only plot a ValuedSubdivision with a discrete fibre function
	my_plot = plot(xlims = xlims, ylims = ylims, legend = true)
	plotting_values = unique(output_values(graph_mesh(SD)))
	plotting_values = sort(plotting_values)
	values_that_have_been_plotted = []
	for i in plotting_values
		current_polygons = filter(x->graph_mesh(SD)[x[1]][2] == i, polygons(SD)) #filters all polygons whose first vertices have the fibre function value i
		color_value = findfirst(x->x == i, plotting_values)/length(plotting_values)
		for j in current_polygons
			if is_complete(j, graph_mesh(SD))
				if (i in values_that_have_been_plotted) == false
					draw_triangle(j, color_value, graph_mesh(SD), label = true, label_text = "$i real solutions")
					push!(values_that_have_been_plotted, i)
				else
					draw_triangle(j, color_value, graph_mesh(SD), label = false)
				end
			else
				continue
			end
		end
	end
	return my_plot
end

#RESTRICTING EPs
#Note that these functions will have to be modified at some point so that a new EP isn't generated
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

#VISUALIZATION
function visualize(GM::GraphMesh)
    scatter(first.(input_points(GM)), last.(input_points(GM)), 
    zcolor = output_values(GM), legend = false, colorbar = true)
end

#TODO: Write this cleanly
function visualize(EP::EnumerativeProblem; xlims = [-1,1], ylims = [-1,1], initial_resolution = 1000, total_resolution = 3*initial_resolution, fibre_function = x->n_real_solutions(x))
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end
	VSD = initialize_valued_subdivision(new_EP, xlims = xlims, ylims = ylims, resolution = initial_resolution, fibre_function = fibre_function)
	remaining_resolution = total_resolution
	while remaining_resolution > 0
		resolution_used = refine!(VSD, new_EP, remaining_resolution, fibre_function = fibre_function)
		remaining_resolution -= resolution_used
	end
	my_plot = draw_valued_subdivision(VSD, xlims = xlims, ylims = ylims)
	return (VSD, my_plot)
end


