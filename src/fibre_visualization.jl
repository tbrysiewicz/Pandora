import Base: getindex, iterate

using Plots: scatter, plot, plot!, Shape, cgrad
using DelaunayTriangulation: triangulate, each_solid_triangle, triangle_vertices, get_point
using LinearAlgebra: qr
using Graphs: SimpleGraph, add_edge!, connected_components



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
	CompletePolygons :: Vector{Vector{Int64}} #Indices
	IncompletePolygons :: Vector{Vector{Int64}} 
end

#GETTERS
graph_mesh(SD::ValuedSubdivision) = SD.GM
complete_polygons(SD::ValuedSubdivision) = SD.CompletePolygons
incomplete_polygons(SD::ValuedSubdivision) = SD.IncompletePolygons
function_cache(GM::GraphMesh) = GM.function_cache

#SETTERS
function set_complete_polygons!(SD::ValuedSubdivision, P::Vector{Vector{Int64}})
	SD.CompletePolygons = P
end

function set_incomplete_polygons!(SD::ValuedSubdivision, P::Vector{Vector{Int64}})
	SD.IncompletePolygons = P
end

function push_to_graph_mesh!(GM::GraphMesh, V::Tuple{Vector{Float64},Float64})
	push!(GM.function_cache, V)
end


function initialize_valued_subdivision(EP::EnumerativeProblem; 
	xlims::Vector = [-1,1], ylims::Vector = [-1,1], 
	fibre_function = x->n_real_solutions(x), initial_resolution = 1000,kwargs...)

	xlength = xlims[2] - xlims[1]
	ylength = ylims[2] - ylims[1]

    #Decide how many x-vals vs y-vals to produce approximately 'resolution' many data pts
	if xlength == ylength
		num_x_divs = Int(floor(sqrt(initial_resolution)))
		num_y_divs = Int(floor(sqrt(initial_resolution)))
	else
		xlength_proportion = (xlength)/(xlength + ylength)
		num_x_divs = sqrt((xlength_proportion*initial_resolution)/(1 - xlength_proportion))
		num_y_divs = initial_resolution/num_x_divs
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
	initial_graphmesh = GraphMesh(function_cache)

    #Create the polygons in the hexagonal lattice
	complete_polygons = []
	incomplete_polygons = []
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
				#these statements add the triangles to complete_polygons or incomplete_polygons depending on whether is_complete returns true
				is_complete([tl_index, tr_index, bl_index], initial_graphmesh) && push!(complete_polygons, [tl_index, tr_index, bl_index])
				is_complete([tl_index, tr_index, bl_index], initial_graphmesh) || push!(incomplete_polygons, [tl_index, tr_index, bl_index])
				
				is_complete([bl_index, tr_index, br_index], initial_graphmesh) && push!(complete_polygons, [bl_index, tr_index, br_index])
				is_complete([bl_index, tr_index, br_index], initial_graphmesh) || push!(incomplete_polygons, [bl_index, tr_index, br_index])
			else
				is_complete([tl_index, tr_index, br_index], initial_graphmesh) && push!(complete_polygons, [tl_index, tr_index, br_index])
				is_complete([tl_index, tr_index, br_index], initial_graphmesh) || push!(incomplete_polygons, [tl_index, tr_index, br_index])

				is_complete([bl_index, tl_index, br_index], initial_graphmesh) && push!(complete_polygons,[bl_index, tl_index, br_index])
				is_complete([bl_index, tl_index, br_index], initial_graphmesh) || push!(incomplete_polygons, [bl_index, tl_index, br_index])
			end
		end
	end
	
	initial_subdivision = ValuedSubdivision(initial_graphmesh,complete_polygons, incomplete_polygons)
	
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
#Have not made any changes to delaunay_triangulation yet
function delaunay_triangulation!(SD::ValuedSubdivision)
    vertices = []
    for v in graph_mesh(SD)
        v[1] in vertices || push!(vertices, v[1]) #ensuring that the triangulation will not contain duplicate points and the package won't give that annoying warning
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	complete_triangles::Vector{Vector{Int64}} = []
	incomplete_triangles::Vector{Vector{Int64}} = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, function_cache(graph_mesh(SD)))
        index_2 = findfirst(x->x[1] == vertex_2, function_cache(graph_mesh(SD)))
        index_3 = findfirst(x->x[1] == vertex_3, function_cache(graph_mesh(SD)))
		if is_complete([index_1, index_2, index_3], graph_mesh(SD))
			push!(complete_triangles, [index_1, index_2, index_3])
		else
			push!(incomplete_triangles, [index_1, index_2, index_3])
		end
	end
    set_complete_polygons!(SD, complete_triangles)
	set_incomplete_polygons!(SD, incomplete_triangles)
	return nothing
end

#have not made any changes to the refine function yet
function refine!(SD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; 
	fibre_function = x->n_real_solutions(x), insertion_method = triforce_point_insertion,kwargs...)
	new_parameters::Vector{Vector{Float64}} = []
	resolution_used = 0
	for T in incomplete_polygons(SD)
		resolution_used >= resolution && break #checking to see if there is any resolution "left" to insert another point
		parameters_to_insert_in_polygon = insertion_method(T, graph_mesh(SD))
		for i in parameters_to_insert_in_polygon 
			push!(new_parameters, i)
		end
		resolution_used += length(parameters_to_insert_in_polygon)
	end
	length(new_parameters) == 0 && error("Did not find any incomplete triangles")
	inserted_point_solutions = solve(EP, new_parameters)
	length(inserted_point_solutions) != length(new_parameters) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		(new_parameters[i] in input_points(graph_mesh(SD))) || push_to_graph_mesh!(graph_mesh(SD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i])))) #checking to make sure we don't add duplicate points to graph_mesh
	end
	delaunay_triangulation!(SD)
	return resolution_used
end

#PLOTTING
function draw_triangle(triangle, color_value, GM::GraphMesh; label = false, 
	label_text = "real solutions", kwargs...)
	triangle = [input_points(GM)[i] for i in triangle]
	triangle = Shape([(t[1], t[2]) for t in triangle])
	c = cgrad(:thermal, rev = false)[color_value]
	if label == true
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, label = label_text;kwargs...)
	else
		plot!(triangle, fillcolor = c, linecolor = c, linewidth = false, label = false; kwargs...)
	end
end

function draw_valued_subdivision(SD::ValuedSubdivision; xlims = [-1,1],	ylims = [-1,1], kwargs...) #currently this function can only plot a ValuedSubdivision with a discrete fibre function
	my_plot = plot(xlims = xlims, ylims = ylims; kwargs...)
	plotting_values = unique(output_values(graph_mesh(SD)))
	plotting_values = sort(plotting_values)
	values_that_have_been_plotted = []
	for i in plotting_values
		current_polygons = filter(x->graph_mesh(SD)[x[1]][2] == i, complete_polygons(SD)) #filters all polygons whose first vertices have the fibre function value i
		color_value = findfirst(x->x == i, plotting_values)/length(plotting_values)
		for j in current_polygons
				if (i in values_that_have_been_plotted) == false
					draw_triangle(j, color_value, graph_mesh(SD), label = true, label_text = "$i real solutions"; kwargs...)
					push!(values_that_have_been_plotted, i)
				else
					draw_triangle(j, color_value, graph_mesh(SD), label = false; kwargs...)
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

#=
xlims = [-1,1], ylims = [-1,1], initial_resolution = 1000, 
total_resolution = 3*initial_resolution, 
fibre_function = x->n_real_solutions(x), kwargs...
=#
#TODO: Write this cleanly
function visualize(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 3*initial_resolution,  
	fibre_function = x-> n_real_solutions(x),kwargs...)
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end
	VSD = initialize_valued_subdivision(new_EP; xlims = xlims, ylims = ylims, 
	fibre_function = fibre_function, kwargs...)
	remaining_resolution = total_resolution
	while remaining_resolution > 0
		resolution_used = refine!(VSD, new_EP, remaining_resolution; fibre_function = fibre_function, kwargs...)
		remaining_resolution -= resolution_used
	end
	my_plot = draw_valued_subdivision(VSD; xlims = xlims, ylims = ylims, kwargs...)
	return (VSD, my_plot)
end

function polygon_adjacency(polygon_1::Vector{Int64}, polygon_2::Vector{Int64}) #returns true if the two inputted polygons share any vertices and are therefore adjacent. Returns false otherwise.
	polygon_1 = Set(polygon_1)
	polygon_2 = Set(polygon_2)
	if isempty(intersect(polygon_1, polygon_2)) == false
		return true
	else
		return false
	end
end

function triangulation_components(triangles::Vector{Vector{Int64}}) #given a list of triangles that may form a triangulation or the union of several triangulations, will return the connected components. E.g. if you gave it a list of all of the complete triangles that have 3 real solutions it will return a list of triangles partitioned by which chamber they belong to.
	n = length(triangles)
	G = SimpleGraph(n)
	for i in 1:n-1
		for j in i+1:n
			polygon_adjacency(triangles[i], triangles[j]) && add_edge!(G, i, j)
		end
	end
	components = connected_components(G)
	components = [[triangles[c] for c in C] for C in components]
	return components
end

function boundary_edges(triangles::Vector{Vector{Int64}}) #given a list of triangles from a triangulation, will return a list of boundary edges of the triangulation (as tuples of vertices)
	edge_counter = Dict()
	for t in triangles
		edge_1 = (min(t[1], t[2]), max(t[1], t[2]))
		edge_2 = (min(t[3], t[2]), max(t[3], t[2]))
		edge_3 = (min(t[1], t[3]), max(t[1], t[3]))
		get!(edge_counter, edge_1, 0)
		get!(edge_counter, edge_2, 0)
		get!(edge_counter, edge_3, 0)
		edge_counter[edge_1] += 1
		edge_counter[edge_2] += 1
		edge_counter[edge_3] += 1
	end
	boundary_edges = [k for (k, v) in edge_counter if v == 1]
	return boundary_edges
end