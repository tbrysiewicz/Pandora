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
	complete_polygons::Vector{Vector{Int}} = []
	incomplete_polygons::Vector{Vector{Int}} = []
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

			if isodd(i) #staggers rows of triangles in the initial mesh
				establish_completeness!([tl_index, tr_index, bl_index], complete_polygons, incomplete_polygons, initial_graphmesh)
				establish_completeness!([bl_index, tr_index, br_index], complete_polygons, incomplete_polygons, initial_graphmesh)
			else
				establish_completeness!([tl_index, tr_index, br_index], complete_polygons, incomplete_polygons, initial_graphmesh)
				establish_completeness!([bl_index, tl_index, br_index], complete_polygons, incomplete_polygons, initial_graphmesh)
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

function establish_completeness!(polygon::Vector{Int}, complete_polygons::Vector{Vector{Int}}, incomplete_polygons::Vector{Vector{Int}}, GM::GraphMesh) #determines whether a polygon is complete and accordingly pushes it to complete_polygons or incomplete_polygons
	is_complete(polygon, GM) ? push!(complete_polygons, polygon) : push!(incomplete_polygons, polygon)
	return nothing
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

function retriangulate_valued_subdivision!(SD::ValuedSubdivision)
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
		establish_completeness!([index_1, index_2, index_3], complete_triangles, incomplete_triangles, graph_mesh(SD))
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
		number_of_points_inserted = 0
		for i in parameters_to_insert_in_polygon 
			if (i in new_parameters) == false && (i in input_points(graph_mesh(SD))) == false #check to make sure that the inserted parameters are not already inserted into another polygon or already in the graph mesh
				push!(new_parameters, i)
				number_of_points_inserted += 1
			end
		end
		resolution_used += number_of_points_inserted
	end
	length(new_parameters) == 0 && error("Did not find any incomplete triangles")
	inserted_point_solutions = solve(EP, new_parameters)
	length(inserted_point_solutions) != length(new_parameters) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		push_to_graph_mesh!(graph_mesh(SD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end
	retriangulate_valued_subdivision!(SD)
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
	my_plot = plot(xlims = xlims, ylims = ylims, aspect_ratio = :equal; kwargs...)
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
	boundary_edges::Vector{Tuple{Int64, Int64}} = [k for (k, v) in edge_counter if v == 1]
	return boundary_edges
end

function find_component(triangles::Vector{Vector{Int}}) 
	found_new_triangles = true
	triangles_in_component = [triangles[1]]
	vertices_in_component = copy(triangles[1])
	classified_indices = [1]
	while found_new_triangles
		found_new_triangles = false
		for i in 1:length(triangles)
			if (i in classified_indices) == false
				current_triangle = triangles[i]
				if length(intersect(current_triangle, vertices_in_component)) > 1
					push!(triangles_in_component, current_triangle)
					for v in current_triangle
						if (v in vertices_in_component) == false
							push!(vertices_in_component, v)
						end
					end
					push!(classified_indices, i)
					found_new_triangles = true
				end
			end
		end
	end
	return triangles_in_component
end

function ordered_boundary(edges::Vector{Tuple{Int, Int}}) #takes as input a list of boundary edges as given by boundary_edges function and returns an ordered list of vertices that make up the boundary
	adjacency_dictionary = Dict()
	for e in edges
		push!(get!(adjacency_dictionary, e[1], Vector()), e[2])
		push!(get!(adjacency_dictionary, e[2], Vector()), e[1])
	end
	start_vertex = first(keys(adjacency_dictionary))
	ordered_vertices = [start_vertex]
	current_vertex = start_vertex
	while true
		neighbors = adjacency_dictionary[current_vertex]
		current_vertex = nothing
		for neighbor in neighbors
			if (neighbor in ordered_vertices) == false
				push!(ordered_vertices, neighbor)
				current_vertex = neighbor
				break
			end
		end
		if isnothing(current_vertex)
			push!(ordered_vertices, ordered_vertices[1])
			break
		else
			continue
		end
	end
	return ordered_vertices
end

function boundary_indices_to_points(vertices::Vector{Int}, VSD::ValuedSubdivision) #takes as input a list of indices of boundary vertices and returns a vector containing the corresponding parameter points
	boundary_points::Vector{Vector{Float64}} = []
	for v in vertices
		point = graph_mesh(VSD)[v][1]
		push!(boundary_points, point)
	end
	return boundary_points
end

function triangle_indices_to_points(triangles::Vector{Vector{Int}}, VSD::ValuedSubdivision) #takes a list of triangles in index form and returns a vector containing the list of triangles as triples of points in the parameter space
	GM = graph_mesh(VSD)
	triangles_as_points::Vector{Vector{Vector{Float64}}} = []
	for t in triangles
		vertex_1 = GM[t[1]][1]
		vertex_2 = GM[t[2]][1]
		vertex_3 = GM[t[3]][1]
		push!(triangles_as_points, [vertex_1, vertex_2, vertex_3])
	end
	return triangles_as_points
end

function signed_area(boundary_vertices::Vector{Vector{Float64}})
	sum = 0
	n = length(boundary_vertices)
	if boundary_vertices[1] == boundary_vertices[n]
		j = n - 1 #checks to see if the first vertex matches the last vertex and changes the iteration of the loop if so. This matters because the ordered_boundary function returns an ordered list of boundary vertices that includes the starting vertex twice (first index and last index)
	else
		j = n
	end
	for i in 1:j
		x_1, y_1 = boundary_vertices[i]
		x_2, y_2 = boundary_vertices[mod1(i+1, n)]
		sum += (x_1*y_2 - x_2*y_1)
	end
	return 0.5*sum
end

function component_boundaries(triangles::Vector{Vector{Int}}, VSD::ValuedSubdivision) #takes as input the triangles of a connnected component (triangles consisting of indices) and returns the boundaries associated with the connected component. This will include the outer boundary (in counter-clockwise orientation) and any holes if there are any (in clockwise orientation)
	boundary_edges1 = boundary_edges(triangles)
	used_edges::Vector{Tuple{Int64, Int64}} = [] #the boundary may consist of several components (e.g. outer boundary, holes) so we need to track which edges are used
	boundary_components = []
	while isempty(setdiff(boundary_edges1, used_edges)) == false
		current_boundary = ordered_boundary(setdiff(boundary_edges1, used_edges))
		push!(boundary_components, current_boundary)
		for e in boundary_edges1
			if (e in used_edges) == false && (e[1] in current_boundary || e[2] in current_boundary)
				push!(used_edges, e)
			end
		end
	end
	areas = []
	signed_areas =[]
	boundary_components_point_form = [] #at this statge the boundary components consist of lists of indices. We need the parameters associated with these indices in order to calculate areas.
	for b in boundary_components
		current_boundary_points = boundary_indices_to_points(b, VSD)
		current_area = signed_area(current_boundary_points)
		push!(signed_areas, current_area)
		push!(areas, abs(current_area))
		push!(boundary_components_point_form, current_boundary_points)
	end
	outer_boundary_index = findfirst(x->x == max(areas...), areas)
	correctly_ordered_boundary_components = [] #while at this current stage each of the boundary components is ordered, we need to correct their orientations (ccw orientation for outer boundary, cw orientation for holes)
	for i in eachindex(boundary_components_point_form)
		if i == outer_boundary_index
			if signed_areas[i] < 0
				push!(correctly_ordered_boundary_components, [reverse(boundary_components_point_form[i])]) #we put the boundary vector into another vector because this is the format required for DelaunayTriangulation.jl
			else
				push!(correctly_ordered_boundary_components, [boundary_components_point_form[i]])
			end
		else
			if signed_areas[i] < 0
				push!(correctly_ordered_boundary_components, [boundary_components_point_form[i]])
			else
				push!(correctly_ordered_boundary_components, [reverse(boundary_components_point_form[i])])
			end
		end
	end
	return correctly_ordered_boundary_components
end