import Base: getindex, iterate

using Plots: scatter, plot, plot!, Shape, cgrad
using DelaunayTriangulation: triangulate, each_solid_triangle, triangle_vertices, get_point, convert_boundary_points_to_indices
using LinearAlgebra: qr



export
	initialize_trihexagonal_mesh,
    initialize_valued_subdivision,
    visualize,
    input_points,
    output_values,
    graph_mesh,
    polygons,
	function_cache,
	complete_polygons,
	incomplete_polygons

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

#delete a specific polygon from IncompletePolygons
function delete_from_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64})
	polygon_index = findfirst(x -> x == P, incomplete_polygons(VSD))
	if polygon_index !== nothing
		deleteat!(VSD.IncompletePolygons, polygon_index)
	else
		error("Polygon not found in IncompletePolygons!")
	end
end

function push_to_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}) #push a single polygon to complete polygons
	push!(VSD.CompletePolygons, P)
end

function push_to_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64})
	push!(VSD.IncompletePolygons, P)
end

function push_to_graph_mesh!(GM::GraphMesh, V::Tuple{Vector{Float64},Float64})
	push!(GM.function_cache, V)
end

#region MESH FUNCTIONS
function trihexagonal_mesh(;xlims::Vector = [-1,1], ylims::Vector = [-1,1], resolution = 1000, kwargs...)

	x_values, y_values = initial_parameter_distribution(xlims = xlims, ylims = ylims, resolution = resolution)

    #Shift x-values by half a step every-other-row to make hexagonal lattice
    shift_amount = (x_values[2] - x_values[1])/2 	
    parameters = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]

	#Making the parameters as a matrix makes it a bit easier to index the triangles
	parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))

	triangles = Vector{Vector{Int}}([])
	for i in 1:length(x_values)-1
		for j in 1:length(y_values)-1
			tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
			tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
			bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
			br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)

			if isodd(i) #staggers rows of triangles in the initial mesh
				push!(triangles,[tl_index, tr_index, bl_index])
				push!(triangles,[bl_index, tr_index, br_index])
			else
				push!(triangles,[tl_index, tr_index, br_index])
				push!(triangles,[bl_index, tl_index, br_index])
			end
		end
	end
	return(triangles,parameters)
end

function rectangular_mesh(;xlims::Vector = [-1,1], ylims::Vector = [-1,1], resolution = 1000, kwargs...)
	x_values, y_values = initial_parameter_distribution(xlims = xlims, ylims = ylims, resolution = resolution)
	parameters = [[i,j] for i in x_values for j in y_values]
	parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))
	rectangles = Vector{Vector{Int}}([])
	for i in 1:length(x_values)-1
		for j in 1:length(y_values)-1
			tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
			tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
			bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
			br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)
			push!(rectangles, [tl_index, tr_index, br_index, bl_index]) #choosing a standardized ordering for rectangle vertices matters
		end
	end
	return (rectangles, parameters)
end

#determines initial parameters based on visualization window and resolution
function initial_parameter_distribution(;xlims = [-1,1], ylims = [-1,1], resolution = 1000)
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

	return x_values, y_values
end

#TODO: Either use kwargs cleverly with the mesh function or eliminate it in the signature
function initialize_valued_subdivision(EP::EnumerativeProblem; 
	xlims::Vector = [-1,1], ylims::Vector = [-1,1], 
	fibre_function = x->n_real_solutions(x), initial_resolution = 1000,
	mesh_function = trihexagonal_mesh, kwargs...)


	(polygons,parameters) = mesh_function(; xlims=xlims, ylims=ylims, resolution=initial_resolution)
	
    #Solve for all pts in the initial hexagonal lattice
    subdivision_solutions = solve(EP, parameters)

    length(subdivision_solutions) != length(parameters) && error("Did not solve for each parameter")

	#Create the GraphMesh via input/output values
    function_cache = []
    for i in eachindex(subdivision_solutions)
        push!(function_cache, (parameters[i], fibre_function(subdivision_solutions[i])))
    end
	initial_graphmesh = GraphMesh(function_cache)

	#Now we scroll through polygons and establish completeness
	complete_polygons::Vector{Vector{Int}} = []
	incomplete_polygons::Vector{Vector{Int}} = []

	for t in polygons
		is_complete(t, initial_graphmesh) ? push!(complete_polygons, t) : push!(incomplete_polygons, t)
	end

	initial_subdivision = ValuedSubdivision(initial_graphmesh, complete_polygons, incomplete_polygons)
	return initial_subdivision
end

#RESTRICTING EPS
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

#REFINEMENT FUNCTIONS

function global_delaunay_refine!(SD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; 
	fibre_function = x->n_real_solutions(x), insertion_method = triforce_point_insertion,kwargs...)
	new_parameters, resolution_used = refinement_parameters(SD, insertion_method, resolution)
	length(new_parameters) == 0 && error("Did not find any incomplete triangles")
	inserted_point_solutions = solve(EP, new_parameters)
	length(inserted_point_solutions) != length(new_parameters) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		push_to_graph_mesh!(graph_mesh(SD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end
	retriangulate_valued_subdivision!(SD)
	return resolution_used
end

function rectangular_refine!(SD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; 
	fibre_function = x->n_real_solutions(x), insertion_method = rectangular_refine!,kwargs...)
	new_parameters, resolution_used, selected_incomplete_rectangles, newly_created_rectangles = refinement_parameters_and_new_rectangles(SD, resolution)
	length(new_parameters) == 0 && error("Did not find any incomplete triangles")
	inserted_point_solutions = solve(EP, new_parameters)
	length(inserted_point_solutions) != length(new_parameters) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		push_to_graph_mesh!(graph_mesh(SD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end
	for R in selected_incomplete_rectangles
		delete_from_incomplete_polygons!(SD, R)
	end
	for R in newly_created_rectangles
		rectangle = [findfirst(x->x[1] == y, function_cache(graph_mesh(SD))) for y in R]
		is_complete(rectangle, graph_mesh(SD)) ? push_to_complete_polygons!(SD, rectangle) : push_to_incomplete_polygons!(SD, rectangle)
	end
	return resolution_used
end

#determines the parameters to be inserted into the graph mesh in a single stage of refinement
#TODO: Need to fix this so that polygons are randomly selected from the incomplete_polygons list. 
function refinement_parameters(VSD::ValuedSubdivision, insertion_method = triforce_point_insertion, resolution = 1000)
	new_parameters::Vector{Vector{Float64}} = []
	resolution_used = 0
	for T in incomplete_polygons(VSD)
		resolution_used >= resolution && break #checking to see if there is any resolution "left" to insert another point
		parameters_to_insert_in_polygon = insertion_method(T, graph_mesh(VSD))
		number_of_points_inserted = 0
		for i in parameters_to_insert_in_polygon 
			if !(i in new_parameters) && !(i in input_points(graph_mesh(VSD))) #check to make sure that the inserted parameters are not already inserted into another polygon or already in the graph mesh
				push!(new_parameters, i)
				number_of_points_inserted += 1
			end
		end
		resolution_used += number_of_points_inserted
	end
	return new_parameters, resolution_used
end

#determines the parameters to be inserted into the graph mesh for rectangular refinement and also tracks which rectangles need to be deleted and added
#TODO: Need to fix this so that polygons are randomly selected from the incomplete_polygons list. 
function refinement_parameters_and_new_rectangles(VSD::ValuedSubdivision, resolution = 1000)
	new_parameters::Vector{Vector{Float64}} = []
	selected_incomplete_rectangles::Vector{Vector{Int64}} = []
	newly_created_rectangles::Vector{Vector{Vector{Float64}}} = []
	resolution_used = 0
	for T in incomplete_polygons(VSD)
		resolution_used >= resolution && break #checking to see if there is any resolution "left" to insert another point
		parameters_to_insert_in_polygon = rectangular_point_insertion(T, graph_mesh(VSD))
		number_of_points_inserted = 0
		for i in parameters_to_insert_in_polygon 
			if !(i in new_parameters) && !(i in input_points(graph_mesh(VSD))) #check to make sure that the inserted parameters are not already inserted into another polygon or already in the graph mesh
				push!(new_parameters, i)
				number_of_points_inserted += 1
			end
		end
		resolution_used += number_of_points_inserted
		number_of_points_inserted > 0 && push!(selected_incomplete_rectangles, T)
		push!(newly_created_rectangles, [graph_mesh(VSD)[T[1]][1], parameters_to_insert_in_polygon[1], parameters_to_insert_in_polygon[5], parameters_to_insert_in_polygon[4]])
		push!(newly_created_rectangles, [parameters_to_insert_in_polygon[1], graph_mesh(VSD)[T[2]][1], parameters_to_insert_in_polygon[2], parameters_to_insert_in_polygon[5]])
		push!(newly_created_rectangles, [parameters_to_insert_in_polygon[5], parameters_to_insert_in_polygon[2], graph_mesh(VSD)[T[3]][1], parameters_to_insert_in_polygon[3]])
		push!(newly_created_rectangles, [parameters_to_insert_in_polygon[4], parameters_to_insert_in_polygon[5], parameters_to_insert_in_polygon[3], graph_mesh(VSD)[T[4]][1]])
	end
	return new_parameters, resolution_used, selected_incomplete_rectangles, newly_created_rectangles
end

function delaunay_valued_subdivision(GM::GraphMesh)
    vertices = []
    for v in GM
        v[1] in vertices || push!(vertices, v[1]) #ensuring that the triangulation will not contain duplicate points and the package won't give that annoying warning
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices) #This comes from DelaunayTriangulations
	triangle_iterator = each_solid_triangle(tri) #Also from DT
	complete_triangles::Vector{Vector{Int64}} = []
	incomplete_triangles::Vector{Vector{Int64}} = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, function_cache(GM))
        index_2 = findfirst(x->x[1] == vertex_2, function_cache(GM))
        index_3 = findfirst(x->x[1] == vertex_3, function_cache(GM))
		polygon = [index_1,index_2,index_3]
		is_complete(polygon, GM) ? push!(complete_triangles, polygon) : push!(incomplete_triangles, polygon)
	end
	VSD = ValuedSubdivision(GM,complete_triangles,incomplete_triangles)
	return(VSD)
end

#TODO: Why did VSD = delaunay_valued_subdivision(graph_mesh(VSD))
#  not overwrite VSD?
function retriangulate_valued_subdivision!(VSD::ValuedSubdivision)
	new_VSD = delaunay_valued_subdivision(graph_mesh(VSD))
	set_complete_polygons!(VSD, complete_polygons(new_VSD))
	set_incomplete_polygons!(VSD, incomplete_polygons(new_VSD))
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

 #returns "triforce" points
function triforce_point_insertion(p::Vector{Int}, GM::GraphMesh)
	polygon_vertices = [GM[x][1] for x in p]
	midpoint_1 = midpoint(polygon_vertices[1], polygon_vertices[2])
	midpoint_2 = midpoint(polygon_vertices[1], polygon_vertices[3])
	midpoint_3 = midpoint(polygon_vertices[2], polygon_vertices[3])
	return [midpoint_1, midpoint_2, midpoint_3]
end

#inserts 4 points into rectangle (center of rectangle and midpoint of each edge)
function rectangular_point_insertion(p::Vector{Int}, GM::GraphMesh)
	rectangle_vertices = [GM[x][1] for x in p]
	center = [sum([x[1] for x in rectangle_vertices])/4, sum([x[2] for x in rectangle_vertices])/4]
	midpoint1 = midpoint(GM[p[1]][1], GM[p[2]][1])
	midpoint2 = midpoint(GM[p[2]][1], GM[p[3]][1])
	midpoint3 = midpoint(GM[p[3]][1], GM[p[4]][1])
	midpoint4 = midpoint(GM[p[4]][1], GM[p[1]][1])
	return [midpoint1, midpoint2, midpoint3, midpoint4, center] #is this the right formatting? Should there be another set of brackets around this?
end

#checks if all polygon vertices share the same fibre function value. 
#TODO: Make it work for non-discrete valued functions
function is_complete(p::Vector{Int}, GM::GraphMesh) 
	v1_value = GM[p[1]][2]
	for i in p[2:end]
		vval = GM[i][2]
		if vval!=v1_value
			return(false)
		end
	end
	return(true)
end

#VISUALIZATION
function visualize(GM::GraphMesh)
    scatter(first.(input_points(GM)), last.(input_points(GM)), 
    zcolor = output_values(GM), legend = false, colorbar = true)
end

#TODO: Really use kwargs... here for passing xlims/ylims to initialize valued subdivision
function visualize(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 4*initial_resolution,  
	fibre_function = x-> n_real_solutions(x), refine_function! = global_delaunay_refine!, mesh_function = trihexagonal_mesh,
	kwargs...)

	#Check enumerative problem is visualizable
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end

	VSD = initialize_valued_subdivision(new_EP; xlims = xlims, ylims = ylims, 
	fibre_function = fibre_function, mesh_function = mesh_function, kwargs...)

	remaining_resolution = total_resolution - initial_resolution
	while remaining_resolution > 0
		resolution_used = refine_function!(VSD, new_EP, remaining_resolution; fibre_function = fibre_function, kwargs...)
		remaining_resolution -= resolution_used
	end

	my_plot = draw_valued_subdivision(VSD; xlims = xlims, ylims = ylims, kwargs...)
	return (VSD, my_plot)
end


#Everything after this point is written for local Delaunay refinement (better name?)

function polygon_adjacency(polygon_1::Vector{Int64}, polygon_2::Vector{Int64}) #returns true if the two inputted polygons share any vertices and are therefore adjacent. Returns false otherwise.
	polygon_1 = Set(polygon_1)
	polygon_2 = Set(polygon_2)
	if isempty(intersect(polygon_1, polygon_2)) == false
		return true
	else
		return false
	end
end

#given a list of triangles from a triangulation, 
#will return a set of boundary edges of the triangulation (as tuples of vertices)
function boundary_edges(triangles::Vector{Vector{Int64}}) 
	edge_counter = Dict()
	for t in triangles
		T = sort(t)
		edge_1 = T[[1,2]]
		edge_2 = T[[1,3]]
		edge_3 = T[[2,3]]
		get!(edge_counter, edge_1, 0)
		get!(edge_counter, edge_2, 0)
		get!(edge_counter, edge_3, 0)
		edge_counter[edge_1] += 1
		edge_counter[edge_2] += 1
		edge_counter[edge_3] += 1
	end
	boundary_edges = Set([k for (k, v) in edge_counter if v == 1])
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

function find_components(triangles::Vector{Vector{Int}})
	components = []
	classified_triangles = Set{Int}() #indices
	while length(classified_triangles) < length(triangles)
		starting_triangle_index = findfirst(x->!(x in classified_triangles), 1:length(triangles))
		current_component = [starting_triangle_index] #indices
		vertices_in_component = Set(triangles[starting_triangle_index])
		push!(classified_triangles, starting_triangle_index)
		found_new_triangles = true
		while found_new_triangles
			found_new_triangles = false
			for i in 1:length(triangles)
				if !(i in classified_triangles)
					shared = 0
					for v in triangles[i]
						shared += v in vertices_in_component ? 1 : 0
						shared > 1 && break
					end
					if shared > 1
						push!(current_component, i)
						for v in triangles[i]
							!(v in vertices_in_component) && push!(vertices_in_component, v)
						end
						push!(classified_triangles, i)
						found_new_triangles = true
					end
				else
					continue
				end
			end
		end
		push!(components, [triangles[i] for i in current_component])
	end
	return components
end

#takes as input a list of boundary edges as given by boundary_edges 
#function and returns an ordered list of vertices that make up the boundary

function ordered_boundary(edges::Set{Vector{Int}}) 
	adjacency_dictionary = Dict()
	for e in edges
		push!(get!(adjacency_dictionary, e[1], Vector()), e[2])
		push!(get!(adjacency_dictionary, e[2], Vector()), e[1])
	end
	start_vertex = first(keys(adjacency_dictionary))
	ordered_vertices = [start_vertex]
	classified_vertices = Set(ordered_vertices) #track which vertices have been added to the ordered boundary
	current_vertex = start_vertex
	while true
		neighbors = adjacency_dictionary[current_vertex]
		current_vertex = nothing
		for neighbor in neighbors
			if !(neighbor in classified_vertices)
				push!(ordered_vertices, neighbor)
				push!(classified_vertices, neighbor)
				current_vertex = neighbor
				break
			end
		end
		if isnothing(current_vertex)
			push!(ordered_vertices, ordered_vertices[1]) #if no more unclassified, adjacent vertices are found, ends boundary by adding the start vertex to the end
			break
		else
			continue
		end
	end
	return ordered_vertices
end

function boundary_indices_to_points(vertices::Vector{Int}, VSD::ValuedSubdivision) #takes as input a list of indices of boundary vertices and returns a vector containing the corresponding parameter points
	return [graph_mesh(VSD)[v][1] for v in vertices]
end

function triangle_indices_to_points(triangles::Vector{Vector{Int}}, VSD::ValuedSubdivision) #takes a list of triangles in index form and returns a vector containing the list of triangles as triples of points in the parameter space
	return [[graph_mesh(VSD)[t][1] for t in T] for T in triangles]
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

function component_boundaries(triangles::Vector{Vector{Int}}, VSD::ValuedSubdivision) #takes as input the triangles of a connected component (triangles consisting of indices) and returns the boundaries associated with the connected component. This will include the outer boundary (in counter-clockwise orientation) and any holes if there are any (in clockwise orientation)
	boundary_edges1 = boundary_edges(triangles)
	used_edges = Set() #the boundary may consist of several components (e.g. outer boundary, holes) so we need to track which edges are used
	boundary_components = []
	while length(used_edges) < length(boundary_edges1)
		current_boundary = ordered_boundary(setdiff(boundary_edges1, used_edges))
		push!(boundary_components, current_boundary)
		for e in boundary_edges1
			if !(e in used_edges) && (e[1] in current_boundary || e[2] in current_boundary)
				push!(used_edges, e)
			end
		end
	end
	areas = []
	signed_areas = []
	boundary_components_point_form = [] #at this stage the boundary components consist of lists of indices. We need the parameters associated with these indices in order to calculate areas.
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

function refine_valued_subdivision!(VSD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int; fibre_function = x->n_real_solutions(x))
	new_parameters::Vector{Vector{Float64}} = []
	resolution_used = 0
	for T in incomplete_polygons(VSD)
		resolution_used >= resolution && break
		parameters_to_insert_in_polygon = triforce_point_insertion(T, graph_mesh(VSD))
		number_of_points_inserted = 0
		for i in parameters_to_insert_in_polygon 
			if !(i in new_parameters) && !(i in input_points(graph_mesh(VSD))) #check to make sure that the inserted parameters are not already inserted into another polygon or already in the graph mesh
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
		push_to_graph_mesh!(graph_mesh(VSD), (new_parameters[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end

	incomplete_vertices = Set(boundary_indices_to_points(vcat(incomplete_polygons(VSD)...), VSD)) #all vertices in incomplete_polygons
	boundary_vertices = Set()
	for p in new_parameters
		push!(incomplete_vertices, p) #the new parameters that have been inserted will not be in the incomplete_vertices list yet so we must manually add them
	end
	incomplete_components = find_components(incomplete_polygons(VSD))
	curves = []
	for c in incomplete_components 
		current_component_boundaries = component_boundaries(c, VSD)
		for x in current_component_boundaries
			push!(curves, x)
			for y in x
				for v in y
					!(v in boundary_vertices) && push!(boundary_vertices, v) #we need to track which vertices are part of the boundary
				end
			end
		end
	end
	interior_points_for_triangulation::Vector{Vector{Float64}} = []
	for p in incomplete_vertices
		!(p in boundary_vertices) && push!(interior_points_for_triangulation, p)
	end
	boundary_nodes, points = convert_boundary_points_to_indices(curves; existing_points = interior_points_for_triangulation)
	tri = triangulate(points; boundary_nodes)
	
	triangle_iterator = each_solid_triangle(tri)
	incomplete_triangles::Vector{Vector{Int64}} = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
        index_1 = findfirst(x->x[1] == vertex_1, function_cache(graph_mesh(VSD)))
        index_2 = findfirst(x->x[1] == vertex_2, function_cache(graph_mesh(VSD)))
        index_3 = findfirst(x->x[1] == vertex_3, function_cache(graph_mesh(VSD)))
		is_complete([index_1, index_2, index_3], graph_mesh(VSD)) ? push_to_complete_polygons!(VSD, [index_1, index_2, index_3]) : push!(incomplete_triangles, [index_1, index_2, index_3])
	end
	set_incomplete_polygons!(VSD, incomplete_triangles) #incomplete_polygons needs to be completely replaced after triangulation, complete_polygons does not

	return resolution_used
end

function draw_visualization(VSD::ValuedSubdivision; xlims = [-1,1],	ylims = [-1,1])
	my_plot = plot(xlims = xlims, ylims = ylims, aspect_ratio = :equal)
	plotting_values = unique(output_values(graph_mesh(VSD)))
	plotting_values = sort(plotting_values)
	values_that_have_been_plotted = []
	complete_components = find_components(complete_polygons(VSD))
	for c in complete_components
		boundaries = component_boundaries(c, VSD)
		if length(boundaries) == 1
			boundary_vertices = boundaries[1][1]
			shape_1 = Shape([(x[1], x[2]) for x in boundary_vertices])
			vertex_one_index = findfirst(x->x==boundary_vertices[1], input_points(graph_mesh(VSD)))
			function_value_for_component = output_values(graph_mesh(VSD))[vertex_one_index]
			color_value = findfirst(x->x == function_value_for_component, plotting_values)/length(plotting_values)
			c = cgrad(:thermal, rev = false)[color_value]
			if (function_value_for_component in values_that_have_been_plotted)
				plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = false)
			else
				plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = "$function_value_for_component real solutions")
				push!(values_that_have_been_plotted, function_value_for_component)
			end
		else
			vertex_one_index = findfirst(x->x==c[1][1], input_points(graph_mesh(VSD)))
			function_value_for_component = output_values(graph_mesh(VSD))[vertex_one_index]
			color_value = findfirst(x->x == function_value_for_component, plotting_values)/length(plotting_values)
			c = cgrad(:thermal, rev = false)[color_value]
			for t in c
				shape_1 = Shape([(x[1], x[2]) for x in t])
				if (function_value_for_component in values_that_have_been_plotted)
				plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = false)
				else
					plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = "$function_value_for_component real solutions")
					push!(values_that_have_been_plotted, function_value_for_component)

				end
			end
		end
	end
	return my_plot
end

function visualize_2(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 3*initial_resolution,  
	fibre_function = x-> n_real_solutions(x))
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end
	VSD = initialize_valued_subdivision(new_EP; xlims = xlims, ylims = ylims, 
	fibre_function = fibre_function)
	remaining_resolution = total_resolution
	while remaining_resolution > 0
		resolution_used = refine_valued_subdivision!(VSD, new_EP, remaining_resolution; fibre_function = fibre_function)
		remaining_resolution -= resolution_used
	end
	my_plot = draw_visualization(VSD; xlims = xlims, ylims = ylims)
	return (VSD, my_plot)
end