import Base: getindex, iterate

using Plots: scatter, scatter!, plot, plot!, Shape, cgrad, Plot
using DelaunayTriangulation: triangulate, each_solid_triangle, triangle_vertices, get_point, convert_boundary_points_to_indices
using LinearAlgebra: qr



export
	ValuedSubdivision,
	refine!,
    visualize,
    input_points,
    output_values,
    polygons,
	fibre_function_cache,
	complete_polygons,
	incomplete_polygons,
	fibre_function_type

#VALUEDSUBDIVISION
mutable struct ValuedSubdivision
    FibreFunctionCache :: Vector{Tuple{Vector{Float64},Float64}}
	CompletePolygons :: Vector{Vector{Int64}} #Indices
	IncompletePolygons :: Vector{Vector{Int64}}
	FibreFunctionType :: Symbol

	function ValuedSubdivision(EP::EnumerativeProblem; 
	xlims::Vector = [-1,1], ylims::Vector = [-1,1], 
	fibre_function = x->n_real_solutions(x), 
	resolution = 1000,
	mesh_function = trihexagonal_mesh)

	VSD = new()

	(polygons,parameters) = mesh_function(;xlims = xlims, ylims = ylims, resolution = resolution)
	
    #Solve for all pts in the initial hexagonal lattice
    subdivision_solutions = solve(EP, parameters)

    length(subdivision_solutions) != length(parameters) && error("Did not solve for each parameter")

	#Create fibre function cache
    fibre_function_cache = []
    for i in eachindex(subdivision_solutions)
        push!(fibre_function_cache, (parameters[i], fibre_function(subdivision_solutions[i])))
    end

	#Checking whether the fibre function is continuous or discrete
	fibre_function_type = 0
	output_values = getindex.(fibre_function_cache, 2)
	fibre_function_value_tally = Dict()
	for v in output_values
		get!(fibre_function_value_tally, v, 0)
		fibre_function_value_tally[v] += 1
	end
	if 	length(keys(fibre_function_value_tally)) < 50 
		fibre_function_type = :discrete 
	else
		fibre_function_type = :continuous
	end

	#Now we scroll through polygons and establish completeness
	complete_polygons::Vector{Vector{Int}} = []
	incomplete_polygons::Vector{Vector{Int}} = []
	
	tol = 0.0
	if fibre_function_type == :continuous
		#determining a tolerance for establishing completeness for continuous fibre function
		polygon_difference_values::Vector{Float64} = []
		for p in polygons
			vertex_values = sort([fibre_function_cache[v][2] for v in p])
			push!(polygon_difference_values, vertex_values[end] - vertex_values[1])
		end
		polygon_difference_values = sort(polygon_difference_values)
		tol_index = Int(ceil(0.5*length(polygon_difference_values)))
		tol = polygon_difference_values[tol_index]
	end

	for p in polygons #determining completeness of polygons using computed tol
		vertex_fibre_function_values = sort([fibre_function_cache[v][2] for v in p])
		if (vertex_fibre_function_values[end] - vertex_fibre_function_values[1]) <= tol
			push!(complete_polygons, p)
		else
			push!(incomplete_polygons, p)
		end
	end

	VSD.FibreFunctionCache = fibre_function_cache
	VSD.CompletePolygons = complete_polygons
	VSD.IncompletePolygons = incomplete_polygons
	VSD.FibreFunctionType = fibre_function_type
	return VSD
	end
end

#GETTERS
fibre_function_cache(VSD::ValuedSubdivision) = VSD.FibreFunctionCache
complete_polygons(VSD::ValuedSubdivision) = VSD.CompletePolygons
incomplete_polygons(VSD::ValuedSubdivision) = VSD.IncompletePolygons
fibre_function_type(VSD::ValuedSubdivision) = VSD.FibreFunctionType
input_points(VSD::ValuedSubdivision) = getindex.(VSD.FibreFunctionCache, 1)
output_values(VSD::ValuedSubdivision) = getindex.(VSD.FibreFunctionCache, 2)

#SETTERS
function set_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
	VSD.CompletePolygons = P
end

function set_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
	VSD.IncompletePolygons = P
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

function push_to_fibre_function_cache!(VSD::ValuedSubdivision, V::Tuple{Vector{Float64},Float64})
	push!(VSD.FibreFunctionCache, V)
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

#RESTRICTING EPS
#Note that these functions will have to be modified at some point so that a new EP isn't generated
function gram_schmidt(basis_vectors::Vector)
	M = hcat(basis_vectors...)
	Q = Matrix(qr(M).Q)
	return(collect(eachcol(Q)))
end



#If the user doesn't care WHICH plane to restrict to, do one through the current base
# parameters so that EP inherits the solutions
function planar_restriction(EP::EnumerativeProblem)
	P = [randn(Float64,n_parameters(EP)) for i in 1:3]
	return(restrict(EP,P))
end

#PLOTTING
function draw_polygon(polygon, color_value, VSD::ValuedSubdivision; label = false, 
	label_text = "real solutions", kwargs...)
	polygon_1 = [fibre_function_cache(VSD)[i][1] for i in polygon]
	polygon_1 = Shape([(t[1], t[2]) for t in polygon_1])
	c = cgrad(:thermal, rev = false)[color_value]
	if label == true
		plot!(polygon_1, fillcolor = c, linecolor = c, linewidth = false, label = label_text; kwargs...)
	else
		plot!(polygon_1, fillcolor = c, linecolor = c, linewidth = false, label = false; kwargs...)
	end
end

function draw_valued_subdivision(VSD::ValuedSubdivision; xlims = [-1,1], ylims = [-1,1], plot_log_transform = false, kwargs...) 
	my_plot = plot(xlims = xlims, ylims = ylims, aspect_ratio = :equal, background_color_inside=:black; kwargs...)
	if fibre_function_type(VSD) == :discrete
		plotting_values = unique(output_values(VSD))
		plotting_values = sort(plotting_values)
		values_that_have_been_plotted = []
		for i in plotting_values
			current_polygons = filter(x->fibre_function_cache(VSD)[x[1]][2] == i, complete_polygons(VSD)) #filters all polygons whose first vertices have the fibre function value i
			color_value = findfirst(x->x == i, plotting_values)/length(plotting_values)
			for j in current_polygons
					if (i in values_that_have_been_plotted) == false
						draw_polygon(j, color_value, VSD; label = true, label_text = "$i real solutions", kwargs...)
						push!(values_that_have_been_plotted, i)
					else
						draw_polygon(j, color_value, VSD; label = false, kwargs...)
					end
			end
		end
		return my_plot
	else
		polygons = vcat(complete_polygons(VSD), incomplete_polygons(VSD))
		values = []
		for p in polygons
			if plot_log_transform == true
				vertex_values = [log(fibre_function_cache(VSD)[v][2]) for v in p]
			else
				vertex_values = [fibre_function_cache(VSD)[v][2] for v in p]
			end
			polygon_value = sum(vertex_values)/length(p)
			push!(values, polygon_value)
		end
		
		max_value = max(values...)
		min_value = min(values...)
		for i in eachindex(polygons)
			color_value = values[i]/(max_value - min_value)
			draw_triangle(polygons[i], color_value, fibre_function_cache(VSD))
		end
		scatter!([0.0], [0.0]; zcolor = [min_value, max_value], color = cgrad(:thermal, rev = false), markersize = 0, colorbar = true, label = false)
		return my_plot
	end
end

#REFINEMENT FUNCTIONS

function global_delaunay_refine!(VSD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; 
	fibre_function = x->n_real_solutions(x), insertion_method = sierpinski_point_insertion, kwargs...)
	new_parameters_to_solve::Vector{Vector{Float64}} = []
	resolution_used = 0
	for T in incomplete_polygons(VSD)
		T_parameters = [graph_mesh(VSD)[v][1] for v in T]
		new_params, new_polygons = insertion_method(T_parameters)
		if length(new_params) > 0
			for p in new_params
				if !(p in new_parameters_to_solve) && !(p in input_points(graph_mesh(VSD)))
					push!(new_parameters_to_solve, p)
					resolution_used += 1
				end
			end
			resolution_used >= resolution && break
		else
			continue
		end
	end
	length(new_parameters_to_solve) == 0 && return resolution_used
	inserted_point_solutions = solve(EP, new_parameters_to_solve)
	length(inserted_point_solutions) != length(new_parameters_to_solve) && error("Did not solve for each parameter")
	for i in eachindex(inserted_point_solutions)
		push_to_graph_mesh!(graph_mesh(VSD), (new_parameters_to_solve[i], Float64(fibre_function(inserted_point_solutions[i]))))
	end
	retriangulate_valued_subdivision!(VSD)
	return resolution_used
end

function refine!(VSD::ValuedSubdivision, EP::EnumerativeProblem, resolution::Int64; 
	fibre_function = x->n_real_solutions(x),
	strategy = :quadtree)
	#Error checking to make sure that the refinement strategy will work with the VSD (e.g. sierpinski refinement won't work on a rectangular mesh)
	local_refinement_method = 0
	if strategy == :quadtree
		if length(complete_polygons(VSD)[1]) == 4
			local_refinement_method = quadtree_insertion
		else
			error("Cannot do quadtree refinement since polygons are not quadrilaterals.")
		end
	elseif strategy == :barycentric
		if length(complete_polygons(VSD)[1]) == 3
			local_refinement_method = barycenter_point_insertion
		else
			error("Cannot do barycentric refinement since polygons are not triangles.")
		end
	elseif strategy == :sierpinski
		if length(complete_polygons(VSD)[1]) == 3
			local_refinement_method = sierpinski_point_insertion
		else
			error("Cannot do sierpinski refinement since polygons are not triangles.")
		end
	elseif strategy == :random
		if length(complete_polygons(VSD)[1]) == 3
			local_refinement_method = random_point_insertion
		else
			error("Cannot do random refinement since polygons are not triangles.")
		end
	else
		error("Invalid strategy inputted. Valid strategies include
				:quadtree
				:barycentric
				:sierpinski
				:random")
	end

	refined_polygons::Vector{Vector{Int64}} = []
	polygons_to_solve_and_sort::Vector{Vector{Vector{Float64}}} = []
	new_parameters_to_solve::Vector{Vector{Float64}} = []
	resolution_used = 0
	for P in incomplete_polygons(VSD)
		P_parameters = [fibre_function_cache(VSD)[v][1] for v in P]
		new_params, new_polygons = local_refinement_method(P_parameters)
		if length(new_params) > 0
			push!(refined_polygons, P)
			push!(polygons_to_solve_and_sort, new_polygons...)
			for p in new_params
				if !(p in new_parameters_to_solve) && !(p in input_points(VSD))
					push!(new_parameters_to_solve, p)
					resolution_used += 1
				end
			end
			resolution_used >= resolution && break
		else
			continue
		end
	end
	length(new_parameters_to_solve) == 0 && return resolution_used
	new_parameters_solutions = solve(EP, new_parameters_to_solve)
	length(new_parameters_solutions) != length(new_parameters_to_solve) && error("Did not solve for each parameter")
	for i in eachindex(new_parameters_solutions)
		push_to_fibre_function_cache!(VSD, (new_parameters_to_solve[i], Float64(fibre_function(new_parameters_solutions[i]))))
	end
	for P in refined_polygons
		delete_from_incomplete_polygons!(VSD, P)
	end
	refinement_tol = 0.0
	if fibre_function_type(VSD) == :continuous
		polygons = vcat(complete_polygons(VSD), incomplete_polygons(VSD)) #this does not include the newly created incomplete rectangles that have been added during the current stage of refinement. Perhaps it should
		refinement_tol = continuous_fibre_function_refinement_threshold(VSD, polygons)
	end
	for P in polygons_to_solve_and_sort
		polygon = [findfirst(x->x[1]==y, fibre_function_cache(VSD)) for y in P]
		is_complete(polygon, VSD; tol = refinement_tol) ? push_to_complete_polygons!(VSD, polygon) : push_to_incomplete_polygons!(VSD, polygon)
	end
	return resolution_used
end

function quadtree_insertion(P::Vector{Vector{Float64}})
	center = [sum([x[1] for x in P])/4, sum([x[2] for x in P])/4]
	midpoint1 = midpoint([P[1], P[2]])
	midpoint2 = midpoint([P[2], P[3]])
	midpoint3 = midpoint([P[3], P[4]])
	midpoint4 = midpoint([P[4], P[1]])
	rectangle_1 = [P[1], midpoint1, center, midpoint4]
	rectangle_2 = [P[2], midpoint2, center, midpoint1]
	rectangle_3 = [P[3], midpoint3, center, midpoint2]
	rectangle_4 = [P[4], midpoint4, center, midpoint3]
	return [midpoint1, midpoint2, midpoint3, midpoint4, center], [rectangle_1, rectangle_2, rectangle_3, rectangle_4]
end

function random_point_insertion(P::Vector{Vector{Float64}})
	c = randn(Float64, length(P))
	c = map(abs, c)
	c = c./sum(c)
	random_convex_combination = [P[i].*c[i] for i in eachindex(P)]
	random_convex_combination = sum(random_convex_combination)
	triangle_1 = [P[1], P[2], random_convex_combination]
	triangle_2 = [P[2], P[3], random_convex_combination]
	triangle_3 = [P[3], P[1], random_convex_combination]
	return [random_convex_combination], [triangle_1, triangle_2, triangle_3]
end

function sierpinski_point_insertion(P::Vector{Vector{Float64}})
	midpoint_1 = midpoint([P[1], P[2]])
	midpoint_2 = midpoint([P[1], P[3]])
	midpoint_3 = midpoint([P[2], P[3]])
	triangle_1 = [P[1], midpoint_1, midpoint_2]
	triangle_2 = [P[2], midpoint_3, midpoint_1]
	triangle_3 = [P[3], midpoint_2, midpoint_3]
	return [midpoint_1, midpoint_2, midpoint_3], [triangle_1, triangle_2, triangle_3]
end

function barycenter_point_insertion(P::Vector{Vector{Float64}})
	barycenter = midpoint([v for v in P])
	triangle_1 = [P[1], P[2], barycenter]
	triangle_2 = [P[2], P[3], barycenter]
	triangle_3 = [P[3], P[1], barycenter]
	return [barycenter], [triangle_1, triangle_2, triangle_3]
end

function delaunay_retriangulate!(VSD::ValuedSubdivision)
	vertices = []
    for v in fibre_function_cache(VSD)
        v[1] in vertices || push!(vertices, v[1]) #ensuring that the triangulation will not contain duplicate points and the package won't give that annoying warning
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices) #This comes from DelaunayTriangulation.jl
	triangle_iterator = each_solid_triangle(tri) #Also from DelaunayTriangulation.jl
	complete_triangles::Vector{Vector{Int64}} = []
	incomplete_triangles::Vector{Vector{Int64}} = []
	if fibre_function_type(VSD) == :continuous
		triangles::Vector{Vector{Int64}} = []
		for T in triangle_iterator
			i,j,k = triangle_vertices(T)
			i,j,k = get_point(tri, i, j, k)
			vertex_1 = [i[1], i[2]]
			vertex_2 = [j[1], j[2]]
			vertex_3 = [k[1], k[2]]
			index_1 = findfirst(x->x[1] == vertex_1, fibre_function_cache(VSD))
			index_2 = findfirst(x->x[1] == vertex_2, fibre_function_cache(VSD))
			index_3 = findfirst(x->x[1] == vertex_3, fibre_function_cache(VSD))
			triangle = [index_1,index_2,index_3]
			push!(triangles, triangle)
		end
		refinement_tol = continuous_fibre_function_refinement_threshold(VSD, triangles)
		for t in triangles
			is_complete(t, VSD; tol = refinement_tol) ? push!(complete_triangles, t) : push!(incomplete_triangles, t)
		end
	else
		for T in triangle_iterator
			i,j,k = triangle_vertices(T)
			i,j,k = get_point(tri, i, j, k)
			vertex_1 = [i[1], i[2]]
			vertex_2 = [j[1], j[2]]
			vertex_3 = [k[1], k[2]]
			index_1 = findfirst(x->x[1] == vertex_1, fibre_function_cache(VSD))
			index_2 = findfirst(x->x[1] == vertex_2, fibre_function_cache(VSD))
			index_3 = findfirst(x->x[1] == vertex_3, fibre_function_cache(VSD))
			triangle = [index_1,index_2,index_3]
			is_complete(triangle, VSD) ? push!(complete_triangles, triangle) : push!(incomplete_triangles, triangle)
		end
	end
	set_complete_polygons!(VSD, complete_triangles)
	set_incomplete_polygons!(VSD, incomplete_triangles)
	return VSD
end

function midpoint(P::Vector{Vector{Float64}})
    return((sum(P))./length(P))
end

#checks if the max difference between vertex fibre function values is less than some tolerance. If so returns true, 
#otherwise returns false
function is_complete(p::Vector{Int}, VSD::ValuedSubdivision; tol = 0.0) 
	vertex_function_values = sort([fibre_function_cache(VSD)[x][2] for x in p])
	vertex_function_values[end] - vertex_function_values[1] <= tol ? (return true) : (return false)
end

#CONTINUOUS FIBRE FUNCTION VISUALIZATION
function is_continuous(VSD::ValuedSubdivision)
	value_tally = Dict()
	for v in output_values(VSD)
		get!(value_tally, v, 0)
		value_tally[v] += 1
	end
	println("n distinct values: ", length(keys(value_tally)))
	if 	length(keys(value_tally)) < 50 
		println("Not continuous")
		return false 
	else
		println("Continuous")
		 (return true) 
	end
	#if there are less than 50 unique values, then the fibre function is discrete
end

#determines a threshold for which polygons should be refined based on the difference in fibre function values across each polygon
function continuous_fibre_function_refinement_threshold(VSD::ValuedSubdivision, polygons::Vector{Vector{Int}})
	polygon_difference_values::Vector{Float64} = []
	for p in polygons
		vertex_values = sort([fibre_function_cache(VSD)[v][2] for v in p])
		push!(polygon_difference_values, vertex_values[end] - vertex_values[1])
	end
	polygon_difference_values = sort(polygon_difference_values)
	threshold_index = Int(ceil(0.5*length(polygon_difference_values)))
	return polygon_difference_values[threshold_index]
end

#VISUALIZATION
function visualize_fibre_function_cache(VSD::ValuedSubdivision)
    scatter(first.(input_points(VSD)), last.(input_points(VSD)), 
    zcolor = output_values(VSD), legend = false, colorbar = true)
end

function visualize(VSD::ValuedSubdivision;
	xlims = [-1.0, 1.0], ylims = [-1.0, 1.0],
	plot_log_transform = false,
	kwargs...)
	my_plot = draw_valued_subdivision(VSD; xlims = xlims, ylims = ylims, plot_log_tranform = plot_log_transform)
	return my_plot
end

#=
#TODO: Really use kwargs... here for passing xlims/ylims to initialize valued subdivision
function visualizer(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 4*initial_resolution,  
	fibre_function = x-> n_real_solutions(x), refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, local_insertion_function = quadtree_insertion,
	kwargs...)

	#Check enumerative problem is visualizable
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end

	VSD = initialize_valued_subdivision(new_EP; xlims = xlims, ylims = ylims, 
	fibre_function = fibre_function, initial_resolution = initial_resolution, mesh_function = mesh_function, kwargs...)

	remaining_resolution = total_resolution - initial_resolution
	while remaining_resolution > 0
		resolution_used = refine_function(VSD, new_EP, remaining_resolution; fibre_function = fibre_function, local_refinement_method = local_insertion_function, kwargs...)
		if resolution_used > 0
			remaining_resolution -= resolution_used
		else
			break
		end
	end

	my_plot = draw_valued_subdivision(VSD; xlims = xlims, ylims = ylims, kwargs...)
	return (VSD, my_plot)
end

function visualize(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 4*initial_resolution,  
	fibre_function = x-> n_real_solutions(x), refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, strategy = nothing,
	kwargs...)
	if isnothing(strategy)
		return visualizer(EP::EnumerativeProblem; xlims=xlims, ylims = ylims, 
		initial_resolution = initial_resolution, total_resolution = total_resolution,  
		fibre_function = fibre_function, refine_function = refine_function, mesh_function = mesh_function,
		kwargs...)
	elseif strategy == :quadtree
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		fibre_function = fibre_function, refine_function = locally_refine!, mesh_function = rectangular_mesh, local_insertion_function = quadtree_insertion, kwargs...)
	elseif strategy == :barycentric
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		fibre_function = fibre_function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh,
		local_insertion_function = barycenter_point_insertion, kwargs...)
	elseif strategy == :sierpinski
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		fibre_function = fibre_function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh, 
		local_insertion_function = sierpinski_point_insertion, kwargs...)
	elseif strategy == :random
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		fibre_function = fibre_function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh, local_insertion_function = random_point_insertion, kwargs...)
	elseif strategy == :delaunay||strategy == :Delaunay
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		fibre_function = fibre_function, refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, kwargs...)
	else
		println("Invalid strategy entered. Valid strategies include: ")
		println(":quadtree")
		println(":barycentric")
		println(":sierpinski")
		println(":random")
		println(":delaunay")
		return nothing
	end
end
=#

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
				plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = "$function_value_for_component")
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
					plot!(shape_1, fillcolor = c, linecolor = c, linewidth = false, label = "$function_value_for_component")
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

