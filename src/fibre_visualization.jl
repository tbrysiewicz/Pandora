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
    function_oracle,
    function_cache,
    complete_polygons,
    incomplete_polygons,
    is_discrete,
    delaunay_retriangulate!

# VALUEDSUBDIVISION
mutable struct ValuedSubdivision
    function_oracle::Function                           #This takes MANY parameters and returns a vector of values
    function_cache::Vector{Tuple{Vector{Float64},Any}}
    complete_polygons::Vector{Vector{Int64}}            # Given as indices of function_cache
    incomplete_polygons::Vector{Vector{Int64}}
    is_discrete::Bool
    is_complete::Function

    function ValuedSubdivision(function_oracle::Function; kwargs...)
        #Pull kwargs
        xlims = get(kwargs, :xlims, [-1,1])
        ylims = get(kwargs, :ylims, [-1,1])
        resolution = get(kwargs, :resolution, 1000)
        strategy = get(kwargs, :strategy, :sierpinski)


        # Determine default mesh's based on strategy
        if strategy == :quadtree
            mesh_function = rectangular_mesh
        elseif strategy == :barycentric || strategy == :sierpinski || strategy == :random
            mesh_function = trihexagonal_mesh
        end
    
        VSD = new()

        # Generate initial mesh
        (polygons, parameters) = mesh_function(; xlims=xlims, ylims=ylims, resolution=resolution, kwargs...)

        # Solve for all pts in the initial mesh
        function_oracle_values = function_oracle(parameters)
        length(function_oracle_values) != length(parameters) && error("Did not solve for each parameter")

        # Create function cache
        function_cache = Vector{Tuple{Vector{Float64},Any}}([])
        for i in eachindex(function_oracle_values)
            push!(function_cache, (parameters[i], function_oracle_values[i]))
        end

        # Checking whether the function is continuous or discrete
        output_values = getindex.(function_cache, 2)
        unique_count = length(unique(output_values))
        is_discrete = unique_count < 50

        # Now we scroll through polygons and establish completeness
        complete_polygons::Vector{Vector{Int}} = []
        incomplete_polygons::Vector{Vector{Int}} = []

        tol = 0.0
        if !is_discrete
            polygon_diffs = [max([function_cache[v][2] for v in p]...) - min([function_cache[v][2] for v in p]...) for p in polygons]
            tol = median(polygon_diffs)
        end

        local_is_complete = (p::Vector{Int}, FC:: Vector{Tuple{Vector{Float64},Any}}) -> is_complete(p, FC; tol=tol)
        if haskey(kwargs, :is_complete)
            ic = kwargs[:is_complete]
            local_is_complete = (p::Vector{Int}, FC:: Vector{Tuple{Vector{Float64},Any}};kwargs...) -> ic(p, FC;kwargs...) # If a custom is_complete function is provided, use it
        else
            local_is_complete = (p::Vector{Int}, FC:: Vector{Tuple{Vector{Float64},Any}};kwargs...) -> is_complete(p, FC; tol=tol, kwargs...) # Default function to determine completeness of polygons
        end

        for p in polygons # determining completeness of polygons using computed tol
            if local_is_complete(p,function_cache)
                push!(complete_polygons, p)
            else
                push!(incomplete_polygons, p)
            end
        end

        VSD.function_oracle = function_oracle
        VSD.function_cache = function_cache
        VSD.complete_polygons = complete_polygons
        VSD.incomplete_polygons = incomplete_polygons
        VSD.is_discrete = is_discrete
        VSD.is_complete = local_is_complete
        return VSD
    end

    function ValuedSubdivision(EP::EnumerativeProblem; kwargs...)
        fibre_function = get(kwargs, :fibre_function, x->valid_real_solution_set(EP,x) ? n_real_solutions(x) : :wild)
        function_oracle = get(kwargs, :function_oracle, x->fibre_function.(EP(x)))

        ValuedSubdivision(function_oracle; kwargs...)
    end
end

# GETTERS
function_oracle(VSD::ValuedSubdivision) = VSD.function_oracle
function_cache(VSD::ValuedSubdivision) = VSD.function_cache
complete_polygons(VSD::ValuedSubdivision) = VSD.complete_polygons
incomplete_polygons(VSD::ValuedSubdivision) = VSD.incomplete_polygons
is_discrete(VSD::ValuedSubdivision) = VSD.is_discrete
input_points(FC::Vector{Tuple{Vector{Float64}, Any}}) = getindex.(FC, 1)
output_values(FC::Vector{Tuple{Vector{Float64}, Any}}) = getindex.(FC, 2)
input_points(VSD::ValuedSubdivision) = input_points(function_cache(VSD))
output_values(VSD::ValuedSubdivision) = output_values(function_cache(VSD))


function establish_completeness!(VSD::ValuedSubdivision)

end

# SETTERS
#TODO: I think you can eliminate the next three functions
function set_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
    VSD.complete_polygons = P
end

function set_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Vector{Int64}})
    VSD.incomplete_polygons = P
end

function delete_from_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    polygon_index = findfirst(x -> x == P, incomplete_polygons(VSD))
    if polygon_index !== nothing
        deleteat!(VSD.incomplete_polygons, polygon_index)
    else
        error("Polygon not found in IncompletePolygons!")
    end
end

#TODO: Probably unnecessary wrappers for push
function push_to_complete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    push!(VSD.complete_polygons, P)
end

function push_to_incomplete_polygons!(VSD::ValuedSubdivision, P::Vector{Int64}; kwargs...)
    push!(VSD.incomplete_polygons, P)
end

function push_to_function_cache!(VSD::ValuedSubdivision, V::Tuple{Vector{Float64},Any}; kwargs...)
    push!(VSD.function_cache, V)
end

#region MESH FUNCTIONS

function initial_parameter_distribution(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)
    xlength = xlims[2] - xlims[1]
    ylength = ylims[2] - ylims[1]

    x_values = range(xlims[1], xlims[2], Int(floor(sqrt((((xlength)/(xlength + ylength))*resolution)/(1 - ((xlength)/(xlength + ylength)))))))
    y_values = range(ylims[1], ylims[2], Int(floor(sqrt((((ylength)/(xlength + ylength))*resolution)/(1 - ((ylength)/(xlength + ylength)))))))

    return x_values, y_values
end


function trihexagonal_mesh(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)

    x_values, y_values = initial_parameter_distribution(; xlims=xlims, ylims=ylims, resolution=resolution)

    shift_amount = (x_values[2] - x_values[1])/2 	
    parameters = [[i - isodd(findfirst(x->x==j, y_values))*shift_amount, j] for i in x_values for j in y_values]

    parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))

    triangles = Vector{Vector{Int}}([])
    for i in 1:length(x_values)-1
        for j in 1:length(y_values)-1
            tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
            tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
            bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
            br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)

            if isodd(i)
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

function rectangular_mesh(; kwargs...)
    xlims = get(kwargs, :xlims, [-1,1])
    ylims = get(kwargs, :ylims, [-1,1])
    resolution = get(kwargs, :resolution, 1000)
    x_values, y_values = initial_parameter_distribution(; xlims=xlims, ylims=ylims, resolution=resolution)
    parameters = [[i,j] for i in x_values for j in y_values]
    parameters_as_matrix = reshape(parameters, length(x_values), length(y_values))
    rectangles = Vector{Vector{Int}}([])
    for i in 1:length(x_values)-1
        for j in 1:length(y_values)-1
            tl_index = findfirst(x->x == parameters_as_matrix[i,j], parameters)
            tr_index = findfirst(x->x == parameters_as_matrix[i+1,j], parameters)
            bl_index = findfirst(x->x == parameters_as_matrix[i, j+1], parameters)
            br_index = findfirst(x->x == parameters_as_matrix[i+1, j+1], parameters)
            push!(rectangles, [tl_index, tr_index, br_index, bl_index])
        end
    end
    return (rectangles, parameters)
end



#REFINEMENT FUNCTIONS


function refine!(VSD::ValuedSubdivision, resolution::Int64;	strategy = nothing)
	local_refinement_method = 0
    if strategy === nothing
        if length(complete_polygons(VSD)[1]) == 4
            strategy = :quadtree
        elseif length(complete_polygons(VSD)[1]) == 3
            strategy = :sierpinski
        end
    end
	if strategy == :quadtree && length(complete_polygons(VSD)[1]) != 4
		error("Cannot do this refinement since polygons are not quadrilaterals.")
	elseif (strategy == :barycentric || strategy == :sierpinski || strategy == :random) && length(complete_polygons(VSD)[1]) != 3
		error("Cannot do this refinement since polygons are not triangles.")
	elseif in(strategy, [:quadtree, :barycentric, :sierpinski, :random]) == false
		error("Invalid strategy inputted. Valid strategies include
				:quadtree
				:barycentric
				:sierpinski
				:random")
	end
	if strategy == :quadtree
		local_refinement_method = quadtree_insertion
	elseif strategy == :barycentric
		local_refinement_method = barycenter_point_insertion
	elseif strategy == :sierpinski
		local_refinement_method = sierpinski_point_insertion
	elseif strategy == :random
		local_refinement_method = random_point_insertion
	end
    FO = function_oracle(VSD)
	refined_polygons::Vector{Vector{Int64}} = []
	polygons_to_solve_and_sort::Vector{Vector{Vector{Float64}}} = []
	new_parameters_to_solve::Vector{Vector{Float64}} = []
	resolution_used = 0
	for P in incomplete_polygons(VSD)
		P_parameters = [function_cache(VSD)[v][1] for v in P]
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
	function_oracle_values = FO(new_parameters_to_solve)
	length(function_oracle_values) != length(new_parameters_to_solve) && error("Did not solve for each parameter")
	for i in eachindex(function_oracle_values)
         push_to_function_cache!(VSD, (new_parameters_to_solve[i], (function_oracle_values[i])))
	end
	for P in refined_polygons
		delete_from_incomplete_polygons!(VSD, P)
	end
	refinement_tol = 0.0
	if is_discrete(VSD) == false
		polygons = vcat(complete_polygons(VSD), incomplete_polygons(VSD)) #this does not include the newly created incomplete rectangles that have been added during the current stage of refinement. Perhaps it should
		refinement_tol = continuous_function_refinement_threshold(VSD, polygons)
	end
	for P in polygons_to_solve_and_sort
		polygon = [findfirst(x->x[1]==y, function_cache(VSD)) for y in P]
		is_complete(polygon, VSD; tol = refinement_tol) ? push_to_complete_polygons!(VSD, polygon) : push_to_incomplete_polygons!(VSD, polygon)
	end
	return resolution_used
end



function refine!(VSD::ValuedSubdivision; kwargs...)
    # Default refinement strategy is quadtree
    return refine!(VSD, 1000000; kwargs...)
end

function quadtree_insertion(P::Vector{Vector{Float64}})
    center = midpoint(P)
    midpoints = [midpoint([P[i], P[mod1(i+1, 4)]]) for i in 1:4]
    rectangles = [[P[i], midpoints[i], center, midpoints[mod1(i-1, 4)]] for i in 1:4]
    return (midpoints..., center), rectangles
end

function random_point_insertion(P::Vector{Vector{Float64}})
    c = abs.(randn(Float64, length(P)))
    c ./= sum(c)
    random_point = sum(P[i] .* c[i] for i in eachindex(P))
    triangles = [[P[1], P[2], random_point],
                 [P[2], P[3], random_point],
                 [P[3], P[1], random_point]]
    return [random_point], triangles
end

function sierpinski_point_insertion(P::Vector{Vector{Float64}})
    m1 = midpoint([P[1], P[2]])
    m2 = midpoint([P[1], P[3]])
    m3 = midpoint([P[2], P[3]])
    triangles = [[P[1], m1, m2],
                 [P[2], m3, m1],
                 [P[3], m2, m3],
                 [m1, m2, m3]]
    return [m1, m2, m3], triangles
end

function barycenter_point_insertion(P::Vector{Vector{Float64}})
    barycenter = midpoint(P)
    triangles = [[P[1], P[2], barycenter],
                 [P[2], P[3], barycenter],
                 [P[3], P[1], barycenter]]
    return [barycenter], triangles
end

function delaunay_retriangulate!(VSD::ValuedSubdivision)
	vertices = []
    for v in function_cache(VSD)
        v[1] in vertices || push!(vertices, v[1]) #ensuring that the triangulation will not contain duplicate points and the package won't give that annoying warning
    end
	vertices = hcat(vertices...)
	tri = triangulate(vertices) #This comes from DelaunayTriangulation.jl
	triangle_iterator = each_solid_triangle(tri) #Also from DelaunayTriangulation.jl
	complete_triangles::Vector{Vector{Int64}} = []
	incomplete_triangles::Vector{Vector{Int64}} = []
	if is_discrete(VSD) == false
		triangles::Vector{Vector{Int64}} = []
		for T in triangle_iterator
			i,j,k = triangle_vertices(T)
			i,j,k = get_point(tri, i, j, k)
			vertex_1 = [i[1], i[2]]
			vertex_2 = [j[1], j[2]]
			vertex_3 = [k[1], k[2]]
			index_1 = findfirst(x->x[1] == vertex_1, function_cache(VSD))
			index_2 = findfirst(x->x[1] == vertex_2, function_cache(VSD))
			index_3 = findfirst(x->x[1] == vertex_3, function_cache(VSD))
			triangle = [index_1,index_2,index_3]
			push!(triangles, triangle)
		end
		refinement_tol = continuous_function_refinement_threshold(VSD, triangles)
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
			index_1 = findfirst(x->x[1] == vertex_1, function_cache(VSD))
			index_2 = findfirst(x->x[1] == vertex_2, function_cache(VSD))
			index_3 = findfirst(x->x[1] == vertex_3, function_cache(VSD))
			triangle = [index_1,index_2,index_3]
			is_complete(triangle, VSD) ? push!(complete_triangles, triangle) : push!(incomplete_triangles, triangle)
		end
	end
	set_complete_polygons!(VSD, complete_triangles)
	set_incomplete_polygons!(VSD, incomplete_triangles)
	return VSD
end

function midpoint(P::Vector{T}) where T <: Any
    if length(P) == 0
        return nothing
    end
    return((sum(P))./length(P))
end

mean(P::Vector{T}) where T <: Any = midpoint(P)

function median(V::AbstractVector)
    Vsorted = sort(V)
    n = length(Vsorted)
    if n == 0
        error("Cannot compute median of empty vector")
    elseif isodd(n)
        return Vsorted[(n+1) ÷ 2]
    else
        return (Vsorted[n ÷ 2] + Vsorted[n ÷ 2 + 1]) / 2
    end
end

#TODO: Abstract this function as a kwarg (later)
function is_complete(polygon::Vector{Int}, FC::Vector{Tuple{Vector{Float64},Any}}; tol = 0.0) 
    vals = [FC[v][2] for v in polygon]
    vertex_function_values = sort(filter(x->isa(x,Number),vals))
    if length(vertex_function_values) == 0
        return false # If there are no vertices, we consider it incomplete
    end
    if (vertex_function_values[end] - vertex_function_values[1]) <= tol
        return true
    else
        return false
    end
end

function is_complete(polygon::Vector{Int}, VSD::ValuedSubdivision; kwargs...) 
    VSD.is_complete(polygon, function_cache(VSD); kwargs...)
end


#determines a threshold for which polygons should be refined based on the difference in function values across each polygon
function continuous_function_refinement_threshold(VSD::ValuedSubdivision, polygons::Vector{Vector{Int}})
	polygon_difference_values::Vector{Float64} = []
	for p in polygons
		vertex_values = sort([function_cache(VSD)[v][2] for v in p])
		push!(polygon_difference_values, vertex_values[end] - vertex_values[1])
	end
	polygon_difference_values = sort(polygon_difference_values)
	threshold_index = Int(ceil(0.5*length(polygon_difference_values)))
	return polygon_difference_values[threshold_index]
end

#VISUALIZATION
function visualize_function_cache(VSD::ValuedSubdivision)
    visualize(function_cache(VSD))
end
function visualize(FC::Vector{Tuple{Vector{Float64},Any}})
    scatter(first.(input_points(FC)), last.(input_points(FC)), 
    zcolor = map(x->isa(x,Number) ? x : -10.0, output_values(FC)), legend = false, colorbar = true)
end


function visualize(VSD::ValuedSubdivision; kwargs...)
    xl = get(kwargs, :xlims, [min(map(x->x[1][1],Pandora.function_cache(VSD))...),max(map(x->x[1][1],Pandora.function_cache(VSD))...)])
    yl = get(kwargs, :ylims, [min(map(x->x[1][2],Pandora.function_cache(VSD))...),max(map(x->x[1][2],Pandora.function_cache(VSD))...)])
    plot_log_transform = get(kwargs, :plot_log_transform, false)
    plot_all_polygons = get(kwargs, :plot_all_polygons, is_discrete(VSD) == false)


    my_plot = plot(xlims = xl, ylims = yl, aspect_ratio = :equal, background_color_inside=:black; kwargs...)
    
    polygon_list = plot_all_polygons == true ? vcat(complete_polygons(VSD), incomplete_polygons(VSD)) : complete_polygons(VSD)
    polygons_to_draw = map(P->map(first,function_cache(VSD)[P]), polygon_list)
    # we remove non-numbers when computing means so that non-numbers function as wild cards essentially 
    polygon_values = map(P-> mean(filter(x->isa(x,Number),map(last,function_cache(VSD)[P]))), polygon_list)
    if plot_log_transform
        # If we are plotting log transformed values, we need to transform the polygon values
        polygon_values = map(x->x==nothing ? nothing : log(x), polygon_values)
    end
    # It is still possible that the mean returns 'nothing' if all values are non-numbers
    real_values = filter(x->x!=nothing, unique(polygon_values))
    max_value = max(real_values...)
    shapes = Shape[]
    colors = []

    for (poly_pts, value) in zip(polygons_to_draw, polygon_values)
        push!(shapes, Shape([(t[1], t[2]) for t in poly_pts]))
        push!(colors, value==nothing ? :white : cgrad(:thermal, rev=false)[value/max_value])
    end
    MyPlot = plot!(shapes; fillcolor=permutedims(colors), linecolor=permutedims(colors), linewidth=0, label = false)

    
    legend_values = sort(real_values, rev=true)
    if in(nothing,polygon_values)
        pushfirst!(legend_values, nothing)
    end
    if length(legend_values) <= 10
        # Add legend entries for unique values
        for val in legend_values
            color = val==nothing ? :white : cgrad(:thermal, rev=false)[val/max_value]
            # Plot an invisible shape with the correct color and label for the legend
            plot!(Shape([(NaN, NaN), (NaN, NaN), (NaN, NaN)]); fillcolor=color, linecolor=:transparent, linewidth=0, label="$val", legend = :outerright)
        end
    else
        # Add a colorbar for the continuous values
        scatter!([NaN], [NaN]; zcolor = [min(real_values...), max(real_values...)], color = cgrad(:thermal, rev = false), markersize = 0, colorbar = true, label = false)
    end
    return(MyPlot)
end


#=
#TODO: Really use kwargs... here for passing xlims/ylims to initialize valued subdivision
function visualizer(EP::EnumerativeProblem; xlims=[-1,1], ylims = [-1,1], 
	initial_resolution = 1000,	total_resolution = 4*initial_resolution,  
	function = x-> n_real_solutions(x), refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, local_insertion_function = quadtree_insertion,
	kwargs...)

	#Check enumerative problem is visualizable
	if n_parameters(EP) > 2
		println("EP consists of more than two parameters. Visualizing a random 2-plane in the parameter space.")
		new_EP = planar_restriction(EP)
	else
		new_EP = EP
	end

	VSD = initialize_valued_subdivision(new_EP; xlims = xlims, ylims = ylims, 
	function = function, initial_resolution = initial_resolution, mesh_function = mesh_function, kwargs...)

	remaining_resolution = total_resolution - initial_resolution
	while remaining_resolution > 0
		resolution_used = refine_function(VSD, new_EP, remaining_resolution; function = function, local_refinement_method = local_insertion_function, kwargs...)
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
	function = x-> n_real_solutions(x), refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, strategy = nothing,
	kwargs...)
	if isnothing(strategy)
		return visualizer(EP::EnumerativeProblem; xlims=xlims, ylims = ylims, 
		initial_resolution = initial_resolution, total_resolution = total_resolution,  
		function = function, refine_function = refine_function, mesh_function = mesh_function,
		kwargs...)
	elseif strategy == :quadtree
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		function = function, refine_function = locally_refine!, mesh_function = rectangular_mesh, local_insertion_function = quadtree_insertion, kwargs...)
	elseif strategy == :barycentric
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		function = function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh,
		local_insertion_function = barycenter_point_insertion, kwargs...)
	elseif strategy == :sierpinski
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		function = function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh, 
		local_insertion_function = sierpinski_point_insertion, kwargs...)
	elseif strategy == :random
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		function = function, refine_function = locally_refine!, mesh_function = trihexagonal_mesh, local_insertion_function = random_point_insertion, kwargs...)
	elseif strategy == :delaunay||strategy == :Delaunay
		return visualizer(EP::EnumerativeProblem; xlims = xlims, ylims = ylims, initial_resolution = initial_resolution, total_resolution = total_resolution,
		function = function, refine_function = global_delaunay_refine!, mesh_function = trihexagonal_mesh, kwargs...)
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
