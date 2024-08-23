using DelaunayTriangulation

export
    restrict_enumerative_problem_to_plane,
    restrict_enumerative_problem,
    visualize_parameterspace







@doc raw"""
    restrict_enumerative_problem_to_plane(EP::EnumerativeProblem)

 Restrict an enumerative problem with >2 parameters to one with 2 parameters by slicing 
   the parameter space with a random 2-plane. This is often used to visualize portions
   of the discriminant of an enumerative problem. 

# Examples
```jldoctest

julia> T = TwentySevenLines()


           X := V(f_1..f_4) ⊆ C^4 x C^20
           |
           |
           | π    ???-to-1
           |
           V
          C^20

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 20 parameters.


julia> restrict_enumerative_problem_to_plane(T)

           X := V(f_1..f_4) ⊆ C^4 x C^2
           |
           |
           | π    ???-to-1
           |
           V
          C^2

An enumerative problem in 4 variable(s) cut out by 4 condition(s) over 2 parameters.

```
"""
function restrict_enumerative_problem_to_plane(EP::EnumerativeProblem)
	P = [randn(Float64,n_parameters(EP)) for i in 1:3]
	return(restrict_enumerative_problem(EP,P))
end

function restrict_enumerative_problem(EP::EnumerativeProblem,P::Vector{Vector{Float64}})
	F = system(EP)
	xv = variables(F)
	xp = parameters(F)
	n = length(P)
	@var t[1:n-1]
	affine_span = P[n]+sum([t[i].*(P[i]-P[n]) for i in 1:n-1])
	NewEquations = [subs(f,xp=>affine_span) for f in expressions(F)]
	return(EnumerativeProblem(System(NewEquations,variables=xv,parameters=t)))
end

mutable struct mesh
    EP::EnumerativeProblem
    fibre_function::Function
    value_dict::Dict
    triangulation::Vector{Vector{Vector{Float64}}}
    xlims::Vector{Float64}
    ylims::Vector{Float64}
    type_of_fibre_function::String
end

#=
function initialize_mesh
    input:
        EP - An enumerative problem
        fibre_function - fibre_function used to output values for each mesh vertex
        xlims - mesh bounds for x parameter
        ylims - mesh bounds for y parameter
        resolution - upper bound on number of vertices in mesh
        continuous - indicates whether the fibre_function is continuous or discrete
    output:
        mesh1 - mesh object
        length(P) - number of points in the mesh; used to keep track of total resolution used
=#
function initialize_mesh(EP::EnumerativeProblem; 
                        xlims = [-2.0, 2.0],
                        ylims = [-2.0, 2.0],
                        fibre_function = x->HomotopyContinuation.nreal(x[1]),
                        resolution = 1000,
                        continuous = false)

    nxy = Int(floor(sqrt(resolution))) #mesh points are distributed equally along x and y axis regardless of size of xlims and ylims producing "square" mesh
    xrange = range(xlims[1],xlims[2],nxy)
    yrange = range(ylims[1],ylims[2],nxy)
    delta = (xrange[2]-xrange[1])/2
    triangles = []
    P = [[i-isodd(findfirst(x->x==j,yrange))*delta,j] for i in xrange for j in yrange] #shifts x values to the left by amount delta in every other row in the mesh
    S = solve_over_params(EP,P; checks = [], retry=true)
    value_dict=Dict{Any,Any}()
    for s in S
        value_dict[s[2]]=fibre_function(s)
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
        	value_dict[s[2]] = false #parameters are given false value if they produce an "error"
        end
    end
    M = reshape(P,length(xrange),length(yrange))
    for i in 1:length(xrange)-1
        for j in 1:length(yrange)-1
            tl = M[i,j]
            tr = M[i+1,j]
            bl = M[i,j+1]
            br = M[i+1,j+1]
            if isodd(i)
                push!(triangles,[tl,tr,bl])
                push!(triangles,[bl,tr,br])
            else
                push!(triangles,[tl,tr,br])
                push!(triangles,[bl,tl,br])
            end
        end
    end
    if continuous == false
        mesh1 = mesh(EP, fibre_function, value_dict, triangles, xlims, ylims, "discrete")
    else
        mesh1 = mesh(EP, fibre_function, value_dict, triangles, xlims, ylims, "continuous")
    end
    return mesh1, length(P)
end

#=function draw_triangle
    input:
        T - triangle to be drawn consisting of a vector containing the three vertices
        v - the value that is used to select a color from the gradient for triangle fillcolor
        label - if label == true, a label is added in the legend for that triangle
        labelText - text for label in the case that label == true
    output: 
        none - a triangle is plotted in current plot
=#
function draw_triangle(T, v; 
                        label = false,
                        labelText = "hello")
    Triangle = Shape([(t[1],t[2]) for t in T])
    c =cgrad(:thermal, rev = false)[v]
	if label == true
		plot!(Triangle, fillcolor = c, linecolor = c, linewidth = false, labels = labelText)
	else
		plot!(Triangle, fillcolor = c, linecolor = c, linewidth = false, labels = false)
	end
end

#=function draw_mesh
    input:
        value_dict
        triangles
        xlims
        ylims
        continuous
        scatter - if scatter is true the triangulation vertices will be plotted
        label - if label is true, a legend will be added to the plot
        plot_logged_values - in the case of a continuous fibre_function if plot_logged_values is true, the logged fibre_function values will be plotted
    output:
        my_plot - plot of mesh
=#
function draw_mesh(value_dict::Dict, triangles::Vector;
                                xlims = [-2.0,2.0],
                                ylims=[-2.0,2.0], 
                                continuous = false, 
                                scatter = false, 
                                label = true, 
                                plot_logged_values = false)
    if continuous == true
		my_plot = plot(xlims=xlims, ylims=ylims, colorbar = true)
		triangle_plotting_values = Dict()
		vertices_to_plot = []
		for T in triangles
            if isa(value_dict[T[1]], Bool) || isa(value_dict[T[2]], Bool) || isa(value_dict[T[3]], Bool)
                non_error_coordinates = findall(x->isa(value_dict[x], Bool) == false, T)
                if length(non_error_coordinates) == 2
                    if plot_logged_values == false
                        triangle_plotting_values[T] = (value_dict[T[non_error_coordinates[1]]] + value_dict[T[non_error_coordinates[2]]])/2
                    else
                        triangle_plotting_values[T] = (log(value_dict[T[non_error_coordinates[1]]]) + log(value_dict[T[non_error_coordinates[2]]]))/2
                    end
                elseif length(non_error_coordinates) == 1
                    if plot_logged_values == false
                        triangle_plotting_values[T] = value_dict[T[non_error_coordinates[1]]]
                    else
                        triangle_plotting_values[T] = log(value_dict[T[non_error_coordinates[1]]])
                    end
                else
                    continue
                end
                for i in non_error_coordinates
                    push!(vertices_to_plot, T[i])
                end
                continue
            end
            if plot_logged_values == false
				triangle_plotting_values[T] = mean([value_dict[x] for x in T])
            else
                triangle_plotting_values[T] = mean([log(value_dict[x]) for x in T])
            end
			for i in T
				push!(vertices_to_plot, i)
			end
		end
        if plot_logged_values == false
			plotting_max = max([value_dict[x] for x in vertices_to_plot]...)
			plotting_min = min([value_dict[x] for x in vertices_to_plot]...)
        else
            plotting_max = max([log(value_dict[x]) for x in vertices_to_plot]...)
            plotting_min = min([log(value_dict[x]) for x in vertices_to_plot]...)
        end
		if scatter == true
            if plot_logged_values == false
				scatter!([a[1] for a in vertices_to_plot], [a[2] for a in vertices_to_plot], zcolor = [value_dict[a] for a in vertices_to_plot], color =:thermal, markershape =:rect, markerstrokewidth = 0.5, labels = false, markersize = 1)
            else
                scatter!([a[1] for a in vertices_to_plot], [a[2] for a in vertices_to_plot], zcolor = [log(value_dict[a]) for a in vertices_to_plot], color =:thermal, markershape =:rect, markerstrokewidth = 0.5, labels = false, markersize = 1) 
            end
		else
            if plot_logged_values == false
				scatter!([a[1] for a in vertices_to_plot], [a[2] for a in vertices_to_plot], zcolor = [value_dict[a] for a in vertices_to_plot], color =:thermal, markershape =:rect, markerstrokewidth = 0, labels = false, markersize = 0) #even if scatter == false, still need to plot "invisible" points to get a colorbar that can be used to generate fillcolors for plotting triangles
            else
                scatter!([a[1] for a in vertices_to_plot], [a[2] for a in vertices_to_plot], zcolor = [log(value_dict[a]) for a in vertices_to_plot], color =:thermal, markershape =:rect, markerstrokewidth = 0, labels = false, markersize = 0)
            end
		end
		for T in keys(triangle_plotting_values)
			color_value = (triangle_plotting_values[T]-plotting_min)/(plotting_max-plotting_min)
            draw_triangle(T, color_value)
		end
    else
        plotting_values = filter(x -> isa(x, Bool) == false, unique(values(value_dict)))
        sorted_plotting_values = sort(plotting_values) #used to assign color values to triangles
        number_of_unique_plotting_values = length(sorted_plotting_values)
        my_plot = plot(xlims=xlims,ylims=ylims,legend=true)
        error_triangles = filter(x-> isa(value_dict[x[1]], Bool) || isa(value_dict[x[2]], Bool) || isa(value_dict[x[3]], Bool), triangles)
        values_already_plotted = []
        for v in sorted_plotting_values
            current_triangles = filter(x->value_dict[x[1]]==v || value_dict[x[2]] == v || value_dict[x[3]] == v, triangles)
            for i in current_triangles
               if i in error_triangles
                    error_coordinates = findall(x->isa(value_dict[x],Bool), i)
                    if length(error_coordinates) == 1
                        non_error_coordinates = findall(x->isa(value_dict[x], Bool)==false, i)
                        if value_dict[i[non_error_coordinates[1]]] == value_dict[i[non_error_coordinates[2]]]
                            draw_triangle(i, findfirst(x->x==value_dict[i[non_error_coordinates[1]]], sorted_plotting_values)/number_of_unique_plotting_values, label = false)
                        else
                            continue
                        end
                    elseif length(error_coordinates) == 2
                        draw_triangle(i, findfirst(x->x == v, sorted_plotting_values)/number_of_unique_plotting_values, label = false)
                    elseif length(error_coordinates) == 3
                        continue
                    end
                end
                if value_dict[i[1]] != value_dict[i[2]] || value_dict[i[1]] != value_dict[i[3]] 
                    continue
                end
                if (value_dict[i[1]] in values_already_plotted) == false && label == true
                    draw_triangle(i, findfirst(x->x == v, sorted_plotting_values)/number_of_unique_plotting_values, label = true, labelText = "$(value_dict[i[1]]) real solutions")
                    push!(values_already_plotted, v)
                else
                    draw_triangle(i, findfirst(x->x == v, sorted_plotting_values)/number_of_unique_plotting_values, label = false)
                end
            end
            if scatter == true
                current_parameters = filter(x->value_dict[x] == v, keys(value_dict))
                c = cgrad(:thermal, rev = false)[findfirst(x->x==v, sorted_plotting_values)/number_of_unique_plotting_values]
                scatter!([A[1] for A in current_parameters], [A[2] for A in current_parameters], markercolor = c, markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, labels = false)
            end
        end
    end
    return(my_plot)
end

#= function triforce_refinement! - finds triangles in mesh over which the number of real solutions changes and inserts triforce points into these triangles then solves for new vertices.
    input:
        mesh1 - mesh object
        resolution
        standard_for_refinement - in the case that fibre_function is continuous, used as the threshold for triangle refinement; this value is produced by standard_for_continuous_function_refinement
    output:
        refinement_termination - indicates whether total resolution has been reached in current round of refinement
        length(S) - number of points solved for in refinement
=#
function triforce_refinement!(mesh1::mesh, resolution = 1000, standard_for_refinement = 0.0)
    value_dict = mesh1.value_dict
    triangles = mesh1.triangulation
    fibre_function = mesh1.fibre_function
	refinement_termination = false
    white_space_triangles = []
    triangles_to_refine = []
    if mesh1.type_of_fibre_function == "continuous"
        triangle_value_ranges = Dict()
        for T in triangles
            if isa(value_dict[T[1]], Bool) || isa(value_dict[T[2]], Bool) || isa(value_dict[T[3]], Bool)
                non_error_coordinates = findall(x->isa(value_dict[x], Bool) == false, T)
                if length(non_error_coordinates) == 2
                    triangle_value_ranges[T] = abs(value_dict[T[non_error_coordinates[1]]] - value_dict[T[non_error_coordinates[2]]])
                else
                    continue
                end
            else
                min_value = min([value_dict[x] for x in T]...)
                max_value = max([value_dict[x] for x in T]...)
                triangle_value_ranges[T] = max_value - min_value
            end
        end
        triangles_to_refine = filter(x->triangle_value_ranges[x]>standard_for_refinement, keys(triangle_value_ranges))
    else
        white_space_triangles = filter(x->value_dict[x[1]]!=value_dict[x[2]]||value_dict[x[1]]!=value_dict[x[3]], triangles)
        for T in white_space_triangles
            if isa(value_dict[T[1]], Bool) == false && isa(value_dict[T[2]], Bool) == false && isa(value_dict[T[3]], Bool) == false
                push!(triangles_to_refine, T)
            else
                non_error_coordinates = findall(x->isa(value_dict[x], Bool) == false, T)
                if length(non_error_coordinates) <= 1
                    continue
                else
                    if value_dict[T[non_error_coordinates[1]]] != value_dict[T[non_error_coordinates[2]]] #triangles with a single error and two differing values on the other two vertices are refined
                        push!(triangles_to_refine, T)
                    else
                        continue
                    end
                end
            end
        end
    end
    restricted_triangles_to_refine = [] #need to restrict triangles_to_refine to meet resolution requirement
	if (length(triangles_to_refine)*3)>resolution #triforce refinement inserts three points in each white space triangle. Therefore, if the length of triangles_to_refine times 3 is greater than the resolution, only resolution/3 many triangles are chosen to be refined
		while length(restricted_triangles_to_refine)*3<resolution
			v = rand(triangles_to_refine,Int(floor(resolution/10)))
			restricted_triangles_to_refine = unique(vcat(restricted_triangles_to_refine,v))
		end
		restricted_triangles_to_refine = restricted_triangles_to_refine[1:Int(floor(resolution/3))]
		refinement_termination = true #when triangles_to_refine needs to be restricted, this indicates that total_resolution has been reached and refinement must be stopped.
	else
		restricted_triangles_to_refine = triangles_to_refine
	end
	new_triangulation = setdiff(triangles, restricted_triangles_to_refine) 
    new_parameters = []
	for T in restricted_triangles_to_refine
		midpoint1 = 0.5*(T[2]-T[1]) + T[1]
		midpoint2 = 0.5*(T[3]-T[2]) + T[2]
		midpoint3 = 0.5*(T[3]-T[1]) + T[1]
		push!(new_triangulation, [midpoint1, midpoint3, T[1]])
		push!(new_triangulation, [T[2], midpoint2, midpoint1])
		push!(new_triangulation, [midpoint2, T[3], midpoint3])
		push!(new_triangulation, [midpoint1, midpoint2, midpoint3])
		push!(new_parameters, midpoint1)
		push!(new_parameters, midpoint2)
		push!(new_parameters, midpoint3)
	end
    mesh1.triangulation = new_triangulation
	S = solve_over_params(mesh1.EP, new_parameters, checks = [], retry=true)
	for i in S
		if length(solutions(i[1]))!=degree(mesh1.EP) || nsingular(i[1])!=0
        	value_dict[i[2]] = false
		else
			value_dict[i[2]] = fibre_function(i)
        end
	end
    return refinement_termination, length(S)
end


#=
function standard_for_continuous_function_refinement - in the case of continuous fibre_function refinement, this function produces the standard for which triangles will be refined by looking at the range of fibre_function values across triangle vertices. 
    input:
        value_dict
        triangles
        plot_proportion - Determines the percentile for the continuous fibre_function refinement standard. For example, if plot proportion is 0.7, the 70th percentile of triangle value ranges is chosen as the standard for refinement, i.e., any triangle with a fibre_function range greater than the 70th percentile of all ranges will be selected for refinement.
    output:
        standard_for_refinement
=#
function standard_for_continuous_function_refinement(value_dict::Dict, triangles::Vector, plot_proportion = 0.7)
    triangle_value_ranges = Dict()
    for T in triangles
        if isa(value_dict[T[1]], Bool) || isa(value_dict[T[2]], Bool) || isa(value_dict[T[3]], Bool)
            non_error_coordinates = findall(x->isa(value_dict[x], Bool) == false, T)
            if length(non_error_coordinates) == 2 #because in the continuous case, triangles that have a single error can be refined, these triangles must also be considered when calculating the standard for refinement.
                triangle_value_ranges[T] = abs(value_dict[T[non_error_coordinates][1]] - value_dict[T[non_error_coordinates][2]])
            else
                continue
            end
        else
            min_value = min([value_dict[x] for x in T]...)
            max_value = max([value_dict[x] for x in T]...)
            triangle_value_ranges[T] = max_value - min_value
        end
    end
	triangle_ranges = []
	for i in values(triangle_value_ranges)
		push!(triangle_ranges, i)
	end
	standard_for_refinement = finite_data_set_percentile(triangle_ranges, plot_proportion*100)
    return standard_for_refinement
end

#= function finite_data_set_percentile
    input:
        data_set - some finite collection of numbers
        percentile - the desired percentile, e.g. 70
    output:
        sorted_data_set[percentile_index] - the smallest value in the data_set that is greater than at least (percentile)% of the values in the data_set
=#
function finite_data_set_percentile(data_set, percentile)
    sorted_data_set = sort(data_set)
    data_set_size = length(data_set)
    percentile_index = findfirst(x->x/data_set_size >= percentile/100, 1:data_set_size)
    return sorted_data_set[percentile_index]
end

#= function preliminary_mesh - computes solutions for 100 points over the visualization window to estimate the number of unique fibre_function values, e.g., number of real solutions, that will be in visualization.
    input:
        EP
        fibre_function
        xlims
        ylims
    output:
        xlims
        ylims
        length(unique_fibre_function_values) - the number of unique fibre_function values found in the mesh
=#
function preliminary_mesh(EP::EnumerativeProblem, fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2.0,2.0], ylims = [-2.0,2.0])
    x_values = range(xlims[1], xlims[2], 10)
    y_values = range(ylims[1], ylims[2], 10)
    preliminary_mesh_points = [[x,y] for x in x_values for y in y_values]
    S = solve_over_params(EP, preliminary_mesh_points, checks =[])
    unique_fibre_function_values = []
    for s in S
        if length(solutions(s[1]))!=degree(EP) || nsingular(s[1])!=0
            continue
        else
            current_fibre_function_value = fibre_function(s)
            if (current_fibre_function_value in unique_fibre_function_values) == false
                push!(unique_fibre_function_values, current_fibre_function_value)
            end
        end
    end
    return (xlims, ylims, length(unique_fibre_function_values))
end

#= function adjust_visualization_window - uses preliminary_mesh to adjust visualization window so that there are roughly no greater than 5 unique fibre_function values in the visualization.
    input:
        EP
        fibre_function
        xlims
        ylims
    output:
        X - newly adjusted x parameter boundary of visualization window
        Y - newly adjusted y parameter boundary of visualization window
=#
function adjust_visualization_window(EP::EnumerativeProblem, fibre_function = x->HomotopyContinuation.nreal(x[1]), xlims = [-2.0, 2.0], ylims = [-2.0, 2.0]) 
    X, Y, N = preliminary_mesh(EP, fibre_function, xlims, ylims)
    while N > 5
        L = X[2] - X[1]
        new_window_length = 2*sqrt((L^2)/N)
        length_to_remove_from_each_end_of_window = (L - new_window_length)/2
        X[1] = X[1] + length_to_remove_from_each_end_of_window
        X[2] = X[2] - length_to_remove_from_each_end_of_window
        Y[1] = Y[1] + length_to_remove_from_each_end_of_window
        Y[2] = Y[2] - length_to_remove_from_each_end_of_window
        X, Y, N = preliminary_mesh(EP, fibre_function, X, Y)
    end
    return X, Y
end

#=function partition_triangles_into_chambers - takes triangles and separates them into their respective chambers
    input:
        triangles - a vector of triangles
    output:
        chambers - a collection of vectors, each containing all of the triangles that make up a particular chamber
=#
function partition_triangles_into_chambers(triangles)
    remaining_triangles = Set(triangles)
    chambers = []
    for i in remaining_triangles
        triangle_queue = [i]
        delete!(remaining_triangles, i)
        current_chamber = []
        for T in triangle_queue
            if (T in current_chamber) == false
                push!(current_chamber, T)
            end
            adjacent_triangles = filter(x->(T[1] in x) || (T[2] in x) || (T[3] in x) , remaining_triangles)
            if length(adjacent_triangles) > 0
                for t in adjacent_triangles
                    push!(triangle_queue, t)
                    if (t in current_chamber) == false
                        push!(current_chamber, t)
                    end
                    delete!(remaining_triangles, t)
                end
            else
                continue
            end
        end
        push!(chambers, current_chamber)
    end
    return chambers
end

#=function collect_boundary_edges - given a collection of triangles corresponding to a chamber, returns the edges that lie on its boundary
    input:
        triangles
    output:
        boundary_edges - a set containing all of the edges that lie on the boundary of the chamber. Edges are represented as sets containing two vertices.
=#
function collect_boundary_edges(triangles)
    boundary_edges = []
    all_edges = []
    for T in triangles
        push!(all_edges, Set([T[1], T[2]]))
        push!(all_edges, Set([T[1], T[3]]))
        push!(all_edges, Set([T[2], T[3]]))
    end
    unique_edges = unique(all_edges)
    for E in unique_edges
        counter = 0
        for E2 in all_edges
            if E == E2
                counter += 1
            end
            if counter > 1
                break
            end
        end
        if counter == 1
            push!(boundary_edges, E)
        end
    end
    return boundary_edges
end

#= function check_if_boundary_has_interior_hole - separates boundary edges into distinct components to determine if a chamber has an interior hole
    input:
        edges - a collection of edges of a chamber
    output:
        returns true if the chamber has an interior hole and returns false if it does not.
=#
function check_if_boundary_has_interior_hole(edges)
    remaining_edges = copy(edges)
    boundary_components = []
    while length(remaining_edges) > 0
        current_boundary_component = []
        edge_queue = [rand(remaining_edges)]
        for e in edge_queue
            deleteat!(remaining_edges, findall(x->x == e, remaining_edges))
            adjacent_edges = filter(x->length(intersect(e,x)) == 1, remaining_edges)
            for E in adjacent_edges
                push!(edge_queue, E)
            end
        end
        push!(boundary_components, current_boundary_component)
    end
    if length(boundary_components) > 1
        return true
    elseif length(boundary_components) == 1
        return false
    end
end

#=function ordered_chamber_boundary - given the boundary edges of a chamber, returns an ordered list of vertices for plotting
    input:
        boundary_edges - boundary edges of a chamber as given by collect_boundary_edges
    output:
        ordered_boundary - ordered list of vertices of chamber boundary
=#
function ordered_chamber_boundary(boundary_edges) 
    remaining_edges = copy(boundary_edges)
    ordered_boundary = []
    E = rand(remaining_edges)
    adjacent_edges = filter(x->length(intersect(E, x))==1, remaining_edges)
    while length(adjacent_edges)!=2
        E = rand(remaining_edges)
        adjacent_edges = filter(x->length(intersect(E, x))==1, remaining_edges)
    end
    while length(remaining_edges) > 0
        current_component = []
        last_added_point = nothing
        while E != nothing
            deleteat!(remaining_edges, findall(x->x == E, remaining_edges))
            if length(current_component) == 0 #just added this
                adjacent_edges = filter(x->length(intersect(E,x)) == 1, remaining_edges) #just added this
            else #just added this
                adjacent_edges = filter(x->last_added_point in x, remaining_edges) #just added this
            end #just added this
            if length(adjacent_edges) > 0
                next_edge = adjacent_edges[1]
                if length(current_component) == 0
                    for v in setdiff(E, intersect(E, next_edge))
                        push!(current_component, v)
                    end
                    for v in intersect(E, next_edge)
                        push!(current_component, v)
                    end
                    for v in setdiff(next_edge, intersect(E, next_edge))
                        push!(current_component, v)
                        last_added_point = v #just added this
                    end
                else
                    for v in setdiff(next_edge, intersect(E, next_edge))
                        push!(current_component, v)
                        last_added_point = v #just added this
                    end
                end
                E = next_edge
            else
                if length(current_component) > 0
                    push!(ordered_boundary, current_component)
                end
                E = nothing
            end
        end 
        if length(remaining_edges) > 0
            E = rand(remaining_edges)
        end
    end
    if length(ordered_boundary) ==1
        return ordered_boundary[1]
    else
        main_boundary_index = findfirst(x->length(x) == max([length(x) for x in ordered_boundary]...), ordered_boundary)
        main_boundary = ordered_boundary[main_boundary_index]
        deleteat!(ordered_boundary, main_boundary_index)
        final_ordered_boundary = []
        for v in main_boundary
            push!(final_ordered_boundary, v)
        end
        while length(ordered_boundary) > 0
            b = rand(ordered_boundary)
            intersection_point = intersect(b, final_ordered_boundary)
            if length(intersection_point ) == 1
                for t in intersection_point
                    intersection_point = t
                end
            else
                continue
            end
            index_to_insert_point_into = findfirst(x->x == intersection_point, final_ordered_boundary)
            index_of_point_to_be_inserted = findfirst(x->x == intersection_point, b)
            for i in 1:length(b)
                insert!(final_ordered_boundary, index_to_insert_point_into, b[index_of_point_to_be_inserted])
                index_to_insert_point_into += 1
                if index_of_point_to_be_inserted == length(b)
                    index_of_point_to_be_inserted = 1
                else
                    index_of_point_to_be_inserted += 1
                end
            end
            deleteat!(ordered_boundary, findall(x->x == b, ordered_boundary))
        end
        return final_ordered_boundary
    end
end

#= function draw_chambers - used to create plot when fibre_function is discrete. Plots chambers as wholes when they do not contain interior holes.
    input:
        mesh1 - mesh object
        scatter - indicates whether triangulation vertices should be plotted
        label - indicates whether a legend should be added to plot
    output:
        my_plot - plot of parameter space
=#
function draw_chambers(mesh1::mesh;
                        scatter = false,
                        label = true)
    fibre_function_values = filter(x->isa(x,Bool) == false, unique(values(mesh1.value_dict)))
    sorted_fibre_function_values = sort(fibre_function_values)
    triangles_containing_errors = filter(x->isa(mesh1.value_dict[x[1]], Bool) || isa(mesh1.value_dict[x[2]], Bool) || isa(mesh1.value_dict[x[3]], Bool), mesh1.triangulation)
    my_plot = plot(xlims = mesh1.xlims, ylims = mesh1.ylims)
    for i in sorted_fibre_function_values
        triangles_to_plot = filter(x-> mesh1.value_dict[x[1]] == i && mesh1.value_dict[x[1]] == mesh1.value_dict[x[2]] && mesh1.value_dict[x[1]] == mesh1.value_dict[x[3]], mesh1.triangulation)
        for j in triangles_containing_errors
            non_error_vertex_indices = findall(x->isa(mesh1.value_dict[x], Bool)==false, j)
            if length(non_error_vertex_indices) == 2
                if mesh1.value_dict[j[non_error_vertex_indices[1]]] == i && mesh1.value_dict[j[non_error_vertex_indices[1]]] == mesh1.value_dict[j[non_error_vertex_indices[2]]]
                    push!(triangles_to_plot, j)
                else 
                    continue
                end
            elseif length(non_error_vertex_indices) == 1
                if mesh1.value_dict[j[non_error_vertex_indices[1]]] == i
                    push!(triangles_to_plot, j)
                else
                    continue
                end
            else
                continue
            end
        end
        current_chambers = partition_triangles_into_chambers(triangles_to_plot)
        has_this_value_been_plotted_yet = false
        for c in current_chambers
            chamber_boundary_edges = collect_boundary_edges(c)
            if check_if_boundary_has_interior_hole(chamber_boundary_edges) == true
                for T in c
                    if has_this_value_been_plotted_yet == false && label == true
                        draw_triangle(T, findfirst(x->x==i, sorted_fibre_function_values)/length(sorted_fibre_function_values), label = true, labelText = "$(i) real solutions")
                        has_this_value_been_plotted_yet = true
                    else
                        draw_triangle(T, findfirst(x->x==i, sorted_fibre_function_values)/length(sorted_fibre_function_values))
                    end
                end
            else
                ordered_boundary_vertices = ordered_chamber_boundary(chamber_boundary_edges)
                shape1 = Shape([(t[1], t[2]) for t in ordered_boundary_vertices])
                color_value = cgrad(:thermal, rev = false)[findfirst(x->x==i, sorted_fibre_function_values)/length(sorted_fibre_function_values)]
                if has_this_value_been_plotted_yet == false && label == true
                    plot!(shape1, fillcolor = color_value, linewidth = false, linecolor = color_value, label = "$(i) real solutions")
                    has_this_value_been_plotted_yet = true
                else
                    plot!(shape1, fillcolor = color_value, linewidth = false, linecolor = color_value, label = false)
                end
            end
        end
        if scatter == true
            color_value = cgrad(:thermal, rev = false)[findfirst(x->x==i, sorted_fibre_function_values)/length(sorted_fibre_function_values)]
            current_parameters = filter(x->mesh1.value_dict[x] == i, keys(mesh1.value_dict))
            scatter!([a[1] for a in current_parameters], [a[2] for a in current_parameters], markersize = 0.8, markershape =:rect, markerstrokewidth = 0.2, markercolor = color_value, labels = false)
        end

    end
    return my_plot
end

#=function Delaunay_triangulation! - takes a mesh and computes a Delaunay triangulation of the vertices, updating the triangulation in the inputted mesh
    input:
        mesh1 - mesh object
    output:
        nothing
=#
function Delaunay_triangulation!(mesh1::mesh)
	vertices = hcat(keys(mesh1.value_dict)...)
	tri = triangulate(vertices)
	triangle_iterator = each_solid_triangle(tri)
	triangles = []
	for T in triangle_iterator
		i,j,k = triangle_vertices(T)
		i,j,k = get_point(tri, i, j, k)
		vertex_1 = [i[1], i[2]]
		vertex_2 = [j[1], j[2]]
		vertex_3 = [k[1], k[2]]
		push!(triangles, [vertex_1,vertex_2,vertex_3])
	end
    mesh1.triangulation = triangles
	return nothing
end

#= function Delaunay_visualization
    input:
        EP 
        xlims
        ylims
        fibre_function
        initial_resolution - maximum number of vertices solved for in initial mesh
        total_resolution - maximum number of vertices solved for over the entirety of the visualizaiton process (initial mesh and subsequent refinement)
        continuous - indicates whether the fibre_function is continuous
        plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
        label - indicates whether a legend will be included in the plot
        scatter - indicates whether triangulation vertices will be plotted
        plot_logged_values - indicates whether, in the case of a continuous fibre_function, if logged fibre_function values are to be plotted
    output:
        mesh1 - mesh object resulting from refinement process
        my_plot - parameter space visualization
=#
function Delaunay_visualization(EP::EnumerativeProblem; 
                                xlims = [-2.0,2.0], 
                                ylims = [-2.0,2.0], 
                                fibre_function = x->HomotopyContinuation.nreal(x[1]), 
                                initial_resolution = 1000, 
                                total_resolution = 4*initial_resolution, 
                                continuous = false, 
                                plot_proportion = 0.7,
                                label = true,
                                scatter = false,
                                plot_logged_values = false)
    xlims, ylims = adjust_visualization_window(EP, fibre_function, xlims, ylims)
	mesh1, L = initialize_mesh(EP, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = initial_resolution, continuous = continuous)
    V = mesh1.value_dict
    T = mesh1.triangulation
	total_resolution -= L
	if continuous == true
		refinement_standard = standard_for_continuous_function_refinement(V, T, plot_proportion)
		while total_resolution > 0
			refinement_termination, L = triforce_refinement!(mesh1, total_resolution, refinement_standard)
			Delaunay_triangulation!(mesh1)
			if refinement_termination == true
				break
			end
			total_resolution -= L
		end
	else
		while total_resolution > 0
            refinement_termination, L = triforce_refinement!(mesh1, total_resolution)
			Delaunay_triangulation!(mesh1)
			if refinement_termination == true
				break
			end
			total_resolution -= L
		end
	end
    my_plot = visualize_mesh(mesh1, scatter = scatter, label = label, plot_logged_values = plot_logged_values)
	return mesh1, my_plot
end

#= function visualize_mesh
    input:
        mesh1 - mesh object
        scatter - indicates whether triangulation vertices should be plotted
        label - indicates whether the visualization should include a legend
        plot_logged_values - indicates whether, in the case of a continuous fibre_function, if logged fibre_function values are to be plotted
    output:
        my_plot - plot of parameter space
=#
function visualize_mesh(mesh1::mesh;
    scatter = false,
    label = true,
    plot_logged_values = false)
    if mesh1.type_of_fibre_function == "continuous"
        my_plot = draw_mesh(mesh1.value_dict, mesh1.triangulation, xlims = mesh1.xlims, ylims = mesh1.ylims, continuous = true, scatter = scatter, label = label, plot_logged_values = plot_logged_values)
        return my_plot
    else
        my_plot = draw_chambers(mesh1, scatter = scatter, label = label)
        return my_plot
    end
end

#= function triforce_visualization
    input:
        EP 
        xlims
        ylims
        fibre_function
        initial_resolution - maximum number of vertices solved for in initial mesh
        total_resolution - maximum number of vertices solved for over the entirety of the visualizaiton process (initial mesh and subsequent refinement)
        continuous - indicates whether the fibre_function is continuous
        plot_proportion - in the case that fibre_function is continuous, plot_proportion is used to calculate the threshold for both refinement. For continuous fibre functions, the
		range of fibre_function values across a triangle is calculated (range across the values of the three vertices) and this range is used to determine if the triangle will be refined. For example,if plot_proportion is 0.7, the 70th percentile of triangle value ranges is used as the threshold for refinement; any triangle with a value range greater than the 70th percentile will be refined. 
        label - indicates whether a legend will be included in the plot
        scatter - indicates whether triangulation vertices will be plotted
        plot_logged_values - indicates whether, in the case of a continuous fibre_function, if logged fibre_function values are to be plotted
    output:
        mesh1 - mesh object resulting from refinement process
        my_plot - parameter space visualization
=#
function triforce_visualization(EP::EnumerativeProblem;
                                    xlims = [-2.0,2.0], 
                                    ylims = [-2.0,2.0], 
                                    fibre_function = x->HomotopyContinuation.nreal(x[1]), 
                                    initial_resolution = 1000, 
                                    total_resolution = 4*initial_resolution, 
                                    continuous = false, 
                                    plot_proportion = 0.7,
                                    label = true,
                                    scatter = false,
                                    plot_logged_values = false)
    xlims, ylims = adjust_visualization_window(EP, fibre_function, xlims, ylims)
    mesh1, L = initialize_mesh(EP, fibre_function = fibre_function, xlims = xlims, ylims = ylims, resolution = initial_resolution, continuous = continuous)
    V = mesh1.value_dict
    T = mesh1.triangulation
    total_resolution -= L
    if continuous == true
        refinement_standard = standard_for_continuous_function_refinement2(V, T, plot_proportion)
        while total_resolution > 0
            refinement_termination, L = triforce_refinement!(mesh1, total_resolution, refinement_standard)
            if refinement_termination == true
                break
            end
            total_resolution -= L
        end
    else
        while total_resolution > 0
            refinement_termination, L = triforce_refinement!(mesh1, total_resolution)
            if refinement_termination == true
                break
            end
            total_resolution -= L
        end
    end
    my_plot = draw_mesh(mesh1.value_dict, mesh1.triangulation, xlims = mesh1.xlims, ylims = mesh1.ylims, continuous = continuous, scatter = scatter, label = label, plot_logged_values = plot_logged_values)

    return mesh1, my_plot
end

function visualize_parameterspace(EP::EnumerativeProblem, strategy = "Delaunay";
                    xlims = [-2.0,2.0], 
                    ylims = [-2.0,2.0], 
                    fibre_function = x->HomotopyContinuation.nreal(x[1]), 
                    initial_resolution = 1000, 
                    total_resolution = 4*initial_resolution, 
                    continuous = false, 
                    plot_proportion = 0.7,
                    label = true,
                    scatter = false,
                    plot_logged_values = false,
                    near = false)
    if isa(near, Bool) == false
        EP = restrict_enumerative_problem(EP, [near + randn(Float64, 20) for i in 1:3])
    end
    if strategy == "Delaunay"
        mesh1, my_plot = Delaunay_visualization(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, initial_resolution = initial_resolution, total_resolution = total_resolution, continuous = continuous, plot_proportion = plot_proportion, label = label, scatter = scatter, plot_logged_values = plot_logged_values)
        return mesh1, my_plot
    elseif strategy == "triforce"
        mesh1, my_plot = triforce_visualization(EP, xlims = xlims, ylims = ylims, fibre_function = fibre_function, initial_resolution = initial_resolution, total_resolution = total_resolution, continuous = continuous, plot_proportion = plot_proportion, label = label, scatter = scatter, plot_logged_values = plot_logged_values)
        return mesh1, my_plot
    else
        throw(ArgumentError("Invalid strategy inputted. Valid strategies include: Delaunay, triforce"))
    end
end

